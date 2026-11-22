#include "PathCut.hpp"

#include "AABBMesh.hpp"
#include "I18N.hpp"
#include "TriangleMesh.hpp"

#include <earcut/earcut.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Slic3r { namespace PathCut {

namespace {

Frame frame_from_normal(const V3 &normal, const V3 &center)
{
    Frame f;
    f.c       = center;
    f.n       = normal.normalized();
    V3 a      = std::abs(f.n.x()) < 0.9 ? V3::UnitX() : V3::UnitY();
    f.u       = f.n.cross(a).normalized();
    f.v       = f.n.cross(f.u);
    return f;
}

// Uniform Catmull-Rom spline through all control points.
template<class P> std::vector<P> catmull_rom(const std::vector<P> &pts, bool closed, int samples_per_segment)
{
    const int n = int(pts.size());
    if (n < 3 || samples_per_segment < 2)
        return pts;
    auto at = [&](int i) -> P {
        if (closed)
            return pts[((i % n) + n) % n];
        if (i < 0) return pts[0] * 2. - pts[1];
        if (i >= n) return pts[n - 1] * 2. - pts[n - 2];
        return pts[i];
    };
    const int      segs = closed ? n : n - 1;
    std::vector<P> out;
    out.reserve(size_t(segs * samples_per_segment + 1));
    for (int s = 0; s < segs; ++s) {
        const P p0 = at(s - 1), p1 = at(s), p2 = at(s + 1), p3 = at(s + 2);
        for (int k = 0; k < samples_per_segment; ++k) {
            const double t = double(k) / samples_per_segment, t2 = t * t, t3 = t2 * t;
            out.push_back(((p1 * 2.) + (p2 - p0) * t + (p0 * 2. - p1 * 5. + p2 * 4. - p3) * t2 + (p1 * 3. - p0 - p2 * 3. + p3) * t3) * 0.5);
        }
    }
    if (!closed)
        out.push_back(pts.back());
    return out;
}

template<class P> std::vector<P> remove_close_points(const std::vector<P> &pts, double eps, bool closed)
{
    std::vector<P> out;
    for (const P &p : pts)
        if (out.empty() || (p - out.back()).norm() > eps)
            out.push_back(p);
    if (closed)
        while (out.size() > 2 && (out.front() - out.back()).norm() <= eps)
            out.pop_back();
    return out;
}


// Point-in-polygon (even-odd) for a closed 2D polygon.
bool point_in_polygon(const std::vector<V2> &poly, const V2 &p)
{
    bool in = false;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const V2 &a = poly[i], &b = poly[j];
        if ((a.y() > p.y()) != (b.y() > p.y()) && p.x() < (b.x() - a.x()) * (p.y() - a.y()) / (b.y() - a.y()) + a.x())
            in = !in;
    }
    return in;
}

double polygon_area(const std::vector<V2> &poly)
{
    double a = 0.;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
        a += poly[j].x() * poly[i].y() - poly[i].x() * poly[j].y();
    return 0.5 * a;
}

// Parameter of the point of segment ab closest to q, and its distance to q.
double closest_on_segment(const V2 &q, const V2 &a, const V2 &b, double &dist)
{
    const V2     ab = b - a;
    const double l2 = ab.squaredNorm();
    const double t  = l2 > 0. ? std::clamp((q - a).dot(ab) / l2, 0., 1.) : 0.;
    dist            = (q - (a + ab * t)).norm();
    return t;
}


uint64_t edge_key(int a, int b) { return a < b ? (uint64_t(uint32_t(a)) << 32) | uint32_t(b) : (uint64_t(uint32_t(b)) << 32) | uint32_t(a); }

// Splits every triangle that has marked edges (mid(a, b) >= 0) into 2, 3 or 4 triangles, keeping
// the winding. Marks are per edge, so neighbouring triangles stay conforming.
template<class MidOf> void split_marked(std::vector<stl_triangle_vertex_indices> &tris, MidOf mid, const std::vector<size_t> *only = nullptr)
{
    const size_t n = only ? only->size() : tris.size();
    for (size_t ii = 0; ii < n; ++ii) {
        const size_t i = only ? (*only)[ii] : ii;
        int v[3] = {tris[i][0], tris[i][1], tris[i][2]};
        int m[3] = {mid(v[0], v[1]), mid(v[1], v[2]), mid(v[2], v[0])};
        const int cnt = int(m[0] >= 0) + int(m[1] >= 0) + int(m[2] >= 0);
        if (cnt == 0)
            continue;
        // Rotate to a canonical pattern: 1 mark on edge 0, 2 marks on edges 0 and 1.
        int r = 0;
        if (cnt == 1) r = m[0] >= 0 ? 0 : m[1] >= 0 ? 1 : 2;
        if (cnt == 2) r = m[2] < 0 ? 0 : m[0] < 0 ? 1 : 2;
        const int a = v[r], b = v[(r + 1) % 3], c = v[(r + 2) % 3];
        const int mab = m[r], mbc = m[(r + 1) % 3], mca = m[(r + 2) % 3];
        if (cnt == 1) {
            tris[i] = {a, mab, c};
            tris.emplace_back(mab, b, c);
        } else if (cnt == 2) {
            tris[i] = {mab, b, mbc};
            tris.emplace_back(a, mab, mbc);
            tris.emplace_back(a, mbc, c);
        } else {
            tris[i] = {a, mab, mca};
            tris.emplace_back(mab, b, mbc);
            tris.emplace_back(mca, mbc, c);
            tris.emplace_back(mab, mbc, mca);
        }
    }
}

// Lawson edge flips: turns a triangulation of a 2D polygon (earcut output, full of slivers) into the
// constrained Delaunay triangulation of the same points. Edges in `fixed` are never flipped.
// `sign` is the orientation of the triangles in the chart (+1 counter-clockwise).
template<class UV, class Fixed>
void delaunay_flips(std::vector<stl_triangle_vertex_indices> &tris, const UV &uv, const Fixed &fixed, double sign)
{
    auto key = [](int a, int b) { return a < b ? (uint64_t(uint32_t(a)) << 32) | uint32_t(b) : (uint64_t(uint32_t(b)) << 32) | uint32_t(a); };
    std::unordered_map<uint64_t, std::array<int, 2>> adj; // edge -> up to two triangles
    adj.reserve(tris.size() * 2);
    auto link = [&](int t) {
        for (int e = 0; e < 3; ++e) {
            auto [it, ins] = adj.emplace(key(tris[t][e], tris[t][(e + 1) % 3]), std::array<int, 2>{t, -1});
            if (!ins && it->second[0] != t && it->second[1] != t)
                (it->second[0] < 0 ? it->second[0] : it->second[1]) = t;
        }
    };
    auto unlink = [&](int t) {
        for (int e = 0; e < 3; ++e) {
            auto it = adj.find(key(tris[t][e], tris[t][(e + 1) % 3]));
            if (it == adj.end()) continue;
            if (it->second[0] == t) it->second[0] = -1;
            if (it->second[1] == t) it->second[1] = -1;
        }
    };
    for (int t = 0; t < int(tris.size()); ++t) link(t);
    auto orient = [&](int a, int b, int c) {
        const auto &A = uv.at(a), &B = uv.at(b), &C = uv.at(c);
        return sign * ((B.x() - A.x()) * (C.y() - A.y()) - (B.y() - A.y()) * (C.x() - A.x()));
    };
    auto in_circle = [&](int a, int b, int c, int d) { // > 0: d inside the circumcircle of (a, b, c)
        const auto &A = uv.at(a), &B = uv.at(b), &C = uv.at(c), &D = uv.at(d);
        const double adx = A.x() - D.x(), ady = A.y() - D.y(), bdx = B.x() - D.x(), bdy = B.y() - D.y(), cdx = C.x() - D.x(), cdy = C.y() - D.y();
        const double det = (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) - (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) +
                           (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
        return sign * det;
    };
    std::vector<uint64_t> stack;
    stack.reserve(adj.size());
    for (const auto &[k, ts] : adj) if (ts[0] >= 0 && ts[1] >= 0 && !fixed.count(k)) stack.push_back(k);
    size_t guard = 0, limit = 50 * tris.size() + 1000;
    while (!stack.empty() && ++guard < limit) {
        const uint64_t k = stack.back();
        stack.pop_back();
        auto it = adj.find(k);
        if (it == adj.end() || it->second[0] < 0 || it->second[1] < 0 || fixed.count(k)) continue;
        const int t1 = it->second[0], t2 = it->second[1];
        // t1 = (a, b, c) with the shared edge a->b, t2 = (b, a, d).
        int e1 = 0;
        while (e1 < 3 && key(tris[t1][e1], tris[t1][(e1 + 1) % 3]) != k) ++e1;
        const int a = tris[t1][e1], b = tris[t1][(e1 + 1) % 3], c = tris[t1][(e1 + 2) % 3];
        int d = -1;
        for (int e = 0; e < 3; ++e) if (tris[t2][e] != a && tris[t2][e] != b) d = tris[t2][e];
        if (d < 0 || c == d || in_circle(a, b, c, d) <= 1e-12) continue;
        if (orient(a, d, c) <= 0. || orient(b, c, d) <= 0.) continue; // quad not convex
        unlink(t1); unlink(t2);
        tris[t1] = {a, d, c};
        tris[t2] = {b, c, d};
        link(t1); link(t2);
        for (uint64_t kk : {key(a, d), key(d, b), key(b, c), key(c, a)}) stack.push_back(kk);
    }
}




// Closes the openings left by a cut: chains the directed cut edges of each level (+1, -1) into loops,
// triangulates them in the cutter's chart (outer loops with their holes), lifts the caps onto the level
// set and refines them so they follow curved walls. Cut edges that receive cap vertices split the
// triangle owning them (cut_owner: cut edge -> index in its.indices), keeping the mesh watertight.
// Returns 0 when done, 1 on an error (in res.error), 2 when canceled.
int close_cut_faces(indexed_triangle_set &its, std::vector<std::pair<int, int>> (&cap_edges)[2], const Cutter &cutter,
                    const std::unordered_map<uint64_t, size_t> &cut_owner, CutResult &res, const std::function<bool()> &canceled)
{
    auto stop = [&]() { return canceled && canceled(); };
    auto pos  = [&](int v) -> V3 { return its.vertices[size_t(v)].template cast<double>(); };
    // 3) Close the openings: chain the cut edges into loops, triangulate them in the chart of their
    //    level set and refine the caps so they follow curved walls.
    const double sag_tol = std::max(0.1, 0.25 * cutter.kerf); // 0.1 mm: half a layer; finer only multiplies the cap triangles
    const double min_len2 = 0.25 * 0.25; // never split cap edges shorter than this (chart mm)
    std::unordered_map<uint64_t, std::vector<int>> inserts; // cut edge -> cap vertices placed on it
    std::unordered_map<uint64_t, uint64_t>         orig_of; // piece of a cut edge -> the cut edge
    for (int li = 0; li < 2; ++li) {
        const int level = li == 0 ? 1 : -1;
        auto     &edges = cap_edges[li];
        if (edges.empty()) continue;
        std::unordered_map<int, std::vector<size_t>> outgoing;
        for (size_t i = 0; i < edges.size(); ++i) outgoing[edges[i].first].push_back(i);
        std::vector<bool>             used(edges.size(), false);
        std::vector<std::vector<int>> loops;
        for (size_t s = 0; s < edges.size(); ++s) {
            if (used[s]) continue;
            std::vector<int> loop;
            size_t           e = s;
            while (true) {
                used[e] = true;
                loop.push_back(edges[e].first);
                const int to = edges[e].second;
                if (to == edges[s].first) break;
                size_t next = size_t(-1);
                for (size_t c : outgoing[to]) if (!used[c]) { next = c; break; }
                if (next == size_t(-1)) {
                    res.error = _u8L("The model has holes along the cut, so the cut faces cannot be closed. Repair the model first.");
                    return 1;
                }
                e = next;
            }
            if (loop.size() >= 3) loops.push_back(std::move(loop));
        }
        // Chart coordinates of the loop vertices; outer loops vs holes by orientation.
        std::unordered_map<int, V2>          uv;
        std::vector<std::vector<V2>>         rings(loops.size());
        std::vector<double>                  area(loops.size());
        size_t                               biggest = 0;
        for (size_t i = 0; i < loops.size(); ++i) {
            for (int v : loops[i]) {
                const V2 c = cutter.chart(pos(v));
                uv[v]      = c;
                rings[i].push_back(c);
            }
            area[i] = polygon_area(rings[i]);
            if (std::abs(area[i]) > std::abs(area[biggest])) biggest = i;
        }
        if (loops.empty()) continue;
        const double outer_sign = area[biggest] > 0. ? 1. : -1.;
        std::vector<int> parent(loops.size(), -1);
        for (size_t i = 0; i < loops.size(); ++i) {
            if (area[i] * outer_sign > 0.) continue;
            double best = std::numeric_limits<double>::max();
            for (size_t j = 0; j < loops.size(); ++j)
                if (area[j] * outer_sign > 0. && std::abs(area[j]) < best && point_in_polygon(rings[j], rings[i].front())) {
                    best      = std::abs(area[j]);
                    parent[i] = int(j);
                }
        }
        for (size_t o = 0; o < loops.size(); ++o) {
            const bool is_outer = area[o] * outer_sign > 0. || parent[o] < 0; // orphan holes are closed on their own
            if (!is_outer) continue;
            std::vector<std::vector<std::array<double, 2>>> polygon;
            std::vector<int>                                ids;
            auto add_ring = [&](size_t r) {
                polygon.emplace_back();
                for (size_t k = 0; k < loops[r].size(); ++k) {
                    polygon.back().push_back({rings[r][k].x(), rings[r][k].y()});
                    ids.push_back(loops[r][k]);
                }
            };
            add_ring(o);
            for (size_t h = 0; h < loops.size(); ++h)
                if (parent[h] == int(o)) add_ring(h);
            const std::vector<uint32_t> tri = mapbox::earcut<uint32_t>(polygon);
            const double sign = area[o] > 0. ? 1. : -1.;
            std::unordered_set<uint64_t> ring_dir; // ring edges, in the direction the cap must walk them
            for (size_t r = 0; r < loops.size(); ++r)
                if (int(r) == int(o) || parent[r] == int(o))
                    for (size_t k = 0; k < loops[r].size(); ++k)
                        ring_dir.insert((uint64_t(uint32_t(loops[r][k])) << 32) | uint32_t(loops[r][(k + 1) % loops[r].size()]));
            std::vector<stl_triangle_vertex_indices> cap;
            for (size_t k = 0; k + 2 < tri.size(); k += 3) {
                int a = ids[tri[k]], b = ids[tri[k + 1]], c = ids[tri[k + 2]];
                const V2    &pa = uv[a], &pb = uv[b], &pc = uv[c];
                const double cr = ((pb - pa).x() * (pc - pa).y() - (pb - pa).y() * (pc - pa).x()) * sign;
                // Earcut emits zero-area triangles along straight runs of the ring, sometimes walking
                // ring edges backwards; leave those out, the gap closing below covers the spot.
                if (cr == 0.) continue;
                if (cr < 0.) std::swap(b, c);
                const int tv[3] = {a, b, c};
                bool      backwards = false;
                for (int e = 0; e < 3; ++e)
                    backwards |= ring_dir.count((uint64_t(uint32_t(tv[(e + 1) % 3])) << 32) | uint32_t(tv[e])) > 0;
                if (backwards) continue;
                cap.emplace_back(a, b, c);
            }
            // The triangulation is flat in the chart but the level set may be curved in 3D (a bent wall,
            // a domed cap): split interior edges whose midpoint is too far from the level set, placing
            // the new vertex on it. Boundary edges are shared with the clipped surface, never split.
            std::unordered_set<uint64_t> boundary;
            for (size_t r = 0; r < loops.size(); ++r)
                if (int(r) == int(o) || parent[r] == int(o))
                    for (size_t k = 0; k < loops[r].size(); ++k)
                        boundary.insert(edge_key(loops[r][k], loops[r][(k + 1) % loops[r].size()]));
            // Earcut skips pieces of rings that touch or cross themselves in the chart (it happens at
            // micrometre scale where the cut grazes the surface). Close whatever it left open: chain the
            // ring edges no triangle uses with the reversed open edges of the triangulation and fan them.
            {
                auto dkey = [](int a, int b) { return (uint64_t(uint32_t(a)) << 32) | uint32_t(b); };
                std::unordered_set<uint64_t> used, ring;
                for (const auto &t : cap) for (int e = 0; e < 3; ++e) used.insert(dkey(t[e], t[(e + 1) % 3]));
                std::vector<std::pair<int, int>> gap;
                for (size_t r = 0; r < loops.size(); ++r)
                    if (int(r) == int(o) || parent[r] == int(o))
                        for (size_t k = 0; k < loops[r].size(); ++k) {
                            const int a = loops[r][k], b = loops[r][(k + 1) % loops[r].size()];
                            ring.insert(dkey(a, b));
                            if (!used.count(dkey(a, b))) gap.emplace_back(a, b);
                        }
                for (const auto &t : cap)
                    for (int e = 0; e < 3; ++e) {
                        const int a = t[e], b = t[(e + 1) % 3];
                        if (!used.count(dkey(b, a)) && !ring.count(dkey(a, b))) gap.emplace_back(b, a);
                    }
                std::unordered_map<int, std::vector<size_t>> from;
                for (size_t i = 0; i < gap.size(); ++i) from[gap[i].first].push_back(i);
                std::vector<bool> done(gap.size(), false);
                // A chain can pass a vertex twice (rings pinched in the chart): split it there into simple
                // cycles, and fan each from a new centre vertex so no fan edge can duplicate an existing one.
                auto close_cycle = [&](const int *c, size_t n) {
                    if (n < 3) return;
                    if (n == 3) { cap.emplace_back(c[0], c[1], c[2]); return; }
                    V2 m = V2::Zero();
                    for (size_t k = 0; k < n; ++k) m += uv[c[k]];
                    m /= double(n);
                    const int x = int(its.vertices.size());
                    uv[x]       = m;
                    its.vertices.emplace_back(cutter.lift(m, level).template cast<float>());
                    for (size_t k = 0; k < n; ++k) cap.emplace_back(x, c[k], c[(k + 1) % n]);
                };
                for (size_t s0 = 0; s0 < gap.size(); ++s0) {
                    if (done[s0]) continue;
                    std::vector<int>                cyc;
                    std::unordered_map<int, size_t> at;
                    size_t e = s0;
                    for (size_t guard = 0; guard <= gap.size(); ++guard) {
                        done[e] = true;
                        at[gap[e].first] = cyc.size();
                        cyc.push_back(gap[e].first);
                        const int w = gap[e].second;
                        if (auto it = at.find(w); it != at.end()) {
                            const size_t p = it->second;
                            close_cycle(cyc.data() + p, cyc.size() - p);
                            for (size_t k = p; k < cyc.size(); ++k) at.erase(cyc[k]);
                            cyc.resize(p);
                        }
                        size_t next = size_t(-1);
                        for (size_t c : from[w]) if (!done[c]) { next = c; break; }
                        if (next == size_t(-1)) break;
                        e = next;
                    }
                }
            }
            // Lifted chart position of each cap vertex (the boundary vertices themselves sit on the
            // linearly interpolated level set, which is not what the flatness test should measure).
            std::unordered_map<int, V3> lifted;
            auto lifted_at = [&](int v) -> const V3 & {
                auto it = lifted.find(v);
                if (it == lifted.end()) it = lifted.emplace(v, cutter.lift(uv[v], level)).first;
                return it->second;
            };
            delaunay_flips(cap, uv, boundary, sign);
            // Slice the cap along the kink lines of the level set, so every piece lies where the level set
            // is planar (or nearly). A slice through a boundary edge also splits the clipped triangle
            // that owns it (see below), keeping the mesh watertight.
            for (const auto &kink : cutter.kinks(sag_tol)) {
                const V2 p0 = kink.first, nrm(-kink.second.y(), kink.second.x());
                auto dist = [&](int v) { const double d = (uv[v] - p0).dot(nrm); return std::abs(d) < 1e-5 ? 0. : d; };
                std::unordered_map<uint64_t, int> xing;
                auto cross_at = [&](int a, int b) {
                    const uint64_t k = edge_key(a, b);
                    auto [it, ins]   = xing.emplace(k, -1);
                    if (!ins) return it->second;
                    const int    lo = std::min(a, b), hi = std::max(a, b);
                    const double t  = dist(lo) / (dist(lo) - dist(hi));
                    const int    x  = int(its.vertices.size());
                    const V2     c  = uv[lo] + (uv[hi] - uv[lo]) * t;
                    uv[x]           = c;
                    if (boundary.erase(k)) {
                        // On the boundary the vertex must stay on the straight cut edge.
                        its.vertices.emplace_back((pos(lo) + (pos(hi) - pos(lo)) * t).template cast<float>());
                        auto      o    = orig_of.find(k);
                        const uint64_t orig = o == orig_of.end() ? k : o->second;
                        inserts[orig].push_back(x);
                        boundary.insert(edge_key(lo, x)); orig_of[edge_key(lo, x)] = orig;
                        boundary.insert(edge_key(x, hi)); orig_of[edge_key(x, hi)] = orig;
                    } else {
                        its.vertices.emplace_back(cutter.lift(c, level).template cast<float>());
                    }
                    it->second = x;
                    return x;
                };
                const size_t n = cap.size();
                for (size_t i = 0; i < n; ++i) {
                    const auto t = cap[i];
                    const double d[3] = {dist(t[0]), dist(t[1]), dist(t[2])};
                    if (!((d[0] > 0. || d[1] > 0. || d[2] > 0.) && (d[0] < 0. || d[1] < 0. || d[2] < 0.))) continue;
                    int P[4], np = 0, N[4], nn = 0;
                    for (int k = 0; k < 3; ++k) {
                        const double da = d[k], db = d[(k + 1) % 3];
                        if (da >= 0.) P[np++] = t[k];
                        if (da <= 0.) N[nn++] = t[k];
                        if ((da > 0. && db < 0.) || (da < 0. && db > 0.)) {
                            const int x = cross_at(t[k], t[(k + 1) % 3]);
                            P[np++] = x; N[nn++] = x;
                        }
                    }
                    cap[i] = {P[0], P[1], P[2]};
                    for (int k = 2; k + 1 < np; ++k) cap.emplace_back(P[0], P[k], P[k + 1]);
                    for (int k = 1; k + 1 < nn; ++k) cap.emplace_back(N[0], N[k], N[k + 1]);
                }
            }
            // Longest interior edge of a cap triangle (in the chart), -1 if all its edges are boundary.
            auto longest = [&](const stl_triangle_vertex_indices &t) {
                int best = -1; double bl = -1.;
                for (int e = 0; e < 3; ++e) {
                    if (boundary.count(edge_key(t[e], t[(e + 1) % 3]))) continue;
                    const double l = (uv[t[e]] - uv[t[(e + 1) % 3]]).squaredNorm();
                    if (l > bl) { bl = l; best = e; }
                }
                return best;
            };
            for (int pass = 0; pass < 60 && !stop(); ++pass) {
                // Triangles with a curved edge get their longest edge split; marks then propagate to the
                // longest edge of every triangle sharing a marked edge, so the splits keep their shape.
                std::unordered_map<uint64_t, int> mid;
                for (const auto &t : cap) {
                    bool curved = false;
                    for (int e = 0; e < 3 && !curved; ++e) {
                        const int a = t[e], b = t[(e + 1) % 3];
                        if (boundary.count(edge_key(a, b)) || (uv[a] - uv[b]).squaredNorm() < min_len2) continue;
                        curved = (cutter.lift((uv[a] + uv[b]) * 0.5, level) - (lifted_at(a) + lifted_at(b)) * 0.5).norm() > sag_tol;
                    }
                    const int e = curved ? longest(t) : -1;
                    if (e >= 0) mid.emplace(edge_key(t[e], t[(e + 1) % 3]), -1);
                }
                if (mid.empty()) break;
                for (bool changed = true; changed && !stop();) {
                    changed = false;
                    for (const auto &t : cap) {
                        const bool marked = mid.count(edge_key(t[0], t[1])) || mid.count(edge_key(t[1], t[2])) || mid.count(edge_key(t[2], t[0]));
                        if (!marked) continue;
                        const int e = longest(t);
                        if (e >= 0 && mid.emplace(edge_key(t[e], t[(e + 1) % 3]), -1).second) changed = true;
                    }
                }
                for (auto &[k, idx] : mid) {
                    const int a = int(k >> 32), b = int(k & 0xffffffffu);
                    const V2  c = (uv[a] + uv[b]) * 0.5;
                    idx         = int(its.vertices.size());
                    uv[idx]     = c;
                    its.vertices.emplace_back(lifted_at(idx).template cast<float>());
                }
                split_marked(cap, [&mid](int a, int b) { auto it = mid.find(edge_key(a, b)); return it == mid.end() ? -1 : it->second; });
                if (cap.size() > 2000000) break;
            }
            res.level.resize(its.vertices.size(), 0);
            for (const auto &t : cap)
                for (int k = 0; k < 3; ++k) res.level[size_t(t[k])] = static_cast<signed char>(level);
            its.indices.insert(its.indices.end(), cap.begin(), cap.end());
        }
        if (stop()) return 2;
    }
    // Split the clipped triangles whose cut edge received cap vertices.
    for (auto &[k, xs] : inserts) {
        const size_t ti = cut_owner.at(k);
        const auto   t  = its.indices[ti];
        int e = 0;
        while (e < 3 && edge_key(t[e], t[(e + 1) % 3]) != k) ++e;
        if (e == 3) continue;
        const int a = t[e], b = t[(e + 1) % 3], o = t[(e + 2) % 3];
        std::sort(xs.begin(), xs.end(), [&](int x, int y) { return (pos(x) - pos(a)).squaredNorm() < (pos(y) - pos(a)).squaredNorm(); });
        std::vector<int> chain{a};
        chain.insert(chain.end(), xs.begin(), xs.end());
        chain.push_back(b);
        its.indices[ti] = {chain[0], chain[1], o};
        for (size_t j = 1; j + 1 < chain.size(); ++j) its.indices.emplace_back(chain[j], chain[j + 1], o);
    }
    return 0;
}



// A limited view wall only cuts what it goes all the way through. Where the wall ends inside the model
// (the leg next to a sword, reached by the limits), the band of triangles it crosses there touches the
// wall's rim: that band is left uncut, and the wall cuts the parts it crosses completely. When every band
// touches the rim nothing is skipped, so the cut reports that it ends inside the object.
std::vector<char> foreign_bands(const indexed_triangle_set &its, const std::vector<double> &f, const std::vector<double> &out, double hk)
{
    const size_t      nt = its.indices.size();
    std::vector<char> skip(nt, 0);
    auto crosses = [&](const stl_triangle_vertex_indices &t) {
        const double f0 = f[t[0]], f1 = f[t[1]], f2 = f[t[2]];
        return !(std::min({f0, f1, f2}) >= hk || std::max({f0, f1, f2}) <= -hk);
    };
    auto inside = [&](const stl_triangle_vertex_indices &t) { return out[t[0]] < 0. && out[t[1]] < 0. && out[t[2]] < 0.; };
    std::unordered_map<uint64_t, std::vector<uint32_t>> by_edge; // edges of the bands
    std::unordered_set<uint64_t>                         rim;     // edges of triangles crossing the wall's rim
    std::vector<int>                                     comp(nt, -1);
    for (size_t ti = 0; ti < nt; ++ti) {
        const auto &t = its.indices[ti];
        if (!crosses(t)) continue;
        const bool in = inside(t);
        if (in) comp[ti] = -2;
        for (int e = 0; e < 3; ++e) {
            const uint64_t k = edge_key(t[e], t[(e + 1) % 3]);
            if (in) by_edge[k].push_back(uint32_t(ti));
            else rim.insert(k);
        }
    }
    std::vector<char>   touches_rim; // per band
    std::vector<size_t> stack;
    for (size_t s = 0; s < nt; ++s) {
        if (comp[s] != -2) continue;
        const int c = int(touches_rim.size());
        bool      r = false;
        comp[s]     = c;
        stack.assign(1, s);
        while (!stack.empty()) {
            const size_t ti = stack.back(); stack.pop_back();
            const auto  &t  = its.indices[ti];
            for (int e = 0; e < 3; ++e) {
                const uint64_t k = edge_key(t[e], t[(e + 1) % 3]);
                r = r || rim.count(k) > 0;
                for (uint32_t n : by_edge[k])
                    if (comp[n] == -2) { comp[n] = c; stack.push_back(n); }
            }
        }
        touches_rim.push_back(r ? 1 : 0);
    }
    if (std::all_of(touches_rim.begin(), touches_rim.end(), [](char r) { return r != 0; }))
        return skip;
    for (size_t ti = 0; ti < nt; ++ti)
        if (comp[ti] >= 0 && touches_rim[size_t(comp[ti])]) skip[ti] = 1;
    return skip;
}

} // namespace

V2 ViewCutter::seg_normal(size_t i) const
{
    const V2 d = (path[i + 1] - path[i]).normalized();
    return {-d.y(), d.x()};
}

// Offset direction at a path vertex: the bisector of its segments, so the wall faces have no gaps or
// overlaps at the bends.
V2 ViewCutter::vertex_normal(size_t k) const
{
    if (k == 0) return seg_normal(0);
    if (k + 1 >= path.size()) return seg_normal(path.size() - 2);
    return (seg_normal(k - 1) + seg_normal(k)).normalized();
}

// Sign of the side of the nearest feature (segment normal, or the mean normal at a shared vertex).
double ViewCutter::signed_from(const V2 &q, size_t seg, double t, double dist) const
{
    V2 n = seg_normal(seg);
    if (t <= 0. && seg > 0) n += seg_normal(seg - 1);
    if (t >= 1. && seg + 2 < path.size()) n += seg_normal(seg + 1);
    const V2 foot = path[seg] + (path[seg + 1] - path[seg]) * t;
    return (q - foot).dot(n) >= 0. ? dist : -dist;
}

double ViewCutter::exact_field(const V2 &q) const
{
    size_t s = 0;
    double t = 0., d = 0.;
    nearest(q, s, t, d);
    return signed_from(q, s, t, d);
}

// Nearest point of the path to a point of the view plane: segment index and parameter.
void ViewCutter::nearest(const V2 &q, size_t &seg, double &t, double &dist) const
{
    dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        double       d;
        const double tt = closest_on_segment(q, path[i], path[i + 1], d);
        if (d < dist) { dist = d; seg = i; t = tt; }
    }
}

// Far from the path the field is clamped to the cell size, which is all the cut needs there.
void ViewCutter::build_grid(const V2 &mn, const V2 &mx)
{
    cell     = std::max(0.1, std::max(mx.x() - mn.x(), mx.y() - mn.y()) / 256.);
    reach    = std::max(5., 4. * cell);
    grid_min = mn - V2(cell, cell);
    nx       = int((mx.x() - mn.x()) / cell) + 3;
    ny       = int((mx.y() - mn.y()) / cell) + 3;
    cell_segs.assign(size_t(nx) * ny, {});
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        const V2 lo = path[i].cwiseMin(path[i + 1]) - V2(reach + cell, reach + cell), hi = path[i].cwiseMax(path[i + 1]) + V2(reach + cell, reach + cell);
        const int x0 = std::max(0, int((lo.x() - grid_min.x()) / cell)), x1 = std::min(nx - 1, int((hi.x() - grid_min.x()) / cell));
        const int y0 = std::max(0, int((lo.y() - grid_min.y()) / cell)), y1 = std::min(ny - 1, int((hi.y() - grid_min.y()) / cell));
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) {
                // keep only cells the (expanded) segment really passes near
                double d;
                closest_on_segment(grid_min + V2((x + 0.5) * cell, (y + 0.5) * cell), path[i], path[i + 1], d);
                if (d <= reach + 0.75 * cell)
                    cell_segs[size_t(y) * nx + x].push_back(int(i));
            }
    }
    cell_side.assign(size_t(nx) * ny, 0);
    // One exact (nearest segment) evaluation per empty cell: rows in parallel, each cell written once.
    tbb::parallel_for(tbb::blocked_range<int>(0, ny), [this](const tbb::blocked_range<int> &rows) {
        for (int y = rows.begin(); y < rows.end(); ++y)
            for (int x = 0; x < nx; ++x)
                if (cell_segs[size_t(y) * nx + x].empty())
                    cell_side[size_t(y) * nx + x] = exact_field(grid_min + V2((x + 0.5) * cell, (y + 0.5) * cell)) >= 0. ? 1 : -1;
    });
}

// Nearest path segment among those listed in the grid cell of q; false when none is near enough to be
// sure no segment outside the cell is nearer.
bool ViewCutter::nearest_in_grid(const V2 &q, size_t &seg, double &t, double &dist) const
{
    const int x = int(std::floor((q.x() - grid_min.x()) / cell)), y = int(std::floor((q.y() - grid_min.y()) / cell));
    if (x < 0 || y < 0 || x >= nx || y >= ny) return false;
    dist = std::numeric_limits<double>::max();
    for (int i : cell_segs[size_t(y) * nx + x]) {
        double       dd;
        const double tt = closest_on_segment(q, path[i], path[i + 1], dd);
        if (dd < dist) { dist = dd; seg = size_t(i); t = tt; }
    }
    return dist <= reach;
}

double ViewCutter::field(const V3 &p) const
{
    const V2  q = frame.to2d(p);
    const int x = int(std::floor((q.x() - grid_min.x()) / cell)), y = int(std::floor((q.y() - grid_min.y()) / cell));
    if (x >= 0 && y >= 0 && x < nx && y < ny && cell_segs[size_t(y) * nx + x].empty())
        return cell_side[size_t(y) * nx + x] * reach;
    size_t seg = 0;
    double t = 0., d = 0.;
    return nearest_in_grid(q, seg, t, d) ? signed_from(q, seg, t, d) : exact_field(q);
}

double ViewCutter::outside(const V3 &p) const
{
    double o = -1.;
    if (limit_depth) {
        const double h = frame.height(p);
        o = std::max(o, std::max(h_min - h, h - h_max));
    }
    if (limit_start || limit_end) {
        const V2     q  = frame.to2d(p);
        const double d0 = limit_start ? (q - path.front()).dot((path[0] - path[1]).normalized()) : -1.;
        const double d1 = limit_end ? (q - path.back()).dot((path.back() - path[path.size() - 2]).normalized()) : -1.;
        if (d0 <= 0. && d1 <= 0.)
            return o; // not beyond a limited end
        size_t seg = 0; double t = 0., d = 0.;
        // Far from the path the grid knows nothing: those points need the exact nearest segment, as the
        // field's sign still flips along the lines through the ends.
        if (!nearest_in_grid(q, seg, t, d))
            nearest(q, seg, t, d);
        if (limit_start && seg == 0 && t <= 0.) o = std::max(o, d0);
        if (limit_end && seg + 2 == path.size() && t >= 1.) o = std::max(o, d1);
    }
    return o;
}

// Only the bends that matter: Douglas-Peucker on the path with tolerance `tol`.
std::vector<std::pair<V2, V2>> ViewCutter::kinks(double tol) const
{
    std::vector<char> keep(path.size(), 0);
    keep.front() = keep.back() = 1;
    std::vector<std::pair<size_t, size_t>> stack{{0, path.size() - 1}};
    while (!stack.empty()) {
        const auto [i, j] = stack.back();
        stack.pop_back();
        const V2 d = (path[j] - path[i]).normalized();
        double   worst = tol; size_t w = 0;
        for (size_t k = i + 1; k < j; ++k) {
            const V2     r  = path[k] - path[i];
            const double dd = std::abs(r.x() * d.y() - r.y() * d.x());
            if (dd > worst) { worst = dd; w = k; }
        }
        if (w) { keep[w] = 1; stack.emplace_back(i, w); stack.emplace_back(w, j); }
    }
    std::vector<std::pair<V2, V2>> out;
    double last = -1e30;
    for (size_t i = 1; i + 1 < path.size(); ++i)
        if (keep[i] && arc[i] - last > 0.05) { // bends closer than this are one bend
            out.emplace_back(V2(arc[i], 0.), V2(0., 1.));
            last = arc[i];
        }
    return out;
}

std::string ViewCutter::rim_error() const
{
    return _u8L("The cut ends inside the object. Increase or remove the limits, or start and end the path beyond the part.");
}

V2 ViewCutter::chart(const V3 &p) const
{
    size_t seg; double t, d;
    nearest(frame.to2d(p), seg, t, d);
    return {arc[seg] + t * (arc[seg + 1] - arc[seg]), frame.height(p)};
}

// The offset direction is blended between the vertex normals so the face follows the bends smoothly.
V3 ViewCutter::lift(const V2 &sh, int level) const
{
    const size_t i  = std::min(size_t(std::upper_bound(arc.begin(), arc.end(), sh.x()) - arc.begin()), arc.size() - 1);
    const size_t i0 = i == 0 ? 0 : i - 1;
    const double l  = arc[i0 + 1] - arc[i0];
    const double t  = l > 0. ? std::clamp((sh.x() - arc[i0]) / l, 0., 1.) : 0.;
    const V2     n  = (vertex_normal(i0) * (1. - t) + vertex_normal(i0 + 1) * t).normalized();
    const V2     q  = path[i0] + (path[i0 + 1] - path[i0]) * t + n * (0.5 * kerf * level);
    return frame.to3d(q, sh.y());
}


// First crossing of two non-adjacent segments of an open polyline, searching from the start: the
// segments i < j and the crossing point. Returns false when the polyline does not cross itself.
bool first_self_crossing(const std::vector<V2> &q, size_t &si, size_t &sj, V2 &x)
{
    auto cross = [](const V2 &a, const V2 &b) { return a.x() * b.y() - a.y() * b.x(); };
    const size_t n = q.size();
    for (size_t i = 0; i + 1 < n; ++i) {
        const V2 a = q[i], b = q[i + 1];
        const V2 lo = a.cwiseMin(b), hi = a.cwiseMax(b);
        for (size_t j = n - 2; j >= i + 2 && j + 1 < n; --j) { // the farthest crossing removes the largest loop
            const V2 c = q[j], d = q[j + 1];
            if (std::max(c.x(), d.x()) < lo.x() || std::min(c.x(), d.x()) > hi.x() || std::max(c.y(), d.y()) < lo.y() ||
                std::min(c.y(), d.y()) > hi.y())
                continue;
            const double d1 = cross(b - a, c - a), d2 = cross(b - a, d - a), d3 = cross(d - c, a - c), d4 = cross(d - c, b - c);
            if (((d1 > 0.) != (d2 > 0.)) && ((d3 > 0.) != (d4 > 0.)) && d1 != 0. && d2 != 0. && d3 != 0. && d4 != 0.) {
                si = i; sj = j;
                x  = c + (d - c) * (d3 / (d3 - d4));
                return true;
            }
        }
    }
    return false;
}

// Cuts the loops out of a polyline that crosses itself (hand jitter between close points, which the
// spline turns into tiny curls). Returns the total length removed.
double remove_self_loops(std::vector<V2> &q)
{
    double removed = 0.;
    size_t i, j;
    V2     x;
    while (q.size() > 3 && first_self_crossing(q, i, j, x)) {
        for (size_t k = i + 1; k <= j; ++k) removed += (q[k] - q[k - 1]).norm();
        std::vector<V2> out(q.begin(), q.begin() + long(i) + 1);
        out.push_back(x);
        out.insert(out.end(), q.begin() + long(j) + 1, q.end());
        q.swap(out);
    }
    return removed;
}


std::vector<V3> smooth_path(const std::vector<V3> &points) { return catmull_rom(points, false, 12); }

ViewCutter make_view_cutter(const std::vector<V3> &points, const V3 &view_dir, const V3 &bb_min, const V3 &bb_max, double kerf,
                            const ViewLimits &limits)
{
    ViewCutter vc;
    vc.kerf           = kerf;
    const double diag = (bb_max - bb_min).norm();
    vc.frame          = frame_from_normal(view_dir, (bb_min + bb_max) * 0.5);

    std::vector<V2> q;
    for (const V3 &p : points)
        q.push_back(vc.frame.to2d(p));
    q = remove_close_points(q, 1e-3, false);
    if (q.size() < 2) {
        vc.error = _u8L("The path needs at least two points.");
        return vc;
    }
    if (q.size() >= 3)
        q = remove_close_points(catmull_rom(q, false, 16), 1e-4, false);
    // Hand jitter between points picked close together curls the spline back on itself: cut the curls
    // out. A path that really doubles back (points picked out of order) is refused instead of guessed at.
    {
        double length = 0.;
        for (size_t i = 1; i < q.size(); ++i) length += (q[i] - q[i - 1]).norm();
        if (remove_self_loops(q) > std::max(3., 0.1 * length)) {
            vc.error = _u8L("The path crosses itself seen from the view. Place the points in order along the cut line.");
            return vc;
        }
        q = remove_close_points(q, 1e-4, false);
        if (q.size() < 2) {
            vc.error = _u8L("The path needs at least two points.");
            return vc;
        }
    }

    // Extend both ends: well outside the object unless limited.
    const V2 t0 = (q[0] - q[1]).normalized(), t1 = (q.back() - q[q.size() - 2]).normalized();
    vc.limit_start = limits.start >= 0.;
    vc.limit_end   = limits.end >= 0.;
    q.insert(q.begin(), q.front() + t0 * (vc.limit_start ? std::max(limits.start, 1e-3) : diag));
    q.push_back(q.back() + t1 * (vc.limit_end ? std::max(limits.end, 1e-3) : diag));
    vc.path = q;

    // The extensions can still cross the path (a path that curls round): the wall would fold over itself.
    size_t ci, cj;
    V2     cx;
    if (first_self_crossing(q, ci, cj, cx)) {
        vc.error = _u8L("The path crosses itself seen from the view. Place the points in order along the cut line.");
        return vc;
    }

    if ((q.back() - q.front()).norm() < 1e-6) {
        vc.error = _u8L("The path must not end where it starts.");
        return vc;
    }

    vc.limit_depth = limits.front >= 0. || limits.back >= 0.;
    if (vc.limit_depth) {
        double lo = std::numeric_limits<double>::max(), hi = -lo;
        for (const V3 &p : points) { lo = std::min(lo, vc.frame.height(p)); hi = std::max(hi, vc.frame.height(p)); }
        vc.h_min = limits.front >= 0. ? lo - limits.front : lo - diag;
        vc.h_max = limits.back >= 0. ? hi + limits.back : hi + diag;
    }

    vc.arc.assign(q.size(), 0.);
    for (size_t i = 1; i < q.size(); ++i)
        vc.arc[i] = vc.arc[i - 1] + (q[i] - q[i - 1]).norm();
    // Grid over the projection of the bounding box.
    V2 mn = V2::Constant(std::numeric_limits<double>::max()), mx = -mn;
    for (int c = 0; c < 8; ++c) {
        const V2 p = vc.frame.to2d(V3(c & 1 ? bb_max.x() : bb_min.x(), c & 2 ? bb_max.y() : bb_min.y(), c & 4 ? bb_max.z() : bb_min.z()));
        mn = mn.cwiseMin(p); mx = mx.cwiseMax(p);
    }
    vc.build_grid(mn, mx);
    return vc;
}


// ------------------------------------------------------------------------------------------------ loop

void LoopCutter::spoke_at(double ang, double &r, double &h) const
{
    const size_t n  = angles.size();
    ang             = angles[0] + std::fmod(std::fmod(ang - angles[0], 2. * PI) + 2. * PI, 2. * PI);
    const size_t i  = size_t(std::upper_bound(angles.begin(), angles.end(), ang) - angles.begin());
    const size_t i0 = (i + n - 1) % n, i1 = i % n;
    double       a0 = angles[i0], a1 = angles[i1];
    if (a1 <= a0) a1 += 2. * PI;
    if (ang < a0) ang += 2. * PI;
    const double t = a1 - a0 > 1e-12 ? (ang - a0) / (a1 - a0) : 0.;
    r = radius[i0] * (1. - t) + radius[i1] * t;
    h = hloop[i0] * (1. - t) + hloop[i1] * t;
}

// Smoothstep from the mean height at the centre to the loop's height at the loop, flat beyond it: no cone
// tip at the centre and no crease along the loop.
double LoopCutter::cap_height(const V2 &q, double *outside) const
{
    const V2     d  = q - center;
    const double rr = d.norm();
    double       r, h;
    spoke_at(std::atan2(d.y(), d.x()), r, h);
    if (outside) *outside = rr - (r + margin);
    const double s = r > 1e-12 ? std::min(1., rr / r) : 1.;
    const double w = s * s * (3. - 2. * s);
    return hmean * (1. - w) + h * w;
}

int LoopCutter::side(const V3 &p) const
{
    double       out = 0.;
    const double dz  = frame.height(p) - cap_height(frame.to2d(p), &out);
    if (out > 0. || std::abs(dz) > 2. * kerf + 1e-3)
        return 0;
    return dz > 0. ? 1 : -1;
}

std::string LoopCutter::rim_error() const
{
    return _u8L("The loop does not go all the way around the part. Place the points around it, or widen the cap margin.");
}

LoopCutter make_loop_cutter(const std::vector<V3> &loop_points, double kerf, double margin, const std::function<V3(const V3 &)> &snap,
                            bool flat)
{
    LoopCutter lc;
    lc.kerf   = kerf;
    lc.margin = margin;
    const std::vector<V3> pts = remove_close_points(loop_points, 1e-3, true);
    if (pts.size() < 3) {
        lc.error = _u8L("The loop needs at least three points.");
        return lc;
    }
    std::vector<V3> dense = catmull_rom(pts, true, 16);
    for (V3 &p : dense) p = snap(p);
    dense   = remove_close_points(dense, 1e-4, true);
    lc.loop = dense;

    // Mean plane of the loop: the direction it encloses the most area seen from (its vector area). Unlike a
    // least-squares plane this stays right for strongly non-planar (saddle) loops.
    V3 c = V3::Zero();
    for (const V3 &p : dense) c += p;
    c /= double(dense.size());
    V3 varea = V3::Zero();
    for (size_t i = 0; i < dense.size(); ++i) varea += (dense[i] - c).cross(dense[(i + 1) % dense.size()] - c);
    const std::string degenerate = _u8L("The loop encloses no area. Place the points around the part.");
    if (varea.norm() < 1e-9) {
        lc.error = degenerate;
        return lc;
    }
    lc.frame = frame_from_normal(varea, c);

    std::vector<V2>     q;
    std::vector<double> h;
    for (const V3 &p : dense) { q.push_back(lc.frame.to2d(p)); h.push_back(lc.frame.height(p)); }
    const double area = polygon_area(q);
    if (std::abs(area) < 1e-9) {
        lc.error = degenerate;
        return lc;
    }
    if (area < 0.) { std::reverse(q.begin(), q.end()); std::reverse(h.begin(), h.end()); }
    V2 cc = V2::Zero();
    for (size_t i = 0; i < q.size(); ++i) {
        const V2 &a = q[i], &b = q[(i + 1) % q.size()];
        cc += (a + b) * (a.x() * b.y() - b.x() * a.y());
    }
    lc.center = cc / (6. * std::abs(area));

    // The cap is a height field over spokes from the centre, so the loop must wind once around it. Snapping
    // the smoothed loop to a bumpy surface makes it wobble back here and there: those few points are dropped.
    const size_t        n = q.size();
    std::vector<double> ang(n), rad(n);
    double              winding = 0.;
    for (size_t i = 0; i < n; ++i) {
        const V2 d = q[i] - lc.center;
        ang[i] = std::atan2(d.y(), d.x());
        rad[i] = d.norm();
    }
    for (size_t i = 0; i < n; ++i)
        winding += std::remainder(ang[(i + 1) % n] - ang[i], 2. * PI);
    const size_t first = size_t(std::min_element(ang.begin(), ang.end()) - ang.begin());
    std::rotate(ang.begin(), ang.begin() + long(first), ang.end());
    std::rotate(rad.begin(), rad.begin() + long(first), rad.end());
    std::rotate(h.begin(), h.begin() + long(first), h.end());
    std::vector<double> a2, r2, h2;
    for (size_t i = 0; i < n; ++i)
        if (rad[i] > 1e-6 && (a2.empty() || ang[i] > a2.back() + 1e-9)) { a2.push_back(ang[i]); r2.push_back(rad[i]); h2.push_back(h[i]); }
    if (std::abs(winding) < PI || a2.size() < 8 || a2.size() * 2 < n) {
        lc.error = _u8L("The loop folds back on itself seen from its mean plane. Place the points more evenly around the part.");
        return lc;
    }
    // Low-pass the radius and height profiles (truncated Fourier series): the snapped loop is noisy, and a
    // smooth cap needs far fewer triangles and cuts more cleanly.
    auto smooth = [&a2](const std::vector<double> &y, int harmonics) {
        const size_t        m = a2.size();
        std::vector<double> w(m), out(m, 0.);
        for (size_t i = 0; i < m; ++i) w[i] = 0.5 * std::fmod(a2[(i + 1) % m] - a2[(i + m - 1) % m] + 4. * PI, 2. * PI); // trapezoid weights
        for (int k = 0; k <= harmonics; ++k) {
            double ck = 0., sk = 0.;
            for (size_t i = 0; i < m; ++i) { ck += w[i] * y[i] * std::cos(k * a2[i]); sk += w[i] * y[i] * std::sin(k * a2[i]); }
            ck /= k == 0 ? 2. * PI : PI;
            sk /= PI;
            for (size_t i = 0; i < m; ++i) out[i] += ck * std::cos(k * a2[i]) + sk * std::sin(k * a2[i]);
        }
        return out;
    };
    lc.radius = smooth(r2, 8);
    lc.hloop  = smooth(h2, 4);
    lc.angles = std::move(a2);
    if (*std::min_element(lc.radius.begin(), lc.radius.end()) < 1e-3) {
        lc.error = _u8L("The loop is too irregular. Place the points more evenly around the part.");
        lc.angles.clear();
        return lc;
    }
    double hm = 0.;
    for (double x : lc.hloop) hm += x;
    lc.hmean = hm / double(lc.hloop.size());
    if (flat)
        std::fill(lc.hloop.begin(), lc.hloop.end(), lc.hmean);
    return lc;
}

bool fit_loop_margin(const indexed_triangle_set &its, LoopCutter &lc, double max_margin)
{
    constexpr double bin = 0.25, gap = 2.;
    const size_t     nbins = size_t(max_margin / bin) + 1;
    std::vector<int> hits(nbins, 0);
    auto record = [&](const V3 &p) {
        const double rr = lc.outside(p) + lc.margin; // distance outside the loop itself
        if (rr >= 0. && rr < max_margin) ++hits[size_t(rr / bin)];
    };
    const double        band = 0.5 * lc.kerf;
    std::vector<double> f(its.vertices.size());
    for (size_t i = 0; i < f.size(); ++i) f[i] = lc.field(its.vertices[i].cast<double>());
    // Where the surface meets the cap: vertices inside the slab, and the crossings of the cap's middle
    // sampled along each triangle (the faces of low-poly models are much larger than a bin).
    for (const auto &t : its.indices) {
        if (std::min({f[t[0]], f[t[1]], f[t[2]]}) > band || std::max({f[t[0]], f[t[1]], f[t[2]]}) < -band) continue;
        V3  x[2];
        int nx = 0;
        for (int e = 0; e < 3; ++e) {
            const double fa = f[t[e]], fb = f[t[(e + 1) % 3]];
            const V3     a = its.vertices[t[e]].cast<double>(), b = its.vertices[t[(e + 1) % 3]].cast<double>();
            if (std::abs(fa) <= band) record(a);
            if ((fa < 0.) != (fb < 0.) && nx < 2) x[nx++] = a + (b - a) * (fa / (fa - fb));
        }
        if (nx == 2) {
            const int steps = std::min(4000, int((x[1] - x[0]).norm() / bin) + 1);
            for (int k = 0; k <= steps; ++k) record(x[0] + (x[1] - x[0]) * (double(k) / steps));
        }
    }
    const size_t need = size_t(gap / bin);
    size_t       run  = 0;
    for (size_t b = size_t(std::max(0., lc.margin - 0.5 * gap) / bin); b < nbins; ++b) {
        run = hits[b] == 0 ? run + 1 : 0;
        if (run >= need) {
            lc.margin = (double(b + 1) - 0.5 * double(need)) * bin;
            return true;
        }
    }
    return false;
}


// ------------------------------------------------------------------------------------------------ cleanup

// Sculpted / AI meshes often carry zero-thickness sheets (the same triangle twice with opposite
// winding) and surfaces pinched along an edge. Removes duplicate faces and detaches the extra faces
// of an edge shared by more than two faces by giving them their own copies of the edge vertices.
// Returns the number of fixes.
size_t make_edge_manifold(indexed_triangle_set &its)
{
    size_t fixes = 0;
    // 1) duplicate faces: same vertex set. Opposite winding -> zero-volume sheet, drop both.
    {
        std::vector<std::pair<std::array<int, 3>, size_t>> keyed(its.indices.size());
        for (size_t i = 0; i < its.indices.size(); ++i) {
            std::array<int, 3> k{its.indices[i][0], its.indices[i][1], its.indices[i][2]};
            std::sort(k.begin(), k.end());
            keyed[i] = {k, i};
        }
        tbb::parallel_sort(keyed.begin(), keyed.end());
        std::vector<bool> drop(its.indices.size(), false);
        auto same_winding = [&](size_t a, size_t b) {
            const auto &A = its.indices[a], &B = its.indices[b];
            for (int r = 0; r < 3; ++r)
                if (A[0] == B[r] && A[1] == B[(r + 1) % 3] && A[2] == B[(r + 2) % 3]) return true;
            return false;
        };
        for (size_t i = 0; i < keyed.size();) {
            size_t j = i + 1;
            while (j < keyed.size() && keyed[j].first == keyed[i].first) ++j;
            if (j - i >= 2) {
                for (size_t k = i + 1; k < j; ++k) drop[keyed[k].second] = true;
                if (!same_winding(keyed[i].second, keyed[i + 1].second)) drop[keyed[i].second] = true;
                ++fixes;
            }
            i = j;
        }
        if (fixes > 0) {
            size_t w = 0;
            for (size_t i = 0; i < its.indices.size(); ++i)
                if (!drop[i]) its.indices[w++] = its.indices[i];
            its.indices.resize(w);
        }
    }
    // 2) Edges shared by more than two faces (bodies touching along an edge). Around such an edge the
    //    faces, in angular order, pair up as the two sides of each solid wedge. Around each vertex of
    //    such an edge, faces linked through ordinary edges or through a wedge pair form one sheet; every
    //    sheet but the first gets its own copy of the vertex. Nothing is opened, the bodies just stop
    //    sharing vertices.
    std::vector<std::pair<uint64_t, uint32_t>> edges;
    edges.reserve(its.indices.size() * 3);
    for (size_t i = 0; i < its.indices.size(); ++i)
        for (int e = 0; e < 3; ++e)
            edges.emplace_back(edge_key(its.indices[i][e], its.indices[i][(e + 1) % 3]), uint32_t(i));
    tbb::parallel_sort(edges.begin(), edges.end());
    auto faces_on_edge = [&edges](uint64_t k) {
        const auto lo = std::lower_bound(edges.begin(), edges.end(), std::make_pair(k, uint32_t(0)));
        size_t     n  = 0;
        for (auto it = lo; it != edges.end() && it->first == k; ++it) ++n;
        return n;
    };
    auto vpos = [&its](int v) { return V3(its.vertices[v].cast<double>()); };
    std::vector<int>                                                   pinched;
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> wedges; // pinched edge -> face pairs
    for (size_t i = 0; i < edges.size();) {
        size_t j = i + 1;
        while (j < edges.size() && edges[j].first == edges[i].first) ++j;
        if (j - i > 2) {
            ++fixes;
            const int a = int(edges[i].first >> 32), b = int(edges[i].first & 0xffffffffu);
            pinched.push_back(a);
            pinched.push_back(b);
            const V3 pa = vpos(a), d = (vpos(b) - pa).normalized();
            const V3 e1 = (std::abs(d.x()) < 0.9 ? V3::UnitX() : V3::UnitY()).cross(d).normalized(), e2 = d.cross(e1);
            struct Side { double angle; uint32_t face; bool solid_ahead; };
            std::vector<Side> sides;
            for (size_t k = i; k < j; ++k) {
                const auto &t = its.indices[edges[k].second];
                int c = t[0];
                for (int q = 0; q < 3; ++q) if (t[q] != a && t[q] != b) c = t[q];
                V3 u = vpos(c) - pa;
                u -= d * u.dot(d);
                const V3 n = (vpos(t[1]) - vpos(t[0])).cross(vpos(t[2]) - vpos(t[0]));
                // The face's normal points away from the wedge ahead of it (larger angle) when that wedge is solid.
                sides.push_back({std::atan2(u.dot(e2), u.dot(e1)), edges[k].second, d.cross(u).dot(n) < 0.});
            }
            std::sort(sides.begin(), sides.end(), [](const Side &l, const Side &r) { return l.angle < r.angle; });
            auto &pairs = wedges[edges[i].first];
            for (size_t k = 0; k < sides.size(); ++k)
                if (sides[k].solid_ahead) {
                    const uint32_t f0 = sides[k].face, f1 = sides[(k + 1) % sides.size()].face;
                    pairs.emplace_back(std::min(f0, f1), std::max(f0, f1));
                }
        }
        i = j;
    }
    if (pinched.empty())
        return fixes;
    std::sort(pinched.begin(), pinched.end());
    pinched.erase(std::unique(pinched.begin(), pinched.end()), pinched.end());
    const int        nv0 = int(its.vertices.size());
    std::vector<int> copy_of; // original vertex of each copy made below
    auto original = [&](int w) { return w < nv0 ? w : copy_of[size_t(w - nv0)]; };
    std::unordered_map<int, std::vector<uint32_t>> around; // pinched vertex -> its faces
    for (int v : pinched) around[v];
    for (size_t f = 0; f < its.indices.size(); ++f)
        for (int e = 0; e < 3; ++e) {
            auto it = around.find(its.indices[f][e]);
            if (it != around.end()) it->second.push_back(uint32_t(f));
        }
    for (int v : pinched) {
        const std::vector<uint32_t> &fs = around[v];
        std::vector<int>             sheet(fs.size());
        for (size_t i = 0; i < fs.size(); ++i) sheet[i] = int(i);
        auto find  = [&sheet](int x) { while (sheet[x] != x) x = sheet[x] = sheet[sheet[x]]; return x; };
        auto unite = [&](int x, int y) { sheet[find(x)] = find(y); };
        // Faces around v grouped by their other vertex on each edge through v (w may already be a copy
        // made for an earlier pinched vertex: edges are looked up by original indices).
        std::unordered_map<int, std::vector<int>> by_edge;
        for (size_t i = 0; i < fs.size(); ++i)
            for (int e = 0; e < 3; ++e) {
                const int w = its.indices[fs[i]][e];
                if (w != v) by_edge[original(w)].push_back(int(i));
            }
        for (const auto &kv : by_edge) {
            const uint64_t          ek    = edge_key(v, kv.first);
            const std::vector<int> &local = kv.second;
            if (faces_on_edge(ek) == 2) {
                for (size_t k = 1; k < local.size(); ++k) unite(local[0], local[k]);
            } else if (auto it = wedges.find(ek); it != wedges.end()) {
                for (const auto &pr : it->second) {
                    int x = -1, y = -1;
                    for (int l : local) {
                        if (fs[l] == pr.first) x = l;
                        if (fs[l] == pr.second) y = l;
                    }
                    if (x >= 0 && y >= 0) unite(x, y);
                }
            }
        }
        std::unordered_map<int, int> copy; // sheet root -> vertex index used by that sheet
        for (size_t i = 0; i < fs.size(); ++i) {
            auto [it, first] = copy.emplace(find(int(i)), v);
            if (first && copy.size() > 1) {
                it->second = int(its.vertices.size());
                its.vertices.push_back(its.vertices[v]);
                copy_of.push_back(v);
            }
            for (int e = 0; e < 3; ++e)
                if (its.indices[fs[i]][e] == v) its.indices[fs[i]][e] = it->second;
        }
    }
    return fixes;
}


CutResult cut_with_slab(indexed_triangle_set its, const Cutter &cutter, double refine, const std::function<bool()> &canceled)
{
    CutResult    res;
    const bool   bounded = cutter.limited(); // the wall only exists where cutter.outside(p) < 0
    const double hk = 0.5 * cutter.kerf;
    auto stop = [&]() { return canceled && canceled(); };
    auto pos  = [&](int v) -> V3 { return its.vertices[v].template cast<double>(); };

    // Per-vertex field and distance outside the limited wall.
    std::vector<double> f, out;
    auto eval = [&](size_t from) {
        f.resize(its.vertices.size());
        out.resize(its.vertices.size(), -1.);
        for (size_t i = from; i < its.vertices.size(); ++i) {
            if ((i & 0xffff) == 0 && stop()) return;
            const V3 p = pos(int(i));
            f[i]       = cutter.field(p);
            if (bounded) out[i] = cutter.outside(p);
        }
    };
    eval(0);
    if (stop()) return {};

    // 1) Refine the triangles that may touch the slab where the cut would not follow the slab closely
    //    enough: where the slab is curved across the triangle (the field is not linear along its edges,
    //    so straight cut edges between the crossings would stray from it) and, for a limited wall, where
    //    the triangle runs across its rim. Elsewhere the model's own triangles already
    //    follow a planar slab exactly; refining them only multiplies the triangles of the result.
    //    The fields grow at most ~1 mm per mm, so a triangle whose closest vertex is more than ~2 edge
    //    lengths away from the slab cannot reach it.
    constexpr double bend_tol = 0.05; // mm: straight cut edges stray at most about this far from the slab (a quarter layer)
    auto needs_split = [&](const stl_triangle_vertex_indices &t) { // longest edge to split, or -1
        double emax = 0.; int ei = 0;
        for (int e = 0; e < 3; ++e) {
            const double l = (pos(t[e]) - pos(t[(e + 1) % 3])).squaredNorm();
            if (l > emax) { emax = l; ei = e; }
        }
        emax = std::sqrt(emax);
        if (emax <= std::min(refine, 0.25)) return -1;
        // (not called `near`: windows.h defines that as a macro)
        double closest = std::numeric_limits<double>::max();
        for (int e = 0; e < 3; ++e)
            closest = std::min(closest, std::max(std::abs(f[t[e]]) - hk, out[t[e]]));
        if (closest >= 2. * emax) return -1;
        if (bounded && emax > refine) {
            // The rim only matters where the triangle meets the wall: decide inside or outside at the
            // crossings with the faces of the wall and at the vertices inside it, not at the vertices
            // alone. A large triangle can cross the wall inside the limits with all its vertices outside
            // them, and the rim can run along the surface far from the wall.
            bool inside = false, outside = false;
            for (int e = 0; e < 3; ++e) {
                const int    a = t[e], b = t[(e + 1) % 3];
                const double fa = f[a], fb = f[b];
                if (std::abs(fa) < hk) (out[a] < 0. ? inside : outside) = true;
                for (const double level : {hk, -hk})
                    if ((fa - level) * (fb - level) < 0.)
                        (cutter.outside(pos(a) + (pos(b) - pos(a)) * ((level - fa) / (fb - fa))) < 0. ? inside : outside) = true;
            }
            if (inside && (outside || out[t[0]] >= 0. || out[t[1]] >= 0. || out[t[2]] >= 0.)) return ei;
        }
        // Only edges a face of the slab crosses get a cut vertex: test those (the fields may also be
        // clamped far from the slab, which is not a bend).
        for (int e = 0; e < 3; ++e) {
            const int    a = t[e], b = t[(e + 1) % 3];
            const double fa = f[a], fb = f[b];
            if (((fa - hk) * (fb - hk) >= 0.) && ((fa + hk) * (fb + hk) >= 0.)) continue;
            if (std::abs(cutter.field((pos(a) + pos(b)) * 0.5) - 0.5 * (fa + fb)) > bend_tol) return ei;
        }
        return -1;
    };
    std::vector<size_t> cand;
    for (size_t i = 0; i < its.indices.size(); ++i) {
        if ((i & 0xffff) == 0 && stop()) return {};
        if (needs_split(its.indices[i]) >= 0) cand.push_back(i);
    }
    // Each pass halves the longest edge of every triangle that still needs it, and no edge is split below
    // the refine size, so this ends; the bound only catches a refinement that runs away all the same.
    for (int pass = 0; !cand.empty(); ++pass) {
        std::unordered_map<uint64_t, int> mid;
        for (size_t i : cand) {
            const auto &t  = its.indices[i];
            const int   ei = needs_split(t);
            if (ei >= 0) mid.emplace(edge_key(t[ei], t[(ei + 1) % 3]), -1);
        }
        if (mid.empty())
            break;
        if (pass == 64) {
            res.error = _u8L("The mesh could not be refined enough along the cut. Repair or remesh the model, or cut it elsewhere.");
            return res;
        }
        const size_t nv0 = its.vertices.size();
        std::vector<char> touched(nv0, 0);
        for (auto &[k, idx] : mid) {
            const int a = int(k >> 32), b = int(k & 0xffffffffu);
            touched[a] = touched[b] = 1;
            idx = int(its.vertices.size());
            its.vertices.emplace_back(((pos(a) + pos(b)) * 0.5).template cast<float>());
        }
        eval(nv0);
        // Only triangles with two touched vertices can hold a marked edge.
        std::vector<size_t> hit;
        for (size_t i = 0; i < its.indices.size(); ++i) {
            const auto &t = its.indices[i];
            if (int(touched[t[0]]) + int(touched[t[1]]) + int(touched[t[2]]) >= 2) hit.push_back(i);
        }
        const size_t n_before = its.indices.size();
        split_marked(its.indices, [&mid](int a, int b) { auto it = mid.find(edge_key(a, b)); return it == mid.end() ? -1 : it->second; }, &hit);
        cand = std::move(hit);
        for (size_t i = n_before; i < its.indices.size(); ++i) cand.push_back(i);
        if (stop()) return {};
    }

    // Vertices lying (almost) on a level set would get crossings within float precision of them, i.e.
    // coincident vertices and zero-area triangles. Nudge them off the level set, by a micrometre or so.
    {
        constexpr double snap = 2e-3, step = 1e-3;
        for (size_t v = 0; v < its.vertices.size(); ++v) {
            if ((v & 0xffff) == 0 && stop()) return {};
            if (bounded && out[v] >= 0.) continue;
            for (int level : {1, -1}) {
                const double x = level > 0 ? f[v] - hk : -f[v] - hk;
                if (std::abs(x) >= snap) continue;
                const V3 p = pos(int(v));
                V3       grad;
                for (int k = 0; k < 3; ++k) {
                    V3 d = V3::Zero();
                    d[k] = step;
                    grad[k] = (cutter.field(p + d) - cutter.field(p - d)) / (2. * step);
                }
                if (level < 0) grad = -grad;
                const double gn = grad.norm();
                if (gn < 1e-6) continue;
                its.vertices[v] = (p + grad * ((snap - x) / (gn * gn))).template cast<float>();
                f[v]            = cutter.field(pos(int(v)));
                if (bounded) out[v] = cutter.outside(pos(int(v)));
            }
        }
    }
    // 2) Clip against f >= hk (level +1) and f <= -hk (level -1). Values exactly on a level are
    //    nudged to the kept side so every crossing lies strictly inside its edge.
    auto g = [&](int v, int level) {
        const double x = level > 0 ? f[v] - hk : -f[v] - hk;
        return std::abs(x) < 1e-5 ? 1e-5 : x;
    };
    // Parts of the model the limited wall ends inside are not cut (see foreign_bands).
    std::vector<char> foreign;
    if (bounded && cutter.skips_parts_it_ends_in()) foreign = foreign_bands(its, f, out, hk);
    auto active_t = [&](size_t ti) {
        const auto &t = its.indices[ti];
        return !bounded || (out[t[0]] < 0. && out[t[1]] < 0. && out[t[2]] < 0. && (foreign.empty() || !foreign[ti]));
    };
    std::unordered_map<uint64_t, int> cross[2];
    auto crossing = [&](int a, int b, int level) {
        auto &m      = cross[level > 0 ? 0 : 1];
        auto [it, ins] = m.emplace(edge_key(a, b), -1);
        if (ins) {
            // Put the crossing on the level set itself rather than where the linear interpolation of the
            // field puts it, so the cut edges agree with the caps built on the exact level set
            // (Illinois regula falsi, the field is cheap to evaluate).
            const int    lo = std::min(a, b), hi = std::max(a, b);
            const V3     pa = pos(lo), pb = pos(hi);
            auto         gx = [&](double t) { const double x = cutter.field(pa + (pb - pa) * t); return level > 0 ? x - hk : -x - hk; };
            double       t0 = 0., t1 = 1., g0 = g(lo, level), g1 = g(hi, level);
            double       t  = g0 / (g0 - g1);
            for (int iter = 0, side = 0; iter < 40; ++iter) {
                t              = std::clamp((t0 * g1 - t1 * g0) / (g1 - g0), 0., 1.);
                const double gt = gx(t);
                if (std::abs(gt) < 1e-7 || (t1 - t0) < 1e-9) break;
                if ((gt > 0.) == (g0 > 0.)) { t0 = t; g0 = gt; if (side == -1) g1 *= 0.5; side = -1; }
                else                        { t1 = t; g1 = gt; if (side == 1) g0 *= 0.5; side = 1; }
            }
            it->second = int(its.vertices.size());
            its.vertices.emplace_back((pa + (pb - pa) * t).template cast<float>());
        }
        return it->second;
    };
    std::vector<stl_triangle_vertex_indices> kept;
    kept.reserve(its.indices.size() + its.indices.size() / 8);
    std::vector<std::pair<int, int>> cap_edges[2]; // directed so the caps get the right winding
    std::vector<size_t>              inactive_touching;
    std::unordered_map<uint64_t, size_t> cut_owner;
    std::unordered_map<uint64_t, int>    removed; // bounded cut only, values unused
    for (size_t ti = 0; ti < its.indices.size(); ++ti) {
        if ((ti & 0xffff) == 0 && stop()) return {};
        const auto &t = its.indices[ti];
        if (!active_t(ti)) {
            kept.push_back(t);
            inactive_touching.push_back(ti);
            continue;
        }
        if (bounded) // edges with both ends inside the slab disappear; an inactive neighbour must not use them
            for (int e = 0; e < 3; ++e) {
                const int a = t[e], b = t[(e + 1) % 3];
                if (std::abs(f[a]) < hk && std::abs(f[b]) < hk) removed.emplace(edge_key(a, b), 0);
            }
        for (int level : {1, -1}) {
            const bool k0 = g(t[0], level) > 0., k1 = g(t[1], level) > 0., k2 = g(t[2], level) > 0.;
            if (!k0 && !k1 && !k2) continue;
            if (k0 && k1 && k2) { kept.push_back(t); continue; }
            const bool keep[3] = {k0, k1, k2};
            int        poly[4], np = 0, first_cross = -1;
            for (int i = 0; i < 3; ++i) {
                const int a = t[i], b = t[(i + 1) % 3];
                if (keep[i]) poly[np++] = a;
                if (keep[i] != keep[(i + 1) % 3]) {
                    if (keep[i]) first_cross = np; // the cut edge starts here
                    poly[np++] = crossing(a, b, level);
                }
            }
            const size_t base = kept.size();
            for (int i = 1; i + 1 < np; ++i)
                kept.emplace_back(poly[0], poly[i], poly[i + 1]);
            // Kept polygon walks X_out -> X_in along the cut; the cap walks it backwards.
            const int x0 = poly[first_cross], x1 = poly[(first_cross + 1) % np];
            cap_edges[level > 0 ? 0 : 1].emplace_back(x1, x0);
            // Fan triangle that holds the cut edge, in case a cap vertex gets inserted on it.
            const int fan = first_cross == 0 ? 1 : first_cross == np - 1 ? np - 2 : first_cross;
            cut_owner.emplace(edge_key(x0, x1), base + size_t(fan - 1));
        }
    }
    // A triangle outside the limits must not border a clipped one: that means the surface runs through
    // the rim of the wall, i.e. the wall ends inside a part it cuts.
    if (bounded)
        for (size_t k = 0; k < inactive_touching.size(); ++k) {
            if ((k & 0xffff) == 0 && stop()) return {};
            const auto &t = its.indices[inactive_touching[k]];
            for (int e = 0; e < 3; ++e)
                for (auto *m : {&cross[0], &cross[1], &removed})
                    if (m->count(edge_key(t[e], t[(e + 1) % 3]))) {
                        res.error = cutter.rim_error();
                        return res;
                    }
        }
    its.indices   = std::move(kept);
    res.first_cap = its.indices.size();

    // 3) Close the openings.
    switch (close_cut_faces(its, cap_edges, cutter, cut_owner, res, canceled)) {
    case 1: return res;
    case 2: return {};
    default: break;
    }
    res.level.resize(its.vertices.size(), 0);
    for (int li = 0; li < 2; ++li)
        for (const auto &kv : cross[li]) res.level[size_t(kv.second)] = static_cast<signed char>(li == 0 ? 1 : -1);
    res.mesh = std::move(its);
    return res;
}


// Preview of a cut: how far, and which way, the detached piece is drawn apart from the piece that stays,
// along the line between their bounding box centres (straight up when they share a centre): 10 % of the
// detached piece's size, between 6 and 30 mm. Zero when either piece is empty.
V3 preview_shift(const indexed_triangle_set &moved, const indexed_triangle_set &stay)
{
    if (moved.indices.empty() || stay.indices.empty())
        return V3::Zero();
    auto box = [](const indexed_triangle_set &its, V3 &mn, V3 &mx) {
        mn = mx = its.vertices.front().cast<double>();
        for (const auto &v : its.vertices) { mn = mn.cwiseMin(v.cast<double>()); mx = mx.cwiseMax(v.cast<double>()); }
    };
    V3 m0, m1, s0, s1;
    box(moved, m0, m1);
    box(stay, s0, s1);
    V3 dir = (m0 + m1) * 0.5 - (s0 + s1) * 0.5;
    if (dir.norm() < 1e-6)
        dir = V3::UnitZ();
    return dir.normalized() * std::clamp(0.1 * (m1 - m0).norm(), 6., 30.);
}

// Preview of a cut: how far both pieces are drawn raised so that the detached piece, drawn `shift` apart,
// does not go below the lowest point of the object, where the object rests on the bed. When it is
// drawn downwards, it then rests there and the piece that stays is drawn above it. Zero otherwise.
double preview_lift(const indexed_triangle_set &moved, const indexed_triangle_set &stay, const V3 &shift)
{
    if (moved.indices.empty() || stay.indices.empty())
        return 0.;
    auto min_z = [](const indexed_triangle_set &its) {
        double z = std::numeric_limits<double>::max();
        for (const auto &v : its.vertices) z = std::min(z, double(v.z()));
        return z;
    };
    const double moved_z = min_z(moved);
    return std::max(0., std::min(moved_z, min_z(stay)) - (moved_z + shift.z()));
}


// The parts are found with a union-find over the vertices instead of its_split(), which would lose the
// cut face marks (CutResult::level) the parts are assigned by.
Grouped group_cut(const CutResult &cut, const Cutter &cutter, bool *joined)
{
    const indexed_triangle_set &its = cut.mesh;
    std::vector<int>            dsu(its.vertices.size());
    for (size_t i = 0; i < dsu.size(); ++i) dsu[i] = int(i);
    auto find = [&dsu](int x) { while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]]; return x; };
    for (const auto &t : its.indices) { dsu[find(t[0])] = find(t[1]); dsu[find(t[1])] = find(t[2]); }
    std::unordered_map<int, size_t>   comp_of_root;
    std::vector<indexed_triangle_set> parts;
    std::vector<int>                  label; // +1 / -1 by the cut face carried, 2 both, 0 none
    std::vector<int>                  remap(its.vertices.size(), -1);
    if (joined) *joined = false;
    for (const auto &t : its.indices) {
        auto [it, ins] = comp_of_root.emplace(find(t[0]), parts.size());
        if (ins) { parts.emplace_back(); label.push_back(0); }
        const size_t c = it->second;
        stl_triangle_vertex_indices nt;
        for (int k = 0; k < 3; ++k) {
            int &m = remap[size_t(t[k])];
            if (m < 0) { m = int(parts[c].vertices.size()); parts[c].vertices.push_back(its.vertices[size_t(t[k])]); }
            nt[k] = m;
            const int lv = size_t(t[k]) < cut.level.size() ? cut.level[size_t(t[k])] : 0;
            if (lv != 0 && label[c] != 2)
                label[c] = label[c] != 0 && label[c] != lv ? 2 : lv;
        }
        parts[c].indices.push_back(nt);
    }
    // Parts untouched by the cut go to the side most of their vertices are on. Where the cutter cannot tell
    // (away from a loop cap), they go with a part of side a they lie in (a pin in the detached foot), else to b.
    for (size_t i = 0; i < parts.size(); ++i) {
        if (label[i] == 2 && joined) *joined = true;
        if (label[i] != 0) continue;
        const auto  &vs   = parts[i].vertices;
        const size_t step = std::max<size_t>(1, vs.size() / 4000);
        long         vote = 0;
        for (size_t k = 0; k < vs.size(); k += step) vote += cutter.side(vs[k].cast<double>());
        label[i] = vote > 0 ? 1 : vote < 0 ? -1 : 0;
    }
    for (size_t i = 0; i < parts.size(); ++i) {
        if (label[i] != 0) continue;
        label[i] = -1;
        const Vec3d p = parts[i].vertices.front().cast<double>();
        for (size_t j = 0; j < parts.size() && label[i] < 0; ++j)
            if (label[j] == 1 && AABBMesh(parts[j]).query_ray_hit(p, Vec3d(1e-4, 2e-4, 1.).normalized()).is_inside()) label[i] = 1;
    }
    Grouped g;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (std::abs(its_volume(parts[i])) < 1.f) ++g.n_small;
        if (label[i] == 1) { its_merge(g.side_a, parts[i]); ++g.n_a; }
        else               { its_merge(g.side_b, parts[i]); ++g.n_b; }
    }
    return g;
}

}} // namespace Slic3r::PathCut
