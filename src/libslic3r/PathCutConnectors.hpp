#pragma once
// Path cut connectors: straight pegs / holes / dowel pins placed on the curved cut face left by
// PathCutGeometry, so the printed parts register (no shift, no twist) when glued.
//
// The cut face of a path cut is curved, so a connector cannot simply stand on a plane as in the planar
// cut. Every connector here is straight and shares one axis: the mean normal of the cut face. Where it
// may go is decided by casting lines parallel to that axis through each side (a grid of the triangles
// projected along the axis keeps that linear in the triangle count): a connector is placed only where
// both sides have material all around it, deep enough, with a wall left, and where the cut face is flat
// enough under it for the peg to be anchored. Among the places that fit, they are spread out.
//
// Like the planar cut connectors, the result is volumes, not a mesh boolean: a peg is a model part of
// one side, a hole a negative volume of the other (subtracted when slicing), a dowel a separate pin.
//
// Material is told from the crossings along each line taken in pairs, so the model parts must not
// overlap each other (the cut core makes the same assumption).
//
// Only depends on PathCutGeometry.hpp (Eigen, admesh), so it can be unit-tested on its own.

#include "PathCutGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace Slic3r { namespace PathCut {

enum class ConnectorType : int { None = 0, Plug = 1, Dowel = 2 };
enum class ConnectorShape : int { Cylinder = 0, Square = 1 };

struct ConnectorParams
{
    ConnectorType  type{ConnectorType::None};
    ConnectorShape shape{ConnectorShape::Cylinder};
    int            count{2};
    double         size{5.};      // diameter, or side of the square [mm]
    double         depth{5.};     // how far the peg (or each half of the pin) goes into a side [mm]
    double         clearance{0.2}; // added to the hole radius and depth [mm]
    double         wall{1.};      // material left around a hole [mm]
};

struct ConnectorResult
{
    std::vector<V3>      positions;      // on the cut face of the `from` side
    V3                   axis{V3::Zero()}; // from `from` into `to`
    indexed_triangle_set pegs;           // Plug: model part of `from`
    indexed_triangle_set holes_to;       // negative volume of `to`
    indexed_triangle_set holes_from;     // Dowel: negative volume of `from`
    indexed_triangle_set pin;            // Dowel: one pin, standing on z = 0 at the origin
};

namespace detail {

// Lines parallel to an axis through a closed mesh: all crossings, via the triangles projected along the
// axis into a uniform grid.
class LineCaster
{
public:
    LineCaster(const indexed_triangle_set &its, const Frame &f) : m_its(its)
    {
        const size_t nv = its.vertices.size(), nt = its.indices.size();
        m_p.resize(nv);
        m_h.resize(nv);
        V2 mn = V2::Constant(std::numeric_limits<double>::max()), mx = -mn;
        for (size_t i = 0; i < nv; ++i) {
            const V3 p = its.vertices[i].cast<double>();
            m_p[i]     = f.to2d(p);
            m_h[i]     = f.height(p);
            mn         = mn.cwiseMin(m_p[i]);
            mx         = mx.cwiseMax(m_p[i]);
        }
        if (nt == 0) return;
        const V2 ext = (mx - mn).cwiseMax(V2::Constant(1e-6));
        m_cell = std::max(std::sqrt(ext.x() * ext.y() / double(std::max<size_t>(nt / 2, 1))), 1e-3);
        while ((ext.x() / m_cell + 1.) * (ext.y() / m_cell + 1.) > 4e6) m_cell *= 1.5;
        m_mn = mn;
        m_nx = int(ext.x() / m_cell) + 1;
        m_ny = int(ext.y() / m_cell) + 1;
        m_start.assign(size_t(m_nx) * size_t(m_ny) + 1, 0);
        auto range = [this](const stl_triangle_vertex_indices &t, int &x0, int &x1, int &y0, int &y1) {
            V2 a = m_p[size_t(t[0])], b = a;
            for (int k = 1; k < 3; ++k) { a = a.cwiseMin(m_p[size_t(t[k])]); b = b.cwiseMax(m_p[size_t(t[k])]); }
            x0 = std::clamp(int((a.x() - m_mn.x()) / m_cell), 0, m_nx - 1);
            x1 = std::clamp(int((b.x() - m_mn.x()) / m_cell), 0, m_nx - 1);
            y0 = std::clamp(int((a.y() - m_mn.y()) / m_cell), 0, m_ny - 1);
            y1 = std::clamp(int((b.y() - m_mn.y()) / m_cell), 0, m_ny - 1);
        };
        for (const auto &t : its.indices) {
            int x0, x1, y0, y1;
            range(t, x0, x1, y0, y1);
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) ++m_start[size_t(y) * size_t(m_nx) + size_t(x) + 1];
        }
        for (size_t i = 1; i < m_start.size(); ++i) m_start[i] += m_start[i - 1];
        m_tri.resize(m_start.back());
        std::vector<size_t> fill(m_start.begin(), m_start.end() - 1);
        for (size_t ti = 0; ti < nt; ++ti) {
            int x0, x1, y0, y1;
            range(its.indices[ti], x0, x1, y0, y1);
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) m_tri[fill[size_t(y) * size_t(m_nx) + size_t(x)]++] = uint32_t(ti);
        }
    }

    // Heights (along the frame normal) where the line through q crosses the mesh, sorted. Returns false
    // when the count is odd (an open mesh, or an unlucky hit exactly on an edge).
    bool crossings(V2 q, std::vector<double> &out) const
    {
        out.clear();
        if (m_tri.empty()) return true;
        q += V2(1.234567e-7, 7.654321e-8) * m_cell; // off the grid of vertices and edges of regular meshes
        const int x = int((q.x() - m_mn.x()) / m_cell), y = int((q.y() - m_mn.y()) / m_cell);
        if (x < 0 || y < 0 || x >= m_nx || y >= m_ny) return true;
        const size_t c = size_t(y) * size_t(m_nx) + size_t(x);
        for (size_t i = m_start[c]; i < m_start[c + 1]; ++i) {
            const auto  &t = m_its.indices[m_tri[i]];
            const V2    &a = m_p[size_t(t[0])], &b = m_p[size_t(t[1])], &cc = m_p[size_t(t[2])];
            const double d = (b - a).x() * (cc - a).y() - (b - a).y() * (cc - a).x();
            if (std::abs(d) < 1e-14) continue; // edge-on to the line
            const double w1 = ((q - a).x() * (cc - a).y() - (q - a).y() * (cc - a).x()) / d;
            const double w2 = ((b - a).x() * (q - a).y() - (b - a).y() * (q - a).x()) / d;
            if (w1 < 0. || w2 < 0. || w1 + w2 > 1.) continue;
            out.push_back((1. - w1 - w2) * m_h[size_t(t[0])] + w1 * m_h[size_t(t[1])] + w2 * m_h[size_t(t[2])]);
        }
        std::sort(out.begin(), out.end());
        return (out.size() & 1) == 0;
    }

private:
    const indexed_triangle_set &m_its;
    std::vector<V2>             m_p;
    std::vector<double>         m_h;
    V2                          m_mn{V2::Zero()};
    double                      m_cell{1.};
    int                         m_nx{0}, m_ny{0};
    std::vector<size_t>         m_start;
    std::vector<uint32_t>       m_tri;
};

// A closed prism (cylinder with `sides` sides, or a square for 4) around the frame axis through q,
// from height h0 to h1 above the frame origin. `radius` is the circumradius.
inline void append_prism(indexed_triangle_set &its, const Frame &f, const V2 &q, double radius, double h0, double h1, int sides)
{
    const int    base = int(its.vertices.size());
    const double rot  = sides == 4 ? PI_ / 4. : 0.;
    for (int k = 0; k < 2; ++k)
        for (int i = 0; i < sides; ++i) {
            const double a = rot + 2. * PI_ * i / sides;
            its.vertices.emplace_back(f.to3d(q + radius * V2(std::cos(a), std::sin(a)), k == 0 ? h0 : h1).cast<float>());
        }
    its.vertices.emplace_back(f.to3d(q, h0).cast<float>());
    its.vertices.emplace_back(f.to3d(q, h1).cast<float>());
    const int c0 = base + 2 * sides, c1 = c0 + 1;
    for (int i = 0; i < sides; ++i) {
        const int j = (i + 1) % sides;
        const int a0 = base + i, b0 = base + j, a1 = base + sides + i, b1 = base + sides + j;
        its.indices.emplace_back(a0, b0, b1); // sides, counter-clockwise seen from outside (u, v, n right handed)
        its.indices.emplace_back(a0, b1, a1);
        its.indices.emplace_back(c0, b0, a0); // bottom, facing -n
        its.indices.emplace_back(c1, a1, b1); // top, facing +n
    }
}

} // namespace detail

// Triangles of the cut face of one side (cut faces with `level`), and the sum of their area vectors
// (outward of that side). Returns the total area.
inline double cut_face(const CutResult &cut, int level, std::vector<size_t> *face, V3 &nsum)
{
    const indexed_triangle_set &m = cut.mesh;
    double                      area = 0.;
    nsum                             = V3::Zero();
    for (size_t i = cut.first_cap; i < cut.first_cap + cut.cap_triangles && i < m.indices.size(); ++i) {
        const auto &t     = m.indices[i];
        int         votes = 0;
        for (int k = 0; k < 3; ++k)
            if (size_t(t[k]) < cut.level.size() && cut.level[size_t(t[k])] == level) ++votes;
        if (votes < 2) continue;
        const V3 a = m.vertices[size_t(t[0])].cast<double>(), b = m.vertices[size_t(t[1])].cast<double>(), c = m.vertices[size_t(t[2])].cast<double>();
        const V3 n = 0.5 * (b - a).cross(c - a);
        nsum += n;
        area += n.norm();
        if (face) face->push_back(i);
    }
    return area;
}

// Which side gets the pegs (Plug): the one whose cut face looks up, so the pegs point up (or level) as
// the parts lie after the cut and never reach below the part they stick out of, which would lift it off
// the bed; for a vertical cut face, the larger side. Returns the cut face level of that side (+1 = side_a).
inline int peg_side(const CutResult &cut, const indexed_triangle_set &side_a, const indexed_triangle_set &side_b)
{
    V3 n;
    if (cut_face(cut, 1, nullptr, n) > 0. && n.norm() > 0.) {
        const double z = n.normalized().z();
        if (std::abs(z) > 0.05) return z > 0. ? 1 : -1;
    }
    return std::abs(signed_volume(side_a)) >= std::abs(signed_volume(side_b)) ? 1 : -1;
}

// Places connectors across the cut face between `from` and `to`, the two sides of `cut` (cut faces with
// level `from_level` belong to `from`). Pegs (Plug) stick out of `from` into holes of `to`.
// `error` is set when not a single connector fits; fewer than asked is not an error.
inline ConnectorResult place_connectors(const CutResult &cut, int from_level, const indexed_triangle_set &from, const indexed_triangle_set &to,
                                        const ConnectorParams &prm, std::string &error, const std::function<bool()> &canceled = {})
{
    ConnectorResult res;
    error.clear();
    if (prm.type == ConnectorType::None || prm.count <= 0) return res;
    if (!(prm.size > 0.) || !(prm.depth > 0.)) {
        error = "No connectors fit: the cut face is too small or too thin for them. Make them smaller or shorter.";
        return res;
    }
    const indexed_triangle_set &m = cut.mesh;

    // Cut face of `from`: its triangles and their mean outward normal (pointing into `to`).
    std::vector<size_t> face;
    V3                  nsum;
    const double        area = cut_face(cut, from_level, &face, nsum);
    if (face.empty() || !(area > 1e-9) || nsum.norm() <= 0.2 * area) {
        error = "No connectors fit: the cut face is too curved for a straight connector.";
        return res;
    }
    const V3    axis = nsum.normalized();
    const Frame f    = frame_from_normal(axis, V3::Zero());
    res.axis         = axis;

    const int    sides = prm.shape == ConnectorShape::Square ? 4 : 32;
    const double r     = 0.5 * prm.size;
    const double rc    = prm.shape == ConnectorShape::Square ? r * std::sqrt(2.) : r; // circumradius
    const double tol   = std::max(prm.clearance, 0.);
    const double rh    = rc + (prm.shape == ConnectorShape::Square ? tol * std::sqrt(2.) : tol); // hole circumradius
    const double wall  = std::max(prm.wall, 0.);
    const double depth = prm.depth;
    const double embed = std::clamp(0.4 * depth, 1., 3.); // peg base sunk into `from` below the face
    const double ro    = rh + wall;                        // must be solid out to here around a hole
    const bool   dowel = prm.type == ConnectorType::Dowel;

    // Candidates: points on the cut face, one per grid cell across the axis. Large cut face triangles
    // (a flat cap is a few big ones) are sampled inside, not only at their centroid, so a spot that fits
    // is not missed for lack of a candidate near it.
    double cell = std::max(0.25 * r, 0.25);
    std::vector<V3> cand;
    for (int attempt = 0; attempt < 8; ++attempt, cell *= 1.6) {
        cand.clear();
        std::unordered_map<uint64_t, size_t> seen;
        std::vector<double> off; // distance of each candidate to its cell centre: keep the closest
        auto add = [&](const V3 &c) {
            const V2       q   = f.to2d(c);
            const double   gx = std::floor(q.x() / cell), gy = std::floor(q.y() / cell);
            const uint64_t key = (uint64_t(uint32_t(int32_t(gx))) << 32) | uint32_t(int32_t(gy));
            const double   d   = (q - V2((gx + 0.5) * cell, (gy + 0.5) * cell)).norm();
            auto [it, ins]     = seen.emplace(key, cand.size());
            if (ins) { cand.push_back(c); off.push_back(d); }
            else if (d < off[it->second]) { cand[it->second] = c; off[it->second] = d; }
        };
        for (size_t i : face) {
            const auto &t = m.indices[i];
            const V3    a = m.vertices[size_t(t[0])].cast<double>(), b = m.vertices[size_t(t[1])].cast<double>(), c = m.vertices[size_t(t[2])].cast<double>();
            add((a + b + c) / 3.);
            const double longest = std::max({(b - a).norm(), (c - b).norm(), (a - c).norm()});
            const int    k       = std::min(64, int(longest / cell));
            for (int u = 1; u < k; ++u)
                for (int v = 1; u + v < k; ++v)
                    add(a + (b - a) * (double(u) / k) + (c - a) * (double(v) / k));
        }
        if (cand.size() <= 6000) break;
    }

    const detail::LineCaster cast_from(from, f), cast_to(to, f);
    std::vector<double>      hits;
    // Solid interval of a side along the line through q, as heights relative to the face point at h:
    // the one whose top (bottom for `to`) is nearest the face.
    auto interval = [&](const detail::LineCaster &lc, const V2 &q, double h, bool top, double &lo, double &hi) {
        if (!lc.crossings(q, hits) || hits.empty()) return false;
        double best = std::numeric_limits<double>::max();
        for (size_t i = 0; i + 1 < hits.size(); i += 2) {
            const double a = hits[i] - h, b = hits[i + 1] - h, d = std::abs(top ? b : a);
            if (d < best) { best = d; lo = a; hi = b; }
        }
        return true;
    };
    // A place that fits, with the faces measured along its axis: heights relative to p of the face of
    // `from` and of `to` on the axis, and the highest point of the face of `from` around it.
    struct Spot { V3 p; double from_face, to_face, from_top; };
    const int ring_n = 24;
    auto fits = [&](const V3 &p, Spot &spot, double &slack) {
        const V2     q = f.to2d(p);
        const double h = f.height(p);
        spot           = {p, 0., 0., -std::numeric_limits<double>::max()};
        slack          = std::numeric_limits<double>::max();
        for (int ring = 0; ring < 3; ++ring) {
            const double rr = ring == 0 ? 0. : ring == 1 ? 0.5 * ro : ro;
            const int    n  = ring == 0 ? 1 : ring == 1 ? ring_n / 2 : ring_n;
            for (int i = 0; i < n; ++i) {
                const double a  = 2. * PI_ * i / n;
                const V2     qq = q + rr * V2(std::cos(a), std::sin(a));
                double a0 = 0., a1 = 0., b0 = 0., b1 = 0.;
                if (!interval(cast_from, qq, h, true, a0, a1) || !interval(cast_to, qq, h, false, b0, b1)) return false;
                if (ring == 0) { spot.from_face = a1; spot.to_face = b0; }
                spot.from_top = std::max(spot.from_top, a1);
                // `from`: the face is near p (not sloping past the sunk base) and solid down to the peg base,
                // or past the bottom of the hole of a pin.
                const double need_from = dowel ? spot.from_face - depth - tol - wall : -embed - 0.2;
                if (a1 < -embed + 0.3 || a1 > 0.5 * depth || a0 > need_from) return false;
                // `to`: its face is above the peg base and the peg goes in at least half way; solid past
                // the bottom of the hole, which is `depth` (+ clearance) into `to`, from its face.
                const double need_to = spot.to_face + depth + tol + wall;
                if (b0 < -embed || b0 > 0.5 * depth || b1 < need_to) return false;
                slack = std::min({slack, need_from - a0, b1 - need_to, a1 + embed - 0.3});
            }
        }
        return true;
    };

    std::vector<Spot>   ok;
    std::vector<double> ok_slack;
    for (size_t i = 0; i < cand.size(); ++i) {
        if ((i & 63) == 0 && canceled && canceled()) return res;
        Spot   spot;
        double s;
        if (fits(cand[i], spot, s)) { ok.push_back(spot); ok_slack.push_back(s); }
    }
    if (ok.empty()) {
        error = "No connectors fit: the cut face is too small or too thin for them. Make them smaller or shorter.";
        return res;
    }
    // First the roomiest spot, then each next one as far as possible from those already placed.
    const double  min_gap = 2. * ro + 1.;
    std::vector<double> dist(ok.size(), std::numeric_limits<double>::max());
    size_t pick = size_t(std::max_element(ok_slack.begin(), ok_slack.end()) - ok_slack.begin());
    std::vector<Spot> placed;
    while (int(placed.size()) < prm.count) {
        placed.push_back(ok[pick]);
        res.positions.push_back(ok[pick].p);
        const V2 last = f.to2d(ok[pick].p);
        double   best = -1.;
        for (size_t i = 0; i < ok.size(); ++i) {
            dist[i] = std::min(dist[i], (f.to2d(ok[i].p) - last).norm());
            if (dist[i] >= min_gap && dist[i] > best) { best = dist[i]; pick = i; }
        }
        if (best < 0.) break;
    }

    // Depths count from each side's own face (the faces are a kerf apart): a pin of 2 x depth fits two
    // holes of depth + clearance, a peg reaching depth into `to` a hole of depth + clearance.
    for (const Spot &s : placed) {
        const V2     q = f.to2d(s.p);
        const double h = f.height(s.p);
        // Holes start below / end above the face they open on everywhere around (negative volumes only
        // remove material of their own part, so overshooting into the gap is harmless).
        detail::append_prism(res.holes_to, f, q, rh, h - embed - 1., h + s.to_face + depth + tol, sides);
        if (dowel)
            detail::append_prism(res.holes_from, f, q, rh, h + s.from_face - depth - tol, h + s.from_top + 1., sides);
        else
            detail::append_prism(res.pegs, f, q, rc, h - embed, h + s.to_face + depth, sides);
    }
    if (dowel) {
        Frame up;
        detail::append_prism(res.pin, up, V2::Zero(), rc, 0., 2. * depth, sides);
    }
    return res;
}

}} // namespace Slic3r::PathCut
