#pragma once
// Path cut ("corte por ruta"): geometry core shared by GLGizmoPathCut and its tests.
//
// Two ways to split a mesh along a user path instead of a plane:
//  * View path:  an open polyline picked from one view direction. It is smoothed, extended past the
//                object and extruded along the view direction into a thin wall.
//  * Loop:       a closed loop picked on the surface around the region to detach (e.g. an arm).
//                A curved cap spanning the loop is thickened into a thin slab that only cuts there.
// In both cases the cutter thickness is the clearance ("kerf") left between the resulting parts.
// The boolean itself (mesh minus cutter) is done by the caller; this file builds the cutters and
// assigns the resulting connected components to one side or the other.
//
// Only depends on Eigen and admesh, so it can be unit-tested without the rest of libslic3r.

#include <admesh/stl.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <map>
#include <cmath>
#include <string>
#include <vector>

namespace Slic3r { namespace PathCut {

using V2 = Eigen::Matrix<double, 2, 1, Eigen::DontAlign>;
using V3 = Eigen::Matrix<double, 3, 1, Eigen::DontAlign>;

constexpr double PI_ = 3.14159265358979323846;

struct Frame
{
    V3 c{V3::Zero()}, u{V3::UnitX()}, v{V3::UnitY()}, n{V3::UnitZ()};
    V2     to2d(const V3 &p) const { const V3 d = p - c; return {d.dot(u), d.dot(v)}; }
    double height(const V3 &p) const { return (p - c).dot(n); }
    V3     to3d(const V2 &q, double h) const { return c + q.x() * u + q.y() * v + h * n; }
};

inline Frame frame_from_normal(const V3 &normal, const V3 &center)
{
    Frame f;
    f.c       = center;
    f.n       = normal.normalized();
    V3 a      = std::abs(f.n.x()) < 0.9 ? V3::UnitX() : V3::UnitY();
    f.u       = f.n.cross(a).normalized();
    f.v       = f.n.cross(f.u);
    return f;
}

// Centripetal-free uniform Catmull-Rom spline through all control points.
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

// ------------------------------------------------------------------------------------------------ mesh helpers
inline int add_vertex(indexed_triangle_set &its, const V3 &p)
{
    its.vertices.emplace_back(p.cast<float>());
    return int(its.vertices.size()) - 1;
}
inline void add_tri(indexed_triangle_set &its, int a, int b, int c) { its.indices.emplace_back(a, b, c); }
inline void add_quad(indexed_triangle_set &its, int a, int b, int c, int d) { add_tri(its, a, b, c); add_tri(its, a, c, d); }

inline double signed_volume(const indexed_triangle_set &its)
{
    double vol = 0.;
    for (const auto &f : its.indices) {
        const V3 a = its.vertices[f[0]].cast<double>(), b = its.vertices[f[1]].cast<double>(), c = its.vertices[f[2]].cast<double>();
        vol += a.dot(b.cross(c));
    }
    return vol / 6.;
}

inline void make_outward(indexed_triangle_set &its)
{
    if (signed_volume(its) < 0.)
        for (auto &f : its.indices)
            std::swap(f[1], f[2]);
}

// Point-in-closed-mesh by ray parity (+Z with a tiny tilt to dodge edges). Brute force with a
// bounding-box prefilter: only used on the small detached parts.
inline bool point_inside(const indexed_triangle_set &its, const V3 &p)
{
    const V3 dir = V3(1e-4, 2e-4, 1.).normalized();
    int      hits = 0;
    for (const auto &f : its.indices) {
        const V3 a = its.vertices[f[0]].cast<double>(), b = its.vertices[f[1]].cast<double>(), c = its.vertices[f[2]].cast<double>();
        if (std::max({a.z(), b.z(), c.z()}) < p.z()) continue;
        if (std::min({a.x(), b.x(), c.x()}) > p.x() + 1e-3 || std::max({a.x(), b.x(), c.x()}) < p.x() - 1e-3) continue;
        if (std::min({a.y(), b.y(), c.y()}) > p.y() + 1e-3 || std::max({a.y(), b.y(), c.y()}) < p.y() - 1e-3) continue;
        // Moeller-Trumbore
        const V3     e1 = b - a, e2 = c - a, pv = dir.cross(e2);
        const double det = e1.dot(pv);
        if (std::abs(det) < 1e-14) continue;
        const double inv = 1. / det;
        const V3     tv  = p - a;
        const double uu  = tv.dot(pv) * inv;
        if (uu < 0. || uu > 1.) continue;
        const V3     qv = tv.cross(e1);
        const double vv = dir.dot(qv) * inv;
        if (vv < 0. || uu + vv > 1.) continue;
        if (e2.dot(qv) * inv > 0.)
            ++hits;
    }
    return (hits & 1) == 1;
}

// Point-in-polygon (even-odd) for a closed 2D polygon.
inline bool point_in_polygon(const std::vector<V2> &poly, const V2 &p)
{
    bool in = false;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const V2 &a = poly[i], &b = poly[j];
        if ((a.y() > p.y()) != (b.y() > p.y()) && p.x() < (b.x() - a.x()) * (p.y() - a.y()) / (b.y() - a.y()) + a.x())
            in = !in;
    }
    return in;
}

// ------------------------------------------------------------------------------------------------ view path
struct ViewCutter
{
    Frame               frame;
    std::vector<V2>     path;      // smoothed + extended path in the view plane
    std::vector<V2>     left_side; // polygon covering everything left of the path
    indexed_triangle_set mesh;
    std::string         error;
    bool valid() const { return error.empty() && !mesh.indices.empty(); }

    // +1 left of the path, -1 right of it.
    int side(const V3 &p) const { return point_in_polygon(left_side, frame.to2d(p)) ? 1 : -1; }
};

// points: picked points in the same coordinate system as the mesh to be cut.
// bb_min/bb_max: bounding box of that mesh. kerf: wall thickness = clearance between parts.
inline ViewCutter make_view_cutter(const std::vector<V3> &points, const V3 &view_dir, const V3 &bb_min, const V3 &bb_max, double kerf)
{
    ViewCutter vc;
    const double diag = (bb_max - bb_min).norm();
    vc.frame          = frame_from_normal(view_dir, (bb_min + bb_max) * 0.5);

    std::vector<V2> q;
    for (const V3 &p : points)
        q.push_back(vc.frame.to2d(p));
    q = remove_close_points(q, 1e-3, false);
    if (q.size() < 2) {
        vc.error = "The path needs at least two points.";
        return vc;
    }
    if (q.size() >= 3)
        q = remove_close_points(catmull_rom(q, false, 16), 1e-4, false);

    // Extend both ends well outside the object so the wall cuts through it completely.
    const V2 t0 = (q[0] - q[1]).normalized(), t1 = (q.back() - q[q.size() - 2]).normalized();
    q.insert(q.begin(), q.front() + t0 * diag);
    q.push_back(q.back() + t1 * diag);
    vc.path = q;

    const V2 chord = q.back() - q.front();
    if (chord.norm() < 1e-6) {
        vc.error = "The path must not end where it starts. Use the loop mode for closed paths.";
        return vc;
    }
    const V2 perp = V2(-chord.y(), chord.x()).normalized();
    vc.left_side  = q;
    vc.left_side.push_back(q.back() + perp * 3. * diag);
    vc.left_side.push_back(q.front() + perp * 3. * diag);

    // Offset polyline by +-kerf/2 (note: no variables named L / R, `L` is a translation macro in Orca) with mitred joints (miter clamped).
    const size_t    n = q.size();
    std::vector<V2> off_l(n), off_r(n);
    for (size_t i = 0; i < n; ++i) {
        const V2 tin  = i > 0 ? V2((q[i] - q[i - 1]).normalized()) : V2((q[1] - q[0]).normalized());
        const V2 tout = i + 1 < n ? V2((q[i + 1] - q[i]).normalized()) : tin;
        V2       t    = tin + tout;
        if (t.norm() < 1e-9) t = tout;
        t.normalize();
        const V2     m     = V2(-t.y(), t.x());
        const double cos_h = std::max(0.5, m.dot(V2(-tout.y(), tout.x())));
        off_l[i] = q[i] + m * (0.5 * kerf / cos_h);
        off_r[i] = q[i] - m * (0.5 * kerf / cos_h);
    }

    // Prism: bottom (h = -diag) and top (h = +diag) copies of the strip L/R.
    indexed_triangle_set &its = vc.mesh;
    std::vector<int>      Lb(n), Lt(n), Rb(n), Rt(n);
    for (size_t i = 0; i < n; ++i) {
        Lb[i] = add_vertex(its, vc.frame.to3d(off_l[i], -diag));
        Lt[i] = add_vertex(its, vc.frame.to3d(off_l[i], diag));
        Rb[i] = add_vertex(its, vc.frame.to3d(off_r[i], -diag));
        Rt[i] = add_vertex(its, vc.frame.to3d(off_r[i], diag));
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        add_quad(its, Lb[i], Rb[i], Rb[i + 1], Lb[i + 1]); // bottom cap
        add_quad(its, Lt[i], Lt[i + 1], Rt[i + 1], Rt[i]); // top cap
        add_quad(its, Lb[i], Lb[i + 1], Lt[i + 1], Lt[i]); // left wall
        add_quad(its, Rb[i], Rt[i], Rt[i + 1], Rb[i + 1]); // right wall
    }
    add_quad(its, Lb[0], Lt[0], Rt[0], Rb[0]);                         // start cap
    add_quad(its, Lb[n - 1], Rb[n - 1], Rt[n - 1], Lt[n - 1]);         // end cap
    make_outward(its);
    return vc;
}

// ------------------------------------------------------------------------------------------------ closed loop
struct LoopCutter
{
    Frame               frame;
    V2                  center2d{V2::Zero()};
    std::vector<double> angles; // spoke angles around center2d, increasing
    std::vector<double> radius; // loop radius per spoke
    std::vector<double> hloop;  // loop height above the mean plane per spoke
    double              hmean{0.};
    double              margin{0.};
    double              kerf{0.};
    std::vector<V3>     loop;   // dense loop (for display)
    indexed_triangle_set mesh;
    std::string         error;
    bool valid() const { return error.empty() && !mesh.indices.empty(); }

    // Interpolate the per-spoke data at an angle.
    void spoke_at(double ang, double &r, double &h) const
    {
        const size_t n = angles.size();
        ang            = angles[0] + std::fmod(std::fmod(ang - angles[0], 2. * PI_) + 2. * PI_, 2. * PI_);
        size_t i       = size_t(std::upper_bound(angles.begin(), angles.end(), ang) - angles.begin());
        size_t i0      = (i + n - 1) % n, i1 = i % n;
        double a0 = angles[i0], a1 = angles[i1];
        if (a1 <= a0) a1 += 2. * PI_;
        if (ang < a0) ang += 2. * PI_;
        const double t = a1 - a0 > 1e-12 ? (ang - a0) / (a1 - a0) : 0.;
        r = radius[i0] * (1. - t) + radius[i1] * t;
        h = hloop[i0] * (1. - t) + hloop[i1] * t;
    }
    // Height of the cap mid-surface above the mean plane at a 2D point; also returns whether the
    // point lies inside the cap footprint (loop + margin).
    double cap_height(const V2 &p, bool *inside = nullptr) const
    {
        const V2     d = p - center2d;
        const double rr = d.norm();
        double       r, h;
        spoke_at(std::atan2(d.y(), d.x()), r, h);
        if (inside) *inside = rr <= r + margin;
        const double s = r > 1e-12 ? std::min(1., rr / r) : 1.;
        return hmean * (1. - s) + h * s;
    }
    // +1 / -1 for points next to the cap on either side, 0 when far from it.
    int side(const V3 &p) const
    {
        bool         inside = false;
        const double dz     = frame.height(p) - cap_height(frame.to2d(p), &inside);
        if (!inside || std::abs(dz) > 2. * kerf + 1e-3)
            return 0;
        return dz > 0. ? 1 : -1;
    }
};

// loop_points: points on the surface around the region to detach (closed, >= 3 points).
// They are smoothed; pass a `snap` functor (V3 -> V3) to put the smoothed curve back on the surface.
template<class Snap>
LoopCutter make_loop_cutter(const std::vector<V3> &loop_points, double kerf, double margin, Snap snap, int rings = 10)
{
    LoopCutter lc;
    lc.kerf   = kerf;
    lc.margin = margin;
    std::vector<V3> pts = remove_close_points(loop_points, 1e-3, true);
    if (pts.size() < 3) {
        lc.error = "The loop needs at least three points.";
        return lc;
    }
    std::vector<V3> dense = catmull_rom(pts, true, 16);
    for (V3 &p : dense)
        p = snap(p);
    dense = remove_close_points(dense, 1e-4, true);
    lc.loop = dense;

    // Mean plane of the loop.
    V3 c = V3::Zero();
    for (const V3 &p : dense) c += p;
    c /= double(dense.size());
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const V3 &p : dense) { const Eigen::Vector3d d = p - c; cov += d * d.transpose(); }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    lc.frame = frame_from_normal(V3(es.eigenvectors().col(0)), c);

    std::vector<V2>     q;
    std::vector<double> h;
    for (const V3 &p : dense) { q.push_back(lc.frame.to2d(p)); h.push_back(lc.frame.height(p)); }

    // Polygon centroid in the plane, then check the loop is star-shaped around it.
    double A = 0.; V2 cc = V2::Zero();
    for (size_t i = 0; i < q.size(); ++i) {
        const V2 &a = q[i], &b = q[(i + 1) % q.size()];
        const double cr = a.x() * b.y() - b.x() * a.y();
        A += cr; cc += (a + b) * cr;
    }
    if (std::abs(A) < 1e-9) {
        lc.error = "The loop is degenerate (it encloses no area).";
        return lc;
    }
    cc /= (3. * A);
    if (A < 0.) { std::reverse(q.begin(), q.end()); std::reverse(h.begin(), h.end()); }
    lc.center2d = cc;

    const size_t n = q.size();
    std::vector<double> ang(n), rad(n);
    for (size_t i = 0; i < n; ++i) {
        const V2 d = q[i] - cc;
        ang[i] = std::atan2(d.y(), d.x());
        rad[i] = d.norm();
    }
    // Start at the smallest angle and require strictly increasing angles (star-shaped loop).
    const size_t i0 = size_t(std::min_element(ang.begin(), ang.end()) - ang.begin());
    std::rotate(ang.begin(), ang.begin() + i0, ang.end());
    std::rotate(rad.begin(), rad.begin() + i0, rad.end());
    std::rotate(h.begin(), h.begin() + i0, h.end());
    std::rotate(q.begin(), q.begin() + i0, q.end());
    for (size_t i = 1; i < n; ++i)
        if (ang[i] <= ang[i - 1] || rad[i] < 1e-6) {
            lc.error = "The loop folds back on itself seen from its mean plane. Place the points more evenly around the part.";
            return lc;
        }
    lc.angles = ang; lc.radius = rad; lc.hloop = h;
    double hm = 0.; for (double x : h) hm += x;
    lc.hmean = hm / double(n);

    // Cap: center + `rings` rings up to the loop + one outer ring (loop pushed out by margin),
    // each vertex duplicated at +-kerf/2 along the plane normal.
    indexed_triangle_set &its = lc.mesh;
    auto pos2d = [&](size_t j, int k) -> V2 {
        const V2 d = (q[j] - cc) / rad[j];
        return k <= rings ? V2(cc + (q[j] - cc) * (double(k) / rings)) : V2(q[j] + d * margin);
    };
    auto hgt = [&](size_t j, int k) { const double s = std::min(1., double(k) / rings); return lc.hmean * (1. - s) + h[j] * s; };
    const int ctop = add_vertex(its, lc.frame.to3d(cc, lc.hmean + kerf / 2));
    const int cbot = add_vertex(its, lc.frame.to3d(cc, lc.hmean - kerf / 2));
    const int K    = rings + 1;
    std::vector<int> top(n * K), bot(n * K);
    for (int k = 1; k <= K; ++k)
        for (size_t j = 0; j < n; ++j) {
            top[(k - 1) * n + j] = add_vertex(its, lc.frame.to3d(pos2d(j, k), hgt(j, k) + kerf / 2));
            bot[(k - 1) * n + j] = add_vertex(its, lc.frame.to3d(pos2d(j, k), hgt(j, k) - kerf / 2));
        }
    for (size_t j = 0; j < n; ++j) {
        const size_t j1 = (j + 1) % n;
        add_tri(its, ctop, top[j], top[j1]);
        add_tri(its, cbot, bot[j1], bot[j]);
        for (int k = 1; k < K; ++k) {
            const size_t a = (k - 1) * n, b = k * n;
            add_quad(its, top[a + j], top[b + j], top[b + j1], top[a + j1]);
            add_quad(its, bot[a + j], bot[a + j1], bot[b + j1], bot[b + j]);
        }
        const size_t o = (K - 1) * n;
        add_quad(its, top[o + j], bot[o + j], bot[o + j1], top[o + j1]);
    }
    make_outward(its);
    return lc;
}

// ------------------------------------------------------------------------------------------------ cleanup
// Booleans (MCUT) reject meshes with non-manifold edges. Typical culprits in sculpted / AI meshes:
// zero-thickness sheets (the same triangle twice with opposite winding) and surfaces pinched along an
// edge. Removes duplicate faces and detaches the extra faces of an edge shared by more than two
// faces by giving them their own copies of the edge vertices. Returns the number of fixes.
inline size_t make_edge_manifold(indexed_triangle_set &its)
{
    size_t fixes = 0;
    // 1) duplicate faces: same vertex set. Opposite winding -> zero-volume sheet, drop both.
    {
        std::map<std::array<int, 3>, std::vector<size_t>> by_set;
        for (size_t i = 0; i < its.indices.size(); ++i) {
            std::array<int, 3> k{its.indices[i][0], its.indices[i][1], its.indices[i][2]};
            std::sort(k.begin(), k.end());
            by_set[k].push_back(i);
        }
        std::vector<bool> drop(its.indices.size(), false);
        for (auto &[k, fs] : by_set) {
            if (fs.size() < 2) continue;
            for (size_t j = 1; j < fs.size(); ++j) drop[fs[j]] = true;
            auto same_winding = [&](size_t a, size_t b) {
                const auto &A = its.indices[a], &B = its.indices[b];
                for (int r = 0; r < 3; ++r)
                    if (A[0] == B[r] && A[1] == B[(r + 1) % 3] && A[2] == B[(r + 2) % 3]) return true;
                return false;
            };
            if (!same_winding(fs[0], fs[1])) drop[fs[0]] = true;
            ++fixes;
        }
        size_t w = 0;
        for (size_t i = 0; i < its.indices.size(); ++i)
            if (!drop[i]) its.indices[w++] = its.indices[i];
        its.indices.resize(w);
    }
    // 2) edges shared by more than two faces: keep one (a,b)/(b,a) pair, detach the rest.
    std::map<std::pair<int, int>, std::vector<size_t>> by_edge;
    for (size_t i = 0; i < its.indices.size(); ++i)
        for (int e = 0; e < 3; ++e) {
            int a = its.indices[i][e], b = its.indices[i][(e + 1) % 3];
            by_edge[{std::min(a, b), std::max(a, b)}].push_back(i);
        }
    for (auto &[edge, fs] : by_edge) {
        if (fs.size() <= 2) continue;
        ++fixes;
        auto dir = [&](size_t f) { // +1 if face walks edge as (min,max)
            for (int e = 0; e < 3; ++e)
                if (its.indices[f][e] == edge.first && its.indices[f][(e + 1) % 3] == edge.second) return 1;
            return -1;
        };
        size_t keep0 = fs[0], keep1 = fs[0];
        for (size_t f : fs) if (f != keep0 && dir(f) != dir(keep0)) { keep1 = f; break; }
        for (size_t f : fs) {
            if (f == keep0 || f == keep1) continue;
            for (int e = 0; e < 3; ++e) {
                int &vi = its.indices[f][e];
                if (vi == edge.first || vi == edge.second) {
                    its.vertices.push_back(its.vertices[vi]);
                    vi = int(its.vertices.size()) - 1;
                }
            }
        }
    }
    return fixes;
}

// ------------------------------------------------------------------------------------------------ grouping
struct Grouped
{
    indexed_triangle_set side_a, side_b; // side_a = +1 side
    size_t               n_a{0}, n_b{0}, n_small{0};
};

inline void its_append(indexed_triangle_set &dst, const indexed_triangle_set &src)
{
    const int off = int(dst.vertices.size());
    dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
    for (auto f : src.indices) dst.indices.emplace_back(f[0] + off, f[1] + off, f[2] + off);
}

// side_of(V3) -> +1 / -1 / 0. Components get the majority vote of their vertices (0 = no vote).
// Components with no vote (untouched by the cutter) go to side_b, unless `untouched_follow_inside`
// and one of their vertices lies inside a component already assigned to side_a
// (e.g. a pin overlapping the detached foot).
template<class SideOf>
Grouped group_components(const std::vector<indexed_triangle_set> &parts, SideOf side_of, bool untouched_follow_inside, double small_volume = 1.)
{
    Grouped          g;
    std::vector<int> label(parts.size(), 0);
    for (size_t i = 0; i < parts.size(); ++i) {
        const auto &vs   = parts[i].vertices;
        const size_t step = std::max<size_t>(1, vs.size() / 4000);
        long         vote = 0;
        for (size_t k = 0; k < vs.size(); k += step) vote += side_of(V3(vs[k].cast<double>()));
        label[i] = vote > 0 ? 1 : vote < 0 ? -1 : 0;
    }
    if (untouched_follow_inside)
        for (size_t i = 0; i < parts.size(); ++i) {
            if (label[i] != 0) continue;
            const auto  &vs   = parts[i].vertices;
            const size_t step = std::max<size_t>(1, vs.size() / 64);
            for (size_t j = 0; j < parts.size() && label[i] == 0; ++j)
                if (label[j] == 1)
                    for (size_t k = 0; k < vs.size(); k += step)
                        if (point_inside(parts[j], V3(vs[k].cast<double>()))) { label[i] = 1; break; }
        }
    for (size_t i = 0; i < parts.size(); ++i) {
        if (std::abs(signed_volume(parts[i])) < small_volume) ++g.n_small;
        if (label[i] == 1) { its_append(g.side_a, parts[i]); ++g.n_a; }
        else               { its_append(g.side_b, parts[i]); ++g.n_b; }
    }
    return g;
}

}} // namespace Slic3r::PathCut
