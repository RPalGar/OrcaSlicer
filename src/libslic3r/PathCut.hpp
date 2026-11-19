#ifndef slic3r_PathCut_hpp_
#define slic3r_PathCut_hpp_

// Cut along a path picked from the view: an open polyline picked on the object is smoothed, extended past
// it and swept along the view direction into a thin wall. Removing the wall splits the object; its
// thickness ("kerf") is the clearance left between the parts.
//
// The cut does not use a general mesh boolean: those are far too slow on the dense meshes this tool is meant
// for (sculpted and scanned figures, often millions of triangles). The wall is described by a signed field
// (distance across it) and is |field| < kerf / 2. Triangles near the wall are refined, clipped against the
// two level sets field = +-kerf / 2 (crossings are computed once per edge, so the result stays watertight)
// and the openings are closed with caps triangulated in a 2D chart of each level set. The work is linear in
// the triangle count, and only the triangles near the wall are touched.
// See docs/HLSD/path-cut.md.

#include <admesh/stl.h>
#include <Eigen/Core>

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r { namespace PathCut {

using V2 = Eigen::Matrix<double, 2, 1, Eigen::DontAlign>;
using V3 = Eigen::Matrix<double, 3, 1, Eigen::DontAlign>;

struct Frame
{
    V3 c{V3::Zero()}, u{V3::UnitX()}, v{V3::UnitY()}, n{V3::UnitZ()};
    V2     to2d(const V3 &p) const { const V3 d = p - c; return {d.dot(u), d.dot(v)}; }
    double height(const V3 &p) const { return (p - c).dot(n); }
    V3     to3d(const V2 &q, double h) const { return c + q.x() * u + q.y() * v + h * n; }
};

// How far the wall reaches, in mm; negative = through the whole object.
struct ViewLimits
{
    double start{-1.}, end{-1.};  // past the first / last picked point, along the path
    double front{-1.}, back{-1.}; // in front of / behind the picked points, along the sweep direction
};

// The wall swept from a path in the view plane (built by make_view_cutter).
struct ViewCutter
{
    Frame               frame;
    std::vector<V2>     path; // smoothed + extended path in the view plane
    std::vector<double> arc;  // arc length at each path point
    double              kerf{0.};
    // Limits (see ViewLimits): the wall may stop short past the first / last picked point and span only
    // heights [h_min, h_max] along the sweep direction.
    bool                limit_start{false}, limit_end{false}, limit_depth{false};
    double              h_min{0.}, h_max{0.};
    std::string         error; // translated
    bool valid() const { return error.empty() && path.size() >= 2; }
    bool limited() const { return limit_start || limit_end || limit_depth; }

    // +1 left of the path, -1 right of it.
    int    side(const V3 &p) const { return exact_field(frame.to2d(p)) >= 0. ? 1 : -1; }
    // Signed distance across the wall, positive on the left of the path (clamped far from it).
    double field(const V3 &p) const;
    // How far a point lies outside the limited wall (negative inside it, always inside when unlimited).
    double outside(const V3 &p) const;
    // Chart of the wall faces: (arc length along the path, height along the view direction), and back
    // to the face of the given level (+1 / -1).
    V2     chart(const V3 &p) const;
    V3     lift(const V2 &chart_pos, int level) const;
    // Chart lines where the wall bends by more than `tol`; between two of them it is planar.
    std::vector<std::pair<V2, V2>> kinks(double tol) const;

    // Acceleration grid over the view plane: segments near each cell, and for cells with none, which
    // side of the path the whole cell is on.
    V2                            grid_min{V2::Zero()};
    double                        cell{1.};
    double                        reach{5.}; // segments closer than this to a cell are listed in it
    int                           nx{0}, ny{0};
    std::vector<std::vector<int>> cell_segs;
    std::vector<signed char>      cell_side;
    void build_grid(const V2 &mn, const V2 &mx);

private:
    V2     seg_normal(size_t i) const;
    V2     vertex_normal(size_t k) const;
    double exact_field(const V2 &q) const;
    double signed_from(const V2 &q, size_t seg, double t, double dist) const;
    void   nearest(const V2 &q, size_t &seg, double &t, double &dist) const;
    bool   nearest_in_grid(const V2 &q, size_t &seg, double &t, double &dist) const;
};

// points: picked points in the coordinate system of the mesh to be cut; view_dir: the sweep direction.
// bb_min/bb_max: bounding box of that mesh. kerf: wall thickness = clearance between the parts.
// On failure the cutter is not valid() and error says why.
ViewCutter make_view_cutter(const std::vector<V3> &points, const V3 &view_dir, const V3 &bb_min, const V3 &bb_max, double kerf,
                            const ViewLimits &limits = {});

// First crossing of two non-adjacent segments of an open polyline, searching from the start: the segments
// i < j and the crossing point. Returns false when the polyline does not cross itself.
bool first_self_crossing(const std::vector<V2> &q, size_t &si, size_t &sj, V2 &x);

// Sculpted / AI meshes often carry zero-thickness sheets (the same triangle twice with opposite winding)
// and surfaces pinched along an edge, which the cut cannot close. Removes duplicate faces and gives the
// extra faces of an edge shared by more than two faces their own copies of the edge vertices.
// Returns the number of fixes.
size_t make_edge_manifold(indexed_triangle_set &its);

struct CutResult
{
    indexed_triangle_set     mesh;  // both sides, closed; split it into components to get the parts
    std::string              error; // translated
    size_t                   first_cap{0}; // index of the first cap triangle in mesh.indices
    std::vector<signed char> level;        // per vertex: +1 / -1 on the cut face of that side, 0 elsewhere
};

// Removes the wall from a closed mesh and closes the openings.
//  refine:   edges longer than this near the wall are split before clipping [mm].
//  canceled: polled now and then; when it returns true the cut stops with an empty result.
CutResult cut_with_slab(indexed_triangle_set its, const ViewCutter &cutter, double refine = 0.5,
                        const std::function<bool()> &canceled = {});

struct Grouped
{
    indexed_triangle_set side_a, side_b; // side_a = +1 side
    size_t               n_a{0}, n_b{0}, n_small{0};
};

// Splits a cut result into its connected parts and assigns each to a side by the cut face it carries.
// Parts without a cut face go to the side of the cutter most of their vertices are on. A part carrying both
// cut faces is still connected around the cut (an arm resting on a leg, say): it goes to side b and
// `joined` is set.
Grouped group_cut(const CutResult &cut, const ViewCutter &cutter, bool *joined = nullptr);

// Preview of a cut: how far the detached piece is drawn apart from the piece that stays, along the line
// between their bounding box centres (straight up when they share a centre): 10 % of the detached piece's
// size, between 6 and 30 mm. Zero when either piece is empty.
V3 preview_shift(const indexed_triangle_set &moved, const indexed_triangle_set &stay);
// How far both pieces are drawn raised so that the detached piece, drawn `shift` apart, does not go below
// the lowest point of the object. Zero otherwise.
double preview_lift(const indexed_triangle_set &moved, const indexed_triangle_set &stay, const V3 &shift);

}} // namespace Slic3r::PathCut

#endif // slic3r_PathCut_hpp_
