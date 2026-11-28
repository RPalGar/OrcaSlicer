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

// The surface cut_with_slab() removes the slab around: a signed field whose zero set is the middle of the
// slab, and a 2D chart of the slab's faces to triangulate the caps in.
struct Cutter
{
    double      kerf{0.}; // slab thickness = clearance between the parts
    std::string error;    // translated; set when the cutter could not be built

    virtual ~Cutter() = default;
    // Whether the slab has a rim inside the object (see outside()).
    virtual bool   limited() const = 0;
    // Signed distance across the slab, positive on the side of level +1 (may be clamped far from it).
    virtual double field(const V3 &p) const = 0;
    // How far a point lies outside a limited slab, negative inside it.
    virtual double outside(const V3 &) const { return -1.; }
    // Which side of the cut a point away from the slab is on: +1, -1, or 0 when that cannot be told.
    virtual int    side(const V3 &p) const = 0;
    // Chart of the slab faces, and back to the face of the given level (+1 / -1).
    virtual V2     chart(const V3 &p) const = 0;
    virtual V3     lift(const V2 &chart_pos, int level) const = 0;
    // Chart lines where the slab bends by more than `tol`; between two of them it is planar.
    virtual std::vector<std::pair<V2, V2>> kinks(double /*tol*/) const { return {}; }
    // Message when the surface runs through the rim of a limited slab.
    virtual std::string rim_error() const = 0;
    // Whether parts the limited slab ends inside are left uncut (see cut_with_slab) instead of refused.
    virtual bool skips_parts_it_ends_in() const { return false; }
};

// The wall swept from a path in the view plane (built by make_view_cutter).
struct ViewCutter final : Cutter
{
    Frame               frame;
    std::vector<V2>     path; // smoothed + extended path in the view plane
    std::vector<double> arc;  // arc length at each path point
    // Limits (see ViewLimits): the wall may stop short past the first / last picked point and span only
    // heights [h_min, h_max] along the sweep direction.
    bool                limit_start{false}, limit_end{false}, limit_depth{false};
    double              h_min{0.}, h_max{0.};
    bool valid() const { return error.empty() && path.size() >= 2; }
    bool limited() const override { return limit_start || limit_end || limit_depth; }

    // +1 left of the path, -1 right of it.
    int    side(const V3 &p) const override { return exact_field(frame.to2d(p)) >= 0. ? 1 : -1; }
    // Distance across the wall, positive on the left of the path.
    double field(const V3 &p) const override;
    double outside(const V3 &p) const override;
    // (arc length along the path, height along the view direction)
    V2     chart(const V3 &p) const override;
    V3     lift(const V2 &chart_pos, int level) const override;
    std::vector<std::pair<V2, V2>> kinks(double tol) const override;
    std::string rim_error() const override;
    bool        skips_parts_it_ends_in() const override { return true; }

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

// The picked points joined by the same smooth curve the cut follows (for drawing the path).
std::vector<V3> smooth_path(const std::vector<V3> &points);

// A cap spanning a loop picked around the region to detach (built by make_loop_cutter): a smooth height field
// over the loop's mean plane, from the loop's mean height at its centre to the loop itself at the loop, flat
// beyond. The cap only reaches its rim, a little past the section of the part around the loop (see
// fit_loop_margin), so only that part is cut.
struct LoopCutter final : Cutter
{
    static constexpr int spokes = 360; // profiles below are sampled at angles 2 pi k / spokes around center

    Frame               frame;
    V2                  center{V2::Zero()};
    std::vector<double> radius; // loop radius per spoke
    std::vector<double> hloop;  // loop height above the mean plane per spoke
    std::vector<double> rim;    // how far the cap reaches per spoke
    double              hmean{0.};
    double              margin{0.}; // requested reach past the loop
    std::vector<V3>     loop;       // the smoothed loop (for drawing it)
    std::vector<V3>     stray;      // when fit_loop_margin() refuses the loop, shown to the user: segments (pairs of
                                    // points) of the section it does not go around, or of the other part in the way
    bool                stray_in_the_way{false}; // stray is another part the cap would run into
    bool valid() const { return error.empty() && radius.size() == size_t(spokes); }
    bool limited() const override { return true; }

    // Height of the cap above the mean plane at a point of the plane; `outside` is how far the point lies
    // outside the cap's rim, negative inside.
    double cap_height(const V2 &q, double *outside = nullptr) const;
    // How much farther a point is above the cap than its distance to it, where the cap slopes: the clearance
    // is measured across the cap, so that it is the same where the cap slopes as where it is level.
    double slope_scale(const V2 &q) const;
    // Distance above the cap, across it.
    double field(const V3 &p) const override;
    double outside(const V3 &p) const override { double o; cap_height(frame.to2d(p), &o); return o; }
    // +1 / -1 next to the cap above / below it, 0 elsewhere.
    int    side(const V3 &p) const override;
    // The cap seen along its normal.
    V2     chart(const V3 &p) const override { return frame.to2d(p); }
    V3     lift(const V2 &q, int level) const override;
    std::string rim_error() const override;
};

// loop_points: points on the surface around the region to detach (closed, at least 3). The smoothed loop is
// put back on the surface with `snap`. flat: the cap is the loop's mean plane, for a flat cut face, instead
// of following the loop's heights. The rim is set `margin` past the loop; fit_loop_margin() fits it to the
// model.
LoopCutter make_loop_cutter(const std::vector<V3> &loop_points, double kerf, double margin, const std::function<V3(const V3 &)> &snap,
                            bool flat);

// Fits the cap's rim to the model: the section of the part the loop goes around (where the surface crosses
// the cap around its centre) often bulges a little past the smoothed loop (ribs, folds, fingers), and other
// parts of the model may come close (a hand resting on a hip). The section is the curve around the cap's
// centre that the loop follows; other curves around it (a fist around a grip) belong to other parts. Per
// spoke, the rim is put `margin` past that section, or halfway to another part when that is nearer. Fails
// with `error` set when the section reaches far past the loop (the loop does not go around the part; `stray`
// then holds the section) or another part leaves no room for the rim.
bool fit_loop_margin(const indexed_triangle_set &its, LoopCutter &lc);

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
CutResult cut_with_slab(indexed_triangle_set its, const Cutter &cutter, double refine = 0.5,
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
Grouped group_cut(const CutResult &cut, const Cutter &cutter, bool *joined = nullptr);

// Where the parts a cut leaves connected around it (see group_cut) still hold together: per such part, the
// shortest line over its surface from one cut face to the other. Cutting across it detaches the part.
std::vector<std::vector<V3>> connections(const CutResult &cut);

// Preview of a cut: how far the detached piece is drawn apart from the piece that stays, along the line
// between their bounding box centres (straight up when they share a centre): 10 % of the detached piece's
// size, between 6 and 30 mm. Zero when either piece is empty. Both meshes must have no unused vertices, as
// group_cut() returns them.
V3 preview_shift(const indexed_triangle_set &moved, const indexed_triangle_set &stay);
// How far both pieces are drawn raised so that the detached piece, drawn `shift` apart, does not go below
// the lowest point of the object. Zero otherwise.
double preview_lift(const indexed_triangle_set &moved, const indexed_triangle_set &stay, const V3 &shift);

}} // namespace Slic3r::PathCut

#endif // slic3r_PathCut_hpp_
