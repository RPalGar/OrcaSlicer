#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <libslic3r/AABBMesh.hpp>
#include <libslic3r/PathCut.hpp>
#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/libslic3r.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

namespace {

// Number of edges used by exactly one face (0 for a closed mesh) and the largest number of faces on one edge.
void edge_counts(const indexed_triangle_set &its, size_t &open, size_t &max_faces)
{
    std::map<std::pair<int, int>, size_t> count;
    for (const auto &f : its.indices)
        for (int e = 0; e < 3; ++e) {
            const int a = f[e], b = f[(e + 1) % 3];
            ++count[{std::min(a, b), std::max(a, b)}];
        }
    open = max_faces = 0;
    for (const auto &[edge, n] : count) {
        if (n == 1) ++open;
        max_faces = std::max(max_faces, n);
    }
}

void require_closed_parts(const std::vector<indexed_triangle_set> &parts)
{
    for (const auto &p : parts) {
        size_t open, max_faces;
        edge_counts(p, open, max_faces);
        CHECK(open == 0);
        CHECK(max_faces == 2);
        CHECK(its_volume(p) > 0.);
    }
}

} // namespace
TEST_CASE("A wavy path picked from the front splits a cube into two closed parts", "[PathCut]")
{
    // 20 mm cube, path drawn looking along +Y, zig-zagging across X at mid height.
    const indexed_triangle_set cube = its_make_cube(20., 20., 20.);
    const double               kerf = 0.2;
    const std::vector<PathCut::V3> pts{{2., 0., 8.}, {7., 0., 13.}, {13., 0., 7.}, {18., 0., 12.}};

    const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3::Constant(20.), kerf);
    REQUIRE(vc.valid());

    const PathCut::CutResult res = PathCut::cut_with_slab(cube, vc);
    REQUIRE(res.error.empty());
    const auto parts = its_split(res.mesh);
    REQUIRE(parts.size() == 2);
    require_closed_parts(parts);

    const PathCut::Grouped g = PathCut::group_cut(res, vc);
    CHECK(g.n_a == 1);
    CHECK(g.n_b == 1);
    // The two halves add up to the cube minus a wall kerf thick, 20 mm deep, as long as the path inside the cube.
    const double total = its_volume(g.side_a) + its_volume(g.side_b);
    CHECK(total < 8000.);
    CHECK(total > 8000. - kerf * 20. * 40.);
}

TEST_CASE("A view path limited in depth cuts the front box and leaves the one behind it whole", "[PathCut]")
{
    // Looking along +Y: a front box (y 0..10) and a back box (y 20..30).
    indexed_triangle_set boxes = its_make_cube(20., 10., 20.), back = its_make_cube(20., 10., 20.);
    its_translate(back, Vec3f(0.f, 20.f, 0.f));
    its_merge(boxes, back);
    const std::vector<PathCut::V3> pts{{2., 0., 10.}, {10., 0., 12.}, {18., 0., 9.}};

    SECTION("15 mm behind the points stops in the gap")
    {
        PathCut::ViewLimits lim;
        lim.back = 15.;
        const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3(20., 30., 20.), 0.2, lim);
        REQUIRE(vc.limited());
        const PathCut::CutResult res = PathCut::cut_with_slab(boxes, vc);
        REQUIRE(res.error.empty());
        const auto parts = its_split(res.mesh);
        CHECK(parts.size() == 3);
        require_closed_parts(parts);
        bool joined = true;
        const PathCut::Grouped g = PathCut::group_cut(res, vc, &joined);
        CHECK_FALSE(joined);
        CHECK(g.n_a + g.n_b == 3);
    }
    SECTION("25 mm behind the points ends inside the back box: the front box is cut, the back box is left whole")
    {
        PathCut::ViewLimits lim;
        lim.back = 25.;
        const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3(20., 30., 20.), 0.2, lim);
        const PathCut::CutResult res = PathCut::cut_with_slab(boxes, vc);
        REQUIRE(res.error.empty());
        const auto parts = its_split(res.mesh);
        CHECK(parts.size() == 3);
        require_closed_parts(parts);
        const bool back_whole = std::any_of(parts.begin(), parts.end(), [](const indexed_triangle_set &p) { return std::abs(its_volume(p) - 4000.) < 1.; });
        CHECK(back_whole);
    }
    SECTION("5 mm behind the points ends inside the front box, the only one it reaches, and is reported")
    {
        PathCut::ViewLimits lim;
        lim.back = 5.;
        const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3(20., 30., 20.), 0.2, lim);
        CHECK_FALSE(PathCut::cut_with_slab(boxes, vc).error.empty());
    }
}

TEST_CASE("The preview draws the detached piece apart from the rest", "[PathCut]")
{
    // A 10 mm cube resting on top of a 20 mm one: drawn 3 mm (the minimum) further up.
    indexed_triangle_set stay  = its_make_cube(20., 20., 20.);
    indexed_triangle_set moved = its_make_cube(10., 10., 10.);
    for (Vec3f &v : moved.vertices) v += Vec3f(5.f, 5.f, 20.f);
    const PathCut::V3 shift = PathCut::preview_shift(moved, stay);
    CHECK_THAT(shift.norm(), WithinAbs(6., 1e-9));
    CHECK(shift.z() > 0.);
    CHECK_THAT(shift.x(), WithinAbs(0., 1e-9));
    // A large detached piece goes 10 % of its size, at most 30 mm; an empty piece is not moved.
    indexed_triangle_set big = its_make_cube(400., 400., 400.);
    for (Vec3f &v : big.vertices) v += Vec3f(0.f, 0.f, 20.f);
    CHECK_THAT(PathCut::preview_shift(big, stay).norm(), WithinAbs(30., 1e-9));
    CHECK(PathCut::preview_shift(indexed_triangle_set{}, stay).norm() == 0.);
    // Drawn above the other piece, nothing has to be raised.
    CHECK_THAT(PathCut::preview_lift(moved, stay, shift), WithinAbs(0., 1e-9));
}

TEST_CASE("The preview never draws a piece below the object", "[PathCut]")
{
    // The detached piece is the 10 mm cube under a 20 mm one: drawn 6 mm further down, it would sink
    // below the bed. Both pieces are raised so that it rests where the object rested, the other above it.
    indexed_triangle_set stay  = its_make_cube(20., 20., 20.);
    indexed_triangle_set moved = its_make_cube(10., 10., 10.);
    for (Vec3f &v : stay.vertices) v += Vec3f(0.f, 0.f, 10.f);
    for (Vec3f &v : moved.vertices) v += Vec3f(5.f, 5.f, 0.f);
    const PathCut::V3 shift = PathCut::preview_shift(moved, stay);
    REQUIRE(shift.z() < 0.);
    const double lift = PathCut::preview_lift(moved, stay, shift);
    CHECK_THAT(lift, WithinAbs(-shift.z(), 1e-9));
    CHECK_THAT(bounding_box(moved).min.z() + shift.z() + lift, WithinAbs(0., 1e-6));
    CHECK_THAT(PathCut::preview_lift(indexed_triangle_set{}, stay, shift), WithinAbs(0., 1e-9));
}

TEST_CASE("Two cubes touching along an edge are separated without opening either", "[PathCut]")
{
    indexed_triangle_set a = its_make_cube(10., 10., 10.), b = its_make_cube(10., 10., 10.);
    its_translate(b, Vec3f(10.f, 0.f, 10.f));
    its_merge(a, b);
    its_merge_vertices(a); // the edge x = 10, z = 10 is now shared by four faces
    size_t open, max_faces;
    edge_counts(a, open, max_faces);
    REQUIRE(max_faces == 4);

    CHECK(PathCut::make_edge_manifold(a) > 0);
    edge_counts(a, open, max_faces);
    CHECK(open == 0);
    CHECK(max_faces == 2);
}

TEST_CASE("Zero-thickness sheets glued to a mesh are removed before cutting", "[PathCut]")
{
    indexed_triangle_set its = its_make_cube(10., 10., 10.);
    // The same triangle twice with opposite winding makes edges shared by four faces.
    const auto f = its.indices.front();
    its.indices.emplace_back(f[0], f[1], f[2]);
    its.indices.emplace_back(f[0], f[2], f[1]);
    size_t open, max_faces;
    edge_counts(its, open, max_faces);
    REQUIRE(max_faces > 2);

    CHECK(PathCut::make_edge_manifold(its) > 0);
    edge_counts(its, open, max_faces);
    CHECK(max_faces == 2);
}

TEST_CASE("Curls from hand jitter are cut out of a view path; a path doubling back is refused", "[PathCut]")
{
    const indexed_triangle_set cube = its_make_cube(20., 20., 20.);
    SECTION("two points picked almost on top of each other, the second slightly behind")
    {
        const std::vector<PathCut::V3> pts{{-2., 0., 10.}, {6., 0., 10.5}, {9., 0., 10.}, {8.5, 0., 10.3}, {10., 0., 9.8}, {14., 0., 9.5}, {22., 0., 10.}};
        const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3::Constant(20.), 0.2);
        REQUIRE(vc.valid());
        size_t i, j;
        PathCut::V2 x;
        CHECK_FALSE(PathCut::first_self_crossing(vc.path, i, j, x));
        const PathCut::CutResult res = PathCut::cut_with_slab(cube, vc);
        REQUIRE(res.error.empty());
        require_closed_parts(its_split(res.mesh));
    }
    SECTION("points picked back and forth")
    {
        const std::vector<PathCut::V3> pts{{-2., 0., 10.}, {15., 0., 10.}, {5., 0., 12.}, {22., 0., 8.}};
        CHECK_FALSE(PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3::Constant(20.), 0.2).valid());
    }
}

TEST_CASE("A limited view wall ending inside a coarse box is reported, not silently skipped", "[PathCut]")
{
    // A box long along the view direction (+Y), two triangles per face: the wall, limited to 10 mm of
    // depth around the points, ends inside it. Every vertex of the long faces lies outside the limits,
    // so the crossings decide, not the vertices.
    const indexed_triangle_set box = its_make_cube(20., 100., 20.);
    const std::vector<PathCut::V3> pts{{-2., 50., 10.}, {10., 50., 11.}, {22., 50., 9.}};
    PathCut::ViewLimits lim;
    lim.front = lim.back = 5.;
    const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), PathCut::V3::Zero(), PathCut::V3(20., 100., 20.), 0.2, lim);
    REQUIRE(vc.limited());
    const PathCut::CutResult res = PathCut::cut_with_slab(box, vc);
    CHECK_FALSE(res.error.empty());
}

TEST_CASE("Limited view walls on coarse oblique faces cut what they go through and report what they end inside", "[PathCut]")
{
    // Large faces oblique to the view, crossed by the wall far from their vertices.
    for (const float deg : {30.f, 45.f, 60.f}) {
        const std::vector<PathCut::V3> pts{{-80., 60., 10.}, {0., 60., 11.}, {80., 60., 9.}};
        // A box long along the view, rotated: 10 mm of depth end inside it.
        indexed_triangle_set box = its_make_cube(20., 100., 20.);
        its_translate(box, Vec3f(-10.f, -50.f, 0.f));
        its_rotate_z(box, float(deg * PI / 180.));
        its_translate(box, Vec3f(0.f, 60.f, 0.f));
        PathCut::ViewLimits narrow;
        narrow.front = narrow.back = 5.;
        const BoundingBoxf3 bb = bounding_box(box);
        const PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), bb.min, bb.max, 0.2, narrow);
        CHECK_FALSE(PathCut::cut_with_slab(box, vc).error.empty());
        // A thin plate, rotated: the depth limits just cover it.
        indexed_triangle_set plate = its_make_cube(40., 2., 20.);
        its_translate(plate, Vec3f(-20.f, -1.f, 0.f));
        its_rotate_z(plate, float(deg * PI / 180.));
        its_translate(plate, Vec3f(0.f, 60.f, 0.f));
        const BoundingBoxf3 pb = bounding_box(plate);
        PathCut::ViewLimits tight;
        tight.front = 60. - pb.min.y() + 0.1;
        tight.back  = pb.max.y() - 60. + 0.1;
        const PathCut::ViewCutter pc = PathCut::make_view_cutter(pts, PathCut::V3::UnitY(), pb.min, pb.max, 0.2, tight);
        REQUIRE(pc.limited());
        const PathCut::CutResult res = PathCut::cut_with_slab(plate, pc);
        REQUIRE(res.error.empty());
        const auto parts = its_split(res.mesh);
        CHECK(parts.size() == 2);
        require_closed_parts(parts);
    }
}

namespace {

std::vector<PathCut::V3> ring(double r, int n, const std::function<double(double)> &z)
{
    std::vector<PathCut::V3> loop;
    for (int i = 0; i < n; ++i) {
        const double a = 2. * PI * i / n;
        loop.emplace_back(r * std::cos(a), r * std::sin(a), z(a));
    }
    return loop;
}

const auto no_snap = [](const PathCut::V3 &p) { return p; };

} // namespace

TEST_CASE("A loop cap leaves other parts of the model next to the loop uncut", "[PathCut]")
{
    // Cylinder r = 5 and a 5 x 10 x 30 box 3 mm beside it, both crossing the plane of a loop around the
    // cylinder at z = 15. The cap stops in the free ring between them: only the cylinder is cut.
    indexed_triangle_set model = its_make_cylinder(5., 30.);
    indexed_triangle_set box   = its_make_cube(5., 10., 30.);
    its_translate(box, Vec3f(8.f, -5.f, 0.f));
    its_merge(model, box);
    for (const bool flat : {true, false}) {
        PathCut::LoopCutter lc = PathCut::make_loop_cutter(ring(5., 12, [](double) { return 15.; }), 0.2, 1.5, no_snap, flat);
        REQUIRE(lc.valid());
        REQUIRE(PathCut::fit_loop_margin(model, lc));
        const PathCut::CutResult res = PathCut::cut_with_slab(model, lc);
        REQUIRE(res.error.empty());
        const auto parts = its_split(res.mesh);
        REQUIRE(parts.size() == 3);
        require_closed_parts(parts);
        CHECK(std::any_of(parts.begin(), parts.end(), [](const indexed_triangle_set &p) { return std::abs(its_volume(p) - 1500.) < 1.; }));
    }
}

TEST_CASE("A tilted loop around a cylinder detaches its end", "[PathCut]")
{
    // The finer cylinders have long thin side faces the cap rim runs across.
    const double r = 5.;
    for (const double segments : {36., 180., 360., 1000.}) {
        CAPTURE(segments);
        const indexed_triangle_set cyl = its_make_cylinder(r, 30., 2. * PI / segments);
        PathCut::LoopCutter lc = PathCut::make_loop_cutter(ring(r, 12, [](double a) { return 15. + 2. * std::cos(a); }), 0.2, 1.5, no_snap, false);
        REQUIRE(lc.valid());
        REQUIRE(PathCut::fit_loop_margin(cyl, lc));
        const PathCut::CutResult res = PathCut::cut_with_slab(cyl, lc);
        REQUIRE(res.error.empty());
        require_closed_parts(its_split(res.mesh));
        bool                   joined = false;
        const PathCut::Grouped g      = PathCut::group_cut(res, lc, &joined);
        CHECK_FALSE(joined);
        REQUIRE(g.n_a == 1);
        REQUIRE(g.n_b == 1);
        const double half = PI * r * r * 15.;
        CHECK_THAT(its_volume(g.side_a), WithinAbs(half, 0.1 * half));
        CHECK_THAT(its_volume(g.side_b), WithinAbs(half, 0.1 * half));
    }
}

TEST_CASE("A tilted flat loop around a long box detaches its end", "[PathCut]")
{
    // All the box's vertices are outside the cap: only its long side faces cross it.
    const indexed_triangle_set box = its_make_cube(10., 10., 100.);
    std::vector<PathCut::V3>   loop;
    const PathCut::V3          corners[4] = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}};
    for (int s = 0; s < 4; ++s)
        for (int k = 0; k < 4; ++k) {
            PathCut::V3 p = corners[s] + (corners[(s + 1) % 4] - corners[s]) * (k / 4.);
            p.z()         = 50. + 3. * (p.x() - 5.) / 5.;
            loop.push_back(p);
        }
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1.5, no_snap, true);
    REQUIRE(lc.valid());
    REQUIRE(PathCut::fit_loop_margin(box, lc));
    const PathCut::CutResult res = PathCut::cut_with_slab(box, lc);
    REQUIRE(res.error.empty());
    const auto parts = its_split(res.mesh);
    CHECK(parts.size() == 2);
    require_closed_parts(parts);
}

TEST_CASE("A flat loop cap cuts along the mean plane of the loop", "[PathCut]")
{
    const double               kerf = 0.2;
    const indexed_triangle_set cyl  = its_make_cylinder(5., 30.);
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(ring(5., 12, [](double a) { return 15. + 1.5 * std::cos(2. * a); }), kerf, 1.5, no_snap, true);
    REQUIRE(lc.valid());
    REQUIRE(PathCut::fit_loop_margin(cyl, lc));
    const PathCut::CutResult res = PathCut::cut_with_slab(cyl, lc);
    REQUIRE(res.error.empty());
    require_closed_parts(its_split(res.mesh));
    for (size_t i = res.first_cap; i < res.mesh.indices.size(); ++i)
        for (int k = 0; k < 3; ++k) {
            const PathCut::V3 p = res.mesh.vertices[size_t(res.mesh.indices[i][k])].cast<double>();
            CHECK_THAT(std::abs(lc.frame.height(p) - lc.hmean), WithinAbs(0.5 * kerf, 1e-3));
        }
}

TEST_CASE("A loop cap stops short of another part close to the loop", "[PathCut]")
{
    // An arm (r = 5) with a body 1.5 mm beside it: the cap is fitted to the arm's section, whatever the
    // requested margin, and the body is left whole.
    indexed_triangle_set model = its_make_cylinder(5., 30.);
    indexed_triangle_set body  = its_make_cube(20., 30., 30.);
    its_translate(body, Vec3f(6.5f, -15.f, 0.f));
    its_merge(model, body);
    for (const bool flat : {true, false})
        for (const double margin : {1., 5.}) {
            CAPTURE(flat, margin);
            PathCut::LoopCutter lc = PathCut::make_loop_cutter(ring(5., 8, [](double a) { return 15. + std::cos(a); }), 0.2, margin, no_snap, flat);
            REQUIRE(lc.valid());
            REQUIRE(PathCut::fit_loop_margin(model, lc));
            const PathCut::CutResult res = PathCut::cut_with_slab(model, lc);
            REQUIRE(res.error.empty());
            const auto parts = its_split(res.mesh);
            REQUIRE(parts.size() == 3);
            require_closed_parts(parts);
            CHECK(std::any_of(parts.begin(), parts.end(), [](const indexed_triangle_set &p) { return std::abs(its_volume(p) - 18000.) < 1.; }));
        }
}

TEST_CASE("A loop cap that would run into another part is refused", "[PathCut]")
{
    // Loop around one of two boxes standing 0.3 mm apart: no room for the rim between them.
    indexed_triangle_set legs = its_make_cube(4., 4., 20.);
    indexed_triangle_set leg2 = its_make_cube(4., 4., 20.);
    its_translate(leg2, Vec3f(4.3f, 0.f, 0.f));
    its_merge(legs, leg2);
    std::vector<PathCut::V3> loop = ring(2.9, 8, [](double) { return 10.; });
    for (auto &p : loop) p += PathCut::V3(2., 2., 0.);
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1., no_snap, true);
    REQUIRE(lc.valid());
    CHECK_FALSE(PathCut::fit_loop_margin(legs, lc));
    CHECK_FALSE(lc.error.empty());
    // Without fitting, a rim through the other box is not cut halfway either.
    PathCut::LoopCutter wide = PathCut::make_loop_cutter(loop, 0.2, 3., no_snap, true);
    CHECK_FALSE(PathCut::cut_with_slab(legs, wide).error.empty());
}

TEST_CASE("A loop on one side of a part, around its edge, is refused", "[PathCut]")
{
    // Three points around the edge of a long box do not go around it: the cap would slice a long sliver off.
    const indexed_triangle_set box = its_make_cube(10., 10., 100.);
    for (const bool flat : {true, false}) {
        CAPTURE(flat);
        PathCut::LoopCutter lc = PathCut::make_loop_cutter({{7., 0., 52.}, {10., 3., 50.}, {8.5, 0., 47.}}, 0.2, 1., no_snap, flat);
        REQUIRE(lc.valid());
        CHECK_FALSE(PathCut::fit_loop_margin(box, lc));
        CHECK_FALSE(lc.error.empty());
    }
}

TEST_CASE("A loop is kept when its points are uneven or snap far away", "[PathCut]")
{
    const indexed_triangle_set cyl = its_make_cylinder(5., 30.);
    // Two points almost on top of each other make the spline curl a little.
    std::vector<PathCut::V3> loop;
    for (const double a : {0.32, 0.34, 0.75, 0.83, 2.78, 2.88, 4.35, 4.91})
        loop.emplace_back(5. * std::cos(a), 5. * std::sin(a), 15. + 3.5 * std::cos(a));
    // A snap that sometimes lands on another part of the model, far from the loop.
    int  calls = 0;
    auto jumpy = [&calls](const PathCut::V3 &p) { return ++calls % 7 == 0 ? PathCut::V3(p + PathCut::V3(0., 20., 0.)) : p; };
    for (const bool flat : {true, false}) {
        CAPTURE(flat);
        PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1., jumpy, flat);
        REQUIRE(lc.valid());
        REQUIRE(PathCut::fit_loop_margin(cyl, lc));
        const PathCut::CutResult res = PathCut::cut_with_slab(cyl, lc);
        REQUIRE(res.error.empty());
        bool joined = false;
        const PathCut::Grouped g = PathCut::group_cut(res, lc, &joined);
        CHECK_FALSE(joined);
        CHECK(g.n_a == 1);
        CHECK(g.n_b == 1);
    }
}

TEST_CASE("A loop that folds back on itself is rejected", "[PathCut]")
{
    std::vector<PathCut::V3> loop; // figure eight
    for (int i = 0; i < 16; ++i) {
        const double t = 2. * PI * i / 16.;
        loop.emplace_back(10. * std::sin(t), 5. * std::sin(2. * t), 0.);
    }
    const PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1.5, no_snap, true);
    CHECK_FALSE(lc.valid());
    CHECK_FALSE(lc.error.empty());
}

TEST_CASE("A loop around a ring leaves it connected, and the connection runs round the ring", "[PathCut]")
{
    // Torus of radius 20 with a tube of radius 4 around the Z axis; a loop around the tube at +X.
    const indexed_triangle_set torus = its_make_torus(20., 4., PI / 30., PI / 12.);
    std::vector<PathCut::V3>   loop;
    for (int i = 0; i < 12; ++i) {
        const double a = 2. * PI * i / 12.;
        loop.emplace_back(20. + 4. * std::cos(a), 0., 4. * std::sin(a));
    }
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1., no_snap, true);
    REQUIRE(lc.valid());
    REQUIRE(PathCut::fit_loop_margin(torus, lc));
    const PathCut::CutResult res = PathCut::cut_with_slab(torus, lc);
    REQUIRE(res.error.empty());
    bool joined = false;
    PathCut::group_cut(res, lc, &joined);
    CHECK(joined);
    const auto lines = PathCut::connections(res);
    REQUIRE(lines.size() == 1);
    // It starts and ends next to the cut and goes past the far side of the ring.
    CHECK(std::abs(lines.front().front().y()) < 1.);
    CHECK(std::abs(lines.front().back().y()) < 1.);
    CHECK(std::any_of(lines.front().begin(), lines.front().end(), [](const PathCut::V3 &p) { return p.x() < -15.; }));
}

namespace {

// Cut with a loop snapped to the model as the tool does; returns the volumes of the two sides, or the error.
std::string loop_cut(const indexed_triangle_set &model, const std::vector<PathCut::V3> &pts, bool flat, double &vol_a, double &vol_b)
{
    const AABBMesh      aabb(model);
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(pts, 0.2, 1.5, [&aabb](const PathCut::V3 &p) {
        int i; Vec3d c; aabb.squared_distance(p, i, c); return PathCut::V3(c); }, flat);
    if (lc.valid()) PathCut::fit_loop_margin(model, lc);
    if (!lc.error.empty()) return lc.error;
    const PathCut::CutResult res = PathCut::cut_with_slab(model, lc);
    if (!res.error.empty()) return res.error;
    const PathCut::Grouped g = PathCut::group_cut(res, lc);
    vol_a = its_volume(g.side_a);
    vol_b = its_volume(g.side_b);
    return {};
}

} // namespace

TEST_CASE("A loop around a thin wide section (a blade) cuts it", "[PathCut]")
{
    for (const double w : {40., 80.})
        for (const double t : {1., 2., 4.})
            for (const bool flat : {true, false}) {
                CAPTURE(w, t, flat);
                indexed_triangle_set blade = its_make_cube(w, t, 100.);
                its_translate(blade, Vec3f(float(-w / 2), float(-t / 2), 0.f));
                const std::vector<PathCut::V3> pts{{w / 2, 0, 50}, {w / 4, t / 2, 50}, {0, t / 2, 50}, {-w / 4, t / 2, 50},
                                                   {-w / 2, 0, 50}, {-w / 4, -t / 2, 50}, {0, -t / 2, 50}, {w / 4, -t / 2, 50}};
                double a = 0., b = 0.;
                REQUIRE(loop_cut(blade, pts, flat, a, b).empty());
                const double half = w * t * (100. - 0.2) / 2.;
                CHECK_THAT(a, WithinAbs(half, 0.01 * half));
                CHECK_THAT(b, WithinAbs(half, 0.01 * half));
            }
}

TEST_CASE("A loop near the end of a bar cuts where the points are", "[PathCut]")
{
    // The spline between the points runs inside the bar, nearer its top face than its sides: snapped straight to
    // the nearest surface it would move the cut towards the top.
    for (const double side : {10., 40.})
        for (const double below : {1., 3.})
            for (const bool flat : {true, false}) {
                CAPTURE(side, below, flat);
                const indexed_triangle_set bar = its_make_cube(side, side, 100.);
                const double z = 100. - below, h = side / 2.;
                double       a = 0., b = 0.;
                REQUIRE(loop_cut(bar, {{h, 0, z}, {side, h, z}, {h, side, z}, {0, h, z}}, flat, a, b).empty());
                CHECK_THAT(std::min(a, b) / (side * side) + 0.1, WithinAbs(below, 0.05));
            }
}

TEST_CASE("A loop over a pin added as a part cuts the pin with the part", "[PathCut]")
{
    // A pin overlapping the arm right where the loop goes round, sticking out 3 mm: the same part, not another one.
    indexed_triangle_set model = its_make_cylinder(5., 40.);
    indexed_triangle_set pin   = its_make_cube(4., 3., 4.);
    its_translate(pin, Vec3f(4.f, -1.5f, 18.f));
    its_merge(model, pin);
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(ring(5., 8, [](double) { return 20.; }), 0.2, 1.5, no_snap, true);
    REQUIRE(lc.valid());
    REQUIRE(PathCut::fit_loop_margin(model, lc));
    const PathCut::CutResult res = PathCut::cut_with_slab(model, lc);
    REQUIRE(res.error.empty());
    CHECK(its_split(res.mesh).size() == 4); // the arm and the pin, each in two
}
