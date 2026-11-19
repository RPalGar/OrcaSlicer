#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <libslic3r/PathCut.hpp>
#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/libslic3r.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
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
