#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <libslic3r/PathCutGeometry.hpp>
#include <libslic3r/TriangleMesh.hpp>

#include <cmath>
#include <map>
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

    const PathCut::CutResult res = PathCut::cut_with_slab(cube, vc, false);
    REQUIRE(res.error.empty());
    const auto parts = its_split(res.mesh);
    REQUIRE(parts.size() == 2);
    require_closed_parts(parts);

    const PathCut::Grouped g = PathCut::group_components(parts, [&vc](const PathCut::V3 &p) { return vc.side(p); }, false);
    CHECK(g.n_a == 1);
    CHECK(g.n_b == 1);
    // The two halves add up to the cube minus a wall kerf thick, 20 mm deep, as long as the path inside the cube.
    const double total = its_volume(g.side_a) + its_volume(g.side_b);
    CHECK(total < 8000.);
    CHECK(total > 8000. - kerf * 20. * 40.);
}

TEST_CASE("A tilted loop around a cylinder detaches its end", "[PathCut]")
{
    // Cylinder r = 5, h = 30, loop around it tilted between z = 13 and z = 17.
    const double               r = 5., kerf = 0.2;
    const indexed_triangle_set cyl = its_make_cylinder(r, 30.);
    std::vector<PathCut::V3>   loop;
    for (int i = 0; i < 12; ++i) {
        const double a = 2. * PI * i / 12.;
        loop.emplace_back(r * std::cos(a), r * std::sin(a), 15. + 2. * std::cos(a));
    }
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, kerf, 1.5, [](const PathCut::V3 &p) { return p; });
    REQUIRE(lc.valid());
    REQUIRE(PathCut::fit_loop_margin(cyl, lc));

    const PathCut::CutResult res = PathCut::cut_with_slab(cyl, lc, true);
    REQUIRE(res.error.empty());
    const auto parts = its_split(res.mesh);
    REQUIRE(parts.size() == 2);
    require_closed_parts(parts);

    const PathCut::Grouped g = PathCut::group_components(parts, [&lc](const PathCut::V3 &p) { return lc.side(p); }, true);
    REQUIRE(g.n_a == 1);
    REQUIRE(g.n_b == 1);
    // Each side holds about half the cylinder.
    const double half = PI * r * r * 15.;
    CHECK_THAT(its_volume(g.side_a), WithinAbs(half, 0.1 * half));
    CHECK_THAT(its_volume(g.side_b), WithinAbs(half, 0.1 * half));
}

TEST_CASE("A loop that misses part of the cross-section is reported instead of cutting halfway", "[PathCut]")
{
    // Loop around one leg of a U made of two boxes joined at the bottom, but with a cap so wide it
    // reaches the other leg: the cap rim runs through it.
    indexed_triangle_set u = its_make_cube(4., 4., 20.);
    indexed_triangle_set leg2 = its_make_cube(4., 4., 20.);
    its_translate(leg2, Vec3f(6.f, 0.f, 0.f));
    its_merge(u, leg2);
    std::vector<PathCut::V3> loop;
    for (int i = 0; i < 8; ++i) {
        const double a = 2. * PI * i / 8.;
        loop.emplace_back(2. + 2.9 * std::cos(a), 2. + 2.9 * std::sin(a), 10.);
    }
    PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 5., [](const PathCut::V3 &p) { return p; });
    REQUIRE(lc.valid());
    const PathCut::CutResult res = PathCut::cut_with_slab(u, lc, true);
    CHECK_FALSE(res.error.empty());
}

TEST_CASE("A loop that folds back on itself is rejected", "[PathCut]")
{
    // Figure-eight seen from its mean plane.
    std::vector<PathCut::V3> loop;
    for (int i = 0; i < 16; ++i) {
        const double t = 2. * PI * i / 16.;
        loop.emplace_back(10. * std::sin(t), 5. * std::sin(2. * t), 0.);
    }
    const PathCut::LoopCutter lc = PathCut::make_loop_cutter(loop, 0.2, 1.5, [](const PathCut::V3 &p) { return p; });
    CHECK_FALSE(lc.valid());
    CHECK_FALSE(lc.error.empty());
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
