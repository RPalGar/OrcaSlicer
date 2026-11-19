# Cut along path — High Level Design

## Purpose and scope

Cut along path splits an object along a line the user draws on it, where the planar Cut tool can only
split along a plane. It is meant for figures and organic models (sculpted, scanned or AI-generated
meshes, often with millions of triangles) that have to be printed in pieces: a head along the hair line,
a figure along the folds of its clothes.

The points are picked from one view direction. The smoothed path is swept along the view direction into
a thin wall that splits the whole object. Optional limits (start, end, depth) stop the wall inside the
object; then it only cuts the parts it goes all the way through.

## Geometry

`src/libslic3r/PathCut.hpp` / `.cpp`, tested by `tests/libslic3r/test_path_cut.cpp`.

A general mesh boolean is far too slow on the meshes this tool is for, so the cut does not use one. The
wall is a signed field `f` over the mesh, the distance across it (`ViewCutter`), and `|f| < kerf / 2` is
the clearance left between the pieces. `cut_with_slab` refines the triangles near the wall where the field
bends or where they cross the rim of a limited wall, clips them against the two level sets
`f = ±kerf / 2`, with each edge crossing computed once so the result stays watertight, and closes the
openings with caps triangulated in a 2D chart of each level set (arc length along the path, height along
the view direction). Cost is linear in the triangle count and only the triangles near the wall are
touched.

Rules that shape the result:

- A limited wall that ends inside a part does not cut it: only the parts the wall goes all the way through
  are split, and when the wall ends inside every part it reaches, the cut is refused with a message.
- `make_edge_manifold` runs first: zero-thickness sheets and edges shared by more than two faces, common
  in sculpted and AI-generated meshes, would leave openings the caps cannot close.
- `group_cut` assigns every resulting part to a side by the cut face it carries. A part carrying both cut
  faces is still connected around the wall; the cut is then reported as not separating the object.

The caps are triangulated with earcut followed by Delaunay edge flips rather than with
`Triangulation`: where the cut grazes the surface the rings of a cap can touch or cross themselves at
micrometre scale, which `Triangulation` does not accept (it needs unique points and non-crossing
constraints). The pieces earcut leaves open there are closed with fans.

## GUI

Cut along path is the Curve mode of the Cut tool, next to Planar and Dovetail in its row of mode buttons; it
has no toolbar icon of its own. In Curve mode `GLGizmoCut3D` hands rendering, the mouse, the keyboard
(Enter: preview, Esc: back to editing), its window and its common data requirements to
`GLGizmoCut3D::Curve`, which works on points picked on the mesh rather than on the cut plane the rest of the
tool is built around.

`Curve` collects the points, runs the cut as a job on the Plater's UI job worker on a copy of the mesh, shows
a preview with the pieces apart, and replaces the object with the pieces only when the user performs the
cut, as one undo step. Volumes that are not cut (negative volumes, modifiers) go with the pieces they are in.
The job is dropped when the tool is closed or cancelled, and the preview when the object changes meanwhile.
The tool keeps no state in the project file.
