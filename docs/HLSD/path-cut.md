# Cut along path — High Level Design

## Purpose and scope

Cut along path splits an object along a line the user draws on it, where the planar Cut tool can only
split along a plane. It is meant for figures and organic models (sculpted, scanned or AI-generated
meshes, often with millions of triangles) that have to be printed in pieces: a head along the hair line,
a figure along the folds of its clothes.

The line is picked in one of two shapes:

- **Open path**: points picked from one view direction. The smoothed path is swept along the view direction
  into a thin wall that splits the whole object. Optional limits (start, end, depth) stop the wall inside the
  object; then it only cuts the parts it goes all the way through.
- **Closed loop**: points picked on the surface around the region to detach (an arm at the shoulder, a head
  along the hair line). A cap spanning the loop, flat or following the loop's heights, detaches only that
  region, so the parts next to it stay whole.

## Geometry

`src/libslic3r/PathCut.hpp` / `.cpp`, tested by `tests/libslic3r/test_path_cut.cpp`.

A general mesh boolean is far too slow on the meshes this tool is for, so the cut does not use one. The
cut surface is a `Cutter`: a signed field `f` over the mesh, the distance across the wall (`ViewCutter`) or
the height above the cap (`LoopCutter`), and `|f| < kerf / 2` is the clearance left between the pieces.
`cut_with_slab` refines the triangles near the slab where the field bends or where they cross the rim of a
limited slab, clips them against the two level sets `f = ±kerf / 2`, with each edge crossing computed once
so the result stays watertight, and closes the openings with caps triangulated in a 2D chart of each level
set (arc length along the path and height along the view direction for the wall, the cap seen along its
normal for the loop). Cost is linear in the triangle count and only the triangles near the slab are
touched.

Rules that shape the result:

- A limited wall that ends inside a part does not cut it: only the parts the wall goes all the way through
  are split, and when the wall ends inside every part it reaches, the cut is refused with a message.
- The loop's cap is a height field over evenly spaced spokes from its centre (the outermost crossing of each
  spoke with the loop, so a loop snapped to a bumpy surface may wobble), smoothed so it needs few triangles.
  A loop that crosses itself seen from its mean plane is refused; small curls from points placed close
  together are tolerated. The smoothing never pulls the loop inwards, so a long thin section (a blade) keeps
  its tips. Snaps that move a point of the loop far (onto a nearby part, or across a thin one)
  are not kept, and only their part along the mean plane is: snapped up or down onto a nearby face, the
  loop would move the cut.
- The loop is the cut: only what it encloses is cut, and the rest of the model is left whole. The cap only
  reaches its rim, which `fit_loop_margin` places from where the surface crosses the cap: the curve around
  the centre that the loop follows is the section of the part the loop is on, and so are curves overlapping
  it (a pin added as a part, overlapping shells of a sculpted mesh); the others belong to other parts, also
  when they go around the centre (a fist around a grip). The crossings with the two faces of the cut count
  too, for surfaces lying almost along the cap (the underside of the fist just above the loop). Per spoke the
  rim is put the cap margin past the section, or halfway to another part when one is nearer. A section
  reaching farther past the loop than the loop's mean radius (measured to the loop itself) means the loop
  does not go around the part: the tool then shows that section in red. No room between the section and
  another part means the cap would run into it: both are refused.
- `make_edge_manifold` runs first: zero-thickness sheets and edges shared by more than two faces, common
  in sculpted and AI-generated meshes, would leave openings the caps cannot close.
- `group_cut` assigns every resulting part to a side by the cut face it carries; a part the cut does not touch
  goes to the side it is on, or, next to a loop, with the piece it lies in. A part carrying both cut
  faces is still connected around the cut (an arm whose hand rests on the hip): `connections` finds the
  shortest line over its surface from one cut face to the other, the tool shows it, and the cut can be
  performed so that a second cut across that line detaches the part.
- Each part of a multi-part object is cut on its own and its pieces keep the part's settings and filament.
  A part the cut does not reach goes to the side it is on or, away from a loop cap, with the piece it lies
  in, else with the nearest one (a hat on the detached head).

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
