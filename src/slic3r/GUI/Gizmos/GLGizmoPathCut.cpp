#include "GLGizmoPathCut.hpp"

#include "slic3r/GUI/3DScene.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GLShader.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/I18N.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp"
#include "slic3r/GUI/MeshUtils.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/Selection.hpp"

#include "libslic3r/AABBMesh.hpp"
#include "libslic3r/Geometry.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/PathCutGeometry.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <boost/log/trivial.hpp>
#include <glad/gl.h>
#include <imgui/imgui.h>
#include <wx/utils.h>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

namespace Slic3r {
namespace GUI {

static const ColorRGBA PATH_COLOR  = {0.f, 0.588f, 0.533f, 1.f}; // Orca teal
static const ColorRGBA POINT_COLOR = {1.f, 0.55f, 0.f, 1.f};

// Messages produced by libslic3r/PathCutGeometry.hpp, listed here so xgettext picks them up.
[[maybe_unused]] static const char *const PATH_CUT_MESSAGES[] = {
    L("The path needs at least two points."),
    L("The path must not end where it starts. Use the loop mode for closed paths."),
    L("The loop needs at least three points."),
    L("The loop is degenerate (it encloses no area)."),
    L("The loop folds back on itself seen from its mean plane. Place the points more evenly around the part."),
    L("The path did not split the object. Make sure it crosses the object completely."),
    L("The loop did not detach anything. Make sure it goes all the way around the part."),
};

GLGizmoPathCut::GLGizmoPathCut(GLCanvas3D &parent, const std::string &icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{}

bool GLGizmoPathCut::on_init() { return true; }

std::string GLGizmoPathCut::on_get_name() const { return _u8L("Cut along path"); }

std::string GLGizmoPathCut::get_tooltip() const { return ""; }

bool GLGizmoPathCut::on_is_activable() const
{
    const Selection &selection = m_parent.get_selection();
    return selection.is_single_full_instance() && !selection.is_wipe_tower();
}

CommonGizmosDataID GLGizmoPathCut::on_get_requirements() const
{
    if (m_c && m_c->raycaster_ptr())
        m_c->raycaster_ptr()->set_only_support_model_part_flag(false);
    return CommonGizmosDataID(int(CommonGizmosDataID::SelectionInfo) | int(CommonGizmosDataID::InstancesHider) |
                              int(CommonGizmosDataID::Raycaster) | int(CommonGizmosDataID::ObjectClipper));
}

void GLGizmoPathCut::on_set_state()
{
    if (m_state == On || m_state == Off)
        reset_path();
    m_message.clear();
}

void GLGizmoPathCut::reset_path()
{
    m_points.clear();
    m_object        = nullptr;
    m_preview_dirty = true;
}

bool GLGizmoPathCut::raycast(const Vec2d &mouse_pos, Vec3d &hit_world) const
{
    const ModelObject *mo = m_c->selection_info() ? m_c->selection_info()->model_object() : nullptr;
    const int          inst = m_parent.get_selection().get_instance_idx();
    if (mo == nullptr || inst < 0 || inst >= int(mo->instances.size()) || m_c->raycaster() == nullptr)
        return false;
    const auto        &raycasters = m_c->raycaster()->raycasters();
    const Camera      &camera     = wxGetApp().plater()->get_camera();
    const Transform3d  inst_m     = mo->instances[inst]->get_transformation().get_matrix();
    double             best       = std::numeric_limits<double>::max();
    bool               found      = false;
    for (size_t i = 0; i < mo->volumes.size() && i < raycasters.size(); ++i) {
        if (!mo->volumes[i]->is_model_part())
            continue;
        const Transform3d trafo = inst_m * mo->volumes[i]->get_matrix();
        Vec3f             hit, normal;
        if (raycasters[i]->unproject_on_mesh(mouse_pos, trafo, camera, hit, normal, m_c->object_clipper()->get_clipping_plane())) {
            const Vec3d  w = trafo * hit.cast<double>();
            const double d = (camera.get_position() - w).squaredNorm();
            if (d < best) {
                best      = d;
                hit_world = w;
                found     = true;
            }
        }
    }
    return found;
}

bool GLGizmoPathCut::on_mouse(const wxMouseEvent &mouse_event)
{
    const Vec2d mouse_pos(mouse_event.GetX(), mouse_event.GetY());
    if (mouse_event.LeftDown() && !mouse_event.CmdDown()) {
        Vec3d hit;
        if (!raycast(mouse_pos, hit))
            return false; // let the canvas rotate the view
        const ModelObject *mo = m_c->selection_info()->model_object();
        if (m_object != mo) {
            m_points.clear();
            m_object = mo;
        }
        if (m_points.empty())
            m_view_dir = wxGetApp().plater()->get_camera().get_dir_forward();
        m_points.push_back(hit);
        m_message.clear();
        m_preview_dirty = true;
        m_parent.set_as_dirty();
        return true;
    }
    if (mouse_event.RightDown() && !m_points.empty()) {
        m_points.pop_back();
        m_preview_dirty = true;
        m_parent.set_as_dirty();
        return true;
    }
    return false;
}

std::vector<Vec3d> GLGizmoPathCut::preview_curve() const
{
    std::vector<PathCut::V3> pts(m_points.begin(), m_points.end());
    const bool closed = m_mode == Mode::Loop && pts.size() >= 3;
    std::vector<PathCut::V3> c = PathCut::catmull_rom(pts, closed, 12);
    std::vector<Vec3d>       out(c.begin(), c.end());
    if (closed && !out.empty())
        out.push_back(out.front());
    return out;
}

void GLGizmoPathCut::rebuild_preview()
{
    m_curve.reset();
    const std::vector<Vec3d> curve = preview_curve();
    if (curve.size() >= 2) {
        GLModel::Geometry g;
        g.format = {GLModel::Geometry::EPrimitiveType::LineStrip, GLModel::Geometry::EVertexLayout::P3};
        for (size_t i = 0; i < curve.size(); ++i) {
            g.add_vertex(Vec3f(curve[i].cast<float>()));
            g.add_index(unsigned(i));
        }
        m_curve.init_from(std::move(g));
    }
    if (!m_sphere.is_initialized())
        m_sphere.init_from(smooth_sphere(16, 1.f));
    m_preview_dirty = false;
}

void GLGizmoPathCut::on_render()
{
    if (m_points.empty())
        return;
    if (m_preview_dirty)
        rebuild_preview();

    const Camera &camera = wxGetApp().plater()->get_camera();
    glsafe(::glDisable(GL_DEPTH_TEST));

    if (m_curve.is_initialized()) {
        if (GLShaderProgram *shader = wxGetApp().get_shader("flat")) {
            shader->start_using();
            shader->set_uniform("view_model_matrix", camera.get_view_matrix());
            shader->set_uniform("projection_matrix", camera.get_projection_matrix());
            shader->set_uniform("width", 0.4f);
            m_curve.set_color(PATH_COLOR);
            m_curve.render();
            shader->stop_using();
        }
    }

    if (GLShaderProgram *shader = wxGetApp().get_shader("gouraud_light")) {
        shader->start_using();
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        shader->set_uniform("emission_factor", 0.5f);
        const BoundingBoxf3 bb     = m_parent.get_selection().get_bounding_box();
        const double        radius = std::clamp(bb.size().norm() * 0.006, 0.25, 2.0);
        const Transform3d  &view   = camera.get_view_matrix();
        m_sphere.set_color(POINT_COLOR);
        for (const Vec3d &p : m_points) {
            const Transform3d model = Geometry::translation_transform(p) * Geometry::scale_transform(radius);
            shader->set_uniform("view_model_matrix", view * model);
            const Matrix3d view_normal = view.matrix().block(0, 0, 3, 3) * model.matrix().block(0, 0, 3, 3).inverse().transpose();
            shader->set_uniform("view_normal_matrix", view_normal);
            m_sphere.render();
        }
        shader->stop_using();
    }
    glsafe(::glEnable(GL_DEPTH_TEST));
}

void GLGizmoPathCut::on_render_input_window(float x, float y, float bottom_limit)
{
    const float scale = m_parent.get_scale();
    ImGuiWrapper::push_toolbar_style(scale);
    GizmoImguiSetNextWIndowPos(x, y, ImGuiCond_Always, 0.0f, 0.0f);
    GizmoImguiBegin(on_get_name(), ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);

    const float input_w = ImGui::CalcTextSize("000.00").x * 2.5f;

    if (m_imgui->bbl_radio_button(_u8L("Path from view").c_str(), m_mode == Mode::ViewPath) && m_mode != Mode::ViewPath) {
        m_mode = Mode::ViewPath;
        reset_path();
    }
    ImGui::SameLine();
    if (m_imgui->bbl_radio_button(_u8L("Closed loop on surface").c_str(), m_mode == Mode::Loop) && m_mode != Mode::Loop) {
        m_mode = Mode::Loop;
        reset_path();
    }

    ImGui::Separator();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26.f);
    if (m_mode == Mode::ViewPath)
        m_imgui->text(_u8L("Click points on the object along the line where it should be split. The path is extended past the "
                           "object and cut straight through it in the direction you were looking when you placed the first point."));
    else
        m_imgui->text(_u8L("Click points on the surface all the way around the part to detach (e.g. an arm). The loop closes "
                           "itself; a curved cap spanning it separates only that part."));
    m_imgui->text(_u8L("Right click removes the last point. Drag on empty space, or Ctrl + drag, to rotate the view."));
    ImGui::PopTextWrapPos();

    ImGui::Separator();
    ImGui::AlignTextToFramePadding();
    m_imgui->text(_u8L("Clearance (mm)"));
    ImGui::SameLine(input_w * 1.6f);
    ImGui::PushItemWidth(input_w);
    ImGui::InputFloat("##pathcut_kerf", &m_kerf, 0.05f, 0.1f, "%.2f");
    m_kerf = std::clamp(m_kerf, 0.02f, 2.f);
    if (m_mode == Mode::Loop) {
        ImGui::AlignTextToFramePadding();
        m_imgui->text(_u8L("Cap margin (mm)"));
        ImGui::SameLine(input_w * 1.6f);
        ImGui::InputFloat("##pathcut_margin", &m_margin, 0.1f, 0.5f, "%.1f");
        m_margin = std::clamp(m_margin, 0.2f, 10.f);
    }
    ImGui::PopItemWidth();

    m_imgui->text(_u8L("Points") + ": " + std::to_string(m_points.size()));

    const size_t min_pts = m_mode == Mode::Loop ? 3 : 2;
    if (m_imgui->button(_L("Undo point"), ImVec2(0, 0), !m_points.empty())) {
        m_points.pop_back();
        m_preview_dirty = true;
    }
    ImGui::SameLine();
    if (m_imgui->button(_L("Clear"), ImVec2(0, 0), !m_points.empty()))
        reset_path();
    ImGui::SameLine();
    if (m_imgui->button(_L("Perform cut"), ImVec2(0, 0), m_points.size() >= min_pts))
        wxGetApp().CallAfter([this]() { perform_cut(); });

    if (!m_message.empty()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26.f);
        if (m_message_is_error)
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.f, 0.35f, 0.35f, 1.f));
        m_imgui->text(m_message);
        if (m_message_is_error)
            ImGui::PopStyleColor();
        ImGui::PopTextWrapPos();
    }

    GizmoImguiEnd();
    ImGuiWrapper::pop_toolbar_style();
}

// MCUT first (fast); CGAL as a fallback when MCUT refuses the input or fails to split it.
static indexed_triangle_set mesh_minus(const indexed_triangle_set &src, const indexed_triangle_set &cutter)
{
    std::vector<TriangleMesh> res;
    try {
        MeshBoolean::mcut::make_boolean(TriangleMesh(src), TriangleMesh(cutter), res, "A_NOT_B");
    } catch (const std::exception &e) {
        BOOST_LOG_TRIVIAL(warning) << "PathCut: MCUT failed: " << e.what();
        res.clear();
    }
    if (!res.empty() && its_split(res.front().its).size() > its_split(src).size())
        return std::move(res.front().its);
    try {
        TriangleMesh a(src);
        MeshBoolean::cgal::minus(a, TriangleMesh(cutter));
        return std::move(a.its);
    } catch (const std::exception &e) {
        BOOST_LOG_TRIVIAL(warning) << "PathCut: CGAL failed: " << e.what();
    }
    return res.empty() ? indexed_triangle_set() : std::move(res.front().its);
}

void GLGizmoPathCut::perform_cut()
{
    const Selection &selection = m_parent.get_selection();
    const int        obj_idx   = selection.get_object_idx();
    const int        inst_idx  = selection.get_instance_idx();
    Plater          *plater    = wxGetApp().plater();
    if (obj_idx < 0 || inst_idx < 0 || obj_idx >= int(plater->model().objects.size()))
        return;
    ModelObject *mo = plater->model().objects[obj_idx];
    if (mo != m_object) {
        reset_path();
        return;
    }
    const ModelInstance *mi = mo->instances[inst_idx];

    // Work in millimetres as printed: bake instance rotation/scale into the mesh, keep the offset.
    const Transform3d inst_no_offset = mi->get_transformation().get_matrix_no_offset();
    const Vec3d       offset         = mi->get_offset();
    const ModelVolume *first_part    = nullptr;
    indexed_triangle_set src;
    for (const ModelVolume *v : mo->volumes) {
        if (!v->is_model_part())
            continue;
        if (!first_part)
            first_part = v;
        indexed_triangle_set its = v->mesh().its;
        its_transform(its, inst_no_offset * v->get_matrix(), true);
        its_merge(src, its);
    }
    if (src.indices.empty())
        return;

    wxBusyCursor wait;
    PathCut::make_edge_manifold(src);
    PathCut::V3 mn = PathCut::V3::Constant(std::numeric_limits<double>::max()), mx = -mn;
    for (const auto &v : src.vertices) {
        mn = mn.cwiseMin(PathCut::V3(v.cast<double>()));
        mx = mx.cwiseMax(PathCut::V3(v.cast<double>()));
    }
    std::vector<PathCut::V3> pts;
    for (const Vec3d &p : m_points)
        pts.emplace_back(p - offset);

    PathCut::Grouped grouped;
    std::string      error;
    if (m_mode == Mode::ViewPath) {
        PathCut::ViewCutter vc = PathCut::make_view_cutter(pts, m_view_dir, mn, mx, m_kerf);
        if (!vc.valid())
            error = vc.error;
        else {
            indexed_triangle_set res = mesh_minus(src, vc.mesh);
            grouped = PathCut::group_components(its_split(res), [&vc](const PathCut::V3 &p) { return vc.side(p); }, false);
        }
    } else {
        AABBMesh aabb(src);
        auto snap = [&aabb](const PathCut::V3 &p) -> PathCut::V3 {
            int   idx;
            Vec3d c;
            aabb.squared_distance(Vec3d(p), idx, c);
            return c;
        };
        PathCut::LoopCutter lc = PathCut::make_loop_cutter(pts, m_kerf, m_margin, snap);
        if (!lc.valid())
            error = lc.error;
        else {
            indexed_triangle_set res = mesh_minus(src, lc.mesh);
            grouped = PathCut::group_components(its_split(res), [&lc](const PathCut::V3 &p) { return lc.side(p); }, true);
        }
    }
    if (error.empty() && (grouped.side_a.indices.empty() || grouped.side_b.indices.empty()))
        error = m_mode == Mode::ViewPath ?
                    "The path did not split the object. Make sure it crosses the object completely." :
                    "The loop did not detach anything. Make sure it goes all the way around the part.";
    if (!error.empty()) {
        m_message          = _u8L(error.c_str());
        m_message_is_error = true;
        m_parent.set_as_dirty();
        return;
    }

    // Detached piece (smaller) first.
    std::vector<indexed_triangle_set *> sides{&grouped.side_a, &grouped.side_b};
    if (std::abs(its_volume(grouped.side_a)) > std::abs(its_volume(grouped.side_b)))
        std::swap(sides[0], sides[1]);

    Model                    tmp;
    const std::vector<std::string> suffix{_u8L("piece"), _u8L("rest")};
    for (size_t i = 0; i < sides.size(); ++i) {
        ModelObject *no = tmp.add_object(*mo);
        no->clear_volumes();
        no->clear_instances();
        no->cut_connectors.clear();
        no->cut_id.invalidate();
        no->name = mo->name + " - " + suffix[i];
        ModelVolume *nv = no->add_volume(TriangleMesh(std::move(*sides[i])));
        nv->name = no->name;
        if (first_part) {
            nv->config.apply(first_part->config);
            nv->set_material_id(first_part->material_id());
        }
        ModelInstance *ni = no->add_instance();
        ni->set_offset(offset);
        no->invalidate_bounding_box();
        no->ensure_on_bed();
    }
    // Both pieces now rest on the bed; move the detached one beside the rest so they don't overlap.
    if (tmp.objects.size() == 2) {
        const BoundingBoxf3 piece = tmp.objects[0]->instance_bounding_box(0);
        const BoundingBoxf3 rest  = tmp.objects[1]->instance_bounding_box(0);
        ModelInstance      *pi    = tmp.objects[0]->instances[0];
        pi->set_offset(pi->get_offset() + Vec3d(rest.max.x() - piece.min.x() + 5., 0., 0.));
        tmp.objects[0]->invalidate_bounding_box();
    }

    m_parent.reset_all_gizmos();
    {
        Plater::TakeSnapshot snapshot(plater, _u8L("Cut along path"));
        plater->apply_cut_object_to_model(size_t(obj_idx), tmp.objects);
    }
    if (grouped.n_small > 0)
        BOOST_LOG_TRIVIAL(info) << "PathCut: " << grouped.n_small << " tiny fragments kept with their side";
}

} // namespace GUI
} // namespace Slic3r
