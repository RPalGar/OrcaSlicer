#pragma once

#include "GLGizmoBase.hpp"
#include "GLGizmosCommon.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/Point.hpp"

#include <string>
#include <vector>

namespace Slic3r {
class ModelObject;
namespace GUI {

// Cut an object along a user-drawn path instead of a plane.
//  * View path: pick points on the object from one view; the smoothed path is extruded along the
//    view direction and splits the whole object (e.g. separate a head along the hair line).
//  * Loop: pick points around a limb on its surface; a curved cap spanning the loop detaches only
//    that region (e.g. an arm or a leg), following a non-planar loop.
// The geometry and the cut itself live in libslic3r/PathCutGeometry.hpp; the cut runs as a background job.
class GLGizmoPathCut : public GLGizmoBase
{
public:
    enum class Mode : int { ViewPath = 0, Loop = 1 };

    GLGizmoPathCut(GLCanvas3D &parent, const std::string &icon_filename, unsigned int sprite_id);

    bool on_mouse(const wxMouseEvent &mouse_event) override;
    std::string get_tooltip() const override;

protected:
    bool               on_init() override;
    std::string        on_get_name() const override;
    bool               on_is_activable() const override;
    void               on_render() override;
    void               on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;
    void               on_render_input_window(float x, float y, float bottom_limit) override;

private:
    bool raycast(const Vec2d &mouse_pos, Vec3d &hit_world) const;
    void reset_path();
    void rebuild_preview();
    void perform_cut();
    std::vector<Vec3d> preview_curve() const;
    Vec3d              sweep_dir() const;
    void               add_point(const Vec3d &hit);
    void               rebuild_cutter_preview();

    Mode               m_mode{Mode::ViewPath};
    std::vector<Vec3d> m_points;                 // picked points, world coordinates
    Vec3d              m_view_dir{Vec3d::UnitY()}; // locked at the first point (view path mode)
    float              m_kerf{0.2f};             // clearance between parts [mm]
    float              m_margin{1.5f};           // minimum reach of the loop cap past the loop, widened as needed [mm]
    enum class Sweep : int { View = 0, Level = 1, Vertical = 2 };
    Sweep              m_sweep{Sweep::View};     // view path: direction the path is swept along
    bool               m_limit_depth{false};     // view path: stop the wall `m_depth` mm behind the picked points
    float              m_depth{20.f};
    bool               m_limit_ends{false};      // view path: stop the wall `m_extend` mm past the end points
    float              m_extend{3.f};
    std::string        m_message;
    bool               m_message_is_error{false};
    const ModelObject *m_object{nullptr};        // object the points were picked on

    GLModel m_curve;
    GLModel m_sphere;
    GLModel m_cutter;                            // translucent preview of the wall / loop cap

    // A left press on the object adds a point on release, unless the mouse was dragged: then it orbits.
    bool  m_press_on_object{false};
    bool  m_orbiting{false};
    Vec2d m_press_pos{Vec2d::Zero()}, m_last_pos{Vec2d::Zero()};
    Vec3d m_press_hit{Vec3d::Zero()};
    bool    m_preview_dirty{true};
    bool    m_busy{false};                     // a cut job is running
};

} // namespace GUI
} // namespace Slic3r
