//! Perspective projection looking at -Z direction

use num_traits::AsPrimitive;

use crate::mat4_col_major::Mat4ColMajor;
#[derive(Debug, Clone, Copy)]
pub struct Perspective<Real> {
    pub lens: Real,
    pub near: Real,
    pub far: Real,
    pub proj_direction: bool,
    pub cam_pos: [Real; 3],
    pub scale: Real,
}

impl<Real> Perspective<Real>
where
    Real: num_traits::Float + num_traits::FloatConst + 'static,
    f64: AsPrimitive<Real>,
{
    /// Proj * Translation * Scale
    pub fn mat4_col_major(&self, aspect_ratio: Real) -> [Real; 16] {
        use crate::vec3::Vec3;
        let cam_projection = crate::mat4_col_major::camera_perspective_blender(
            aspect_ratio,
            self.lens,
            self.near,
            self.far,
            self.proj_direction,
        );
        let transl = crate::mat4_col_major::from_translate(&self.cam_pos.scale(-Real::one()));
        let scale = crate::mat4_col_major::from_scale_uniform(self.scale);
        cam_projection.mult_mat(&transl.mult_mat(&scale))
    }

    pub fn camera_translation(&mut self, asp: Real, cursor_dx: Real, cursor_dy: Real) {
        let mp = self.mat4_col_major(asp);
        let sx = (mp[3 + 4 * 3] - mp[12]) / mp[0]; // (mp[3 + 4 * 3] - mp[0 + 4 * 3]) / mp[0 + 4 * 0];
        let sy = (mp[3 + 4 * 3] - mp[13]) / mp[5]; // (mp[3 + 4 * 3] - mp[1 + 4 * 3]) / mp[1 + 4 * 1];
        self.cam_pos[0] = self.cam_pos[0] - sx * cursor_dx;
        self.cam_pos[1] = self.cam_pos[1] - sy * cursor_dy;
    }
}
