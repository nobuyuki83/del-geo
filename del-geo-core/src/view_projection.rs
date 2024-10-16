

pub struct Perspective {
    pub lens: f32,
    pub near: f32,
    pub far: f32,
    pub proj_direction: bool,
    pub cam_pos: [f32;3],
    pub scale: f32,
}

impl Perspective {

    /// Proj * Translation * Scale
    pub fn mat4_col_major(&self, aspect_ratio: f32) -> [f32;16] {
        let cam_projection = crate::mat4_col_major::camera_perspective_blender(
            aspect_ratio,
            self.lens, self.near, self.far, self.proj_direction);
        let transl =
            crate::mat4_col_major::from_translate(&[-self.cam_pos[0], -self.cam_pos[1], -self.cam_pos[2]]);
        let scale = crate::mat4_col_major::from_scale_uniform(self.scale);
        let ts = crate::mat4_col_major::mult_mat(&transl, &scale);
        crate::mat4_col_major::mult_mat(&cam_projection, &ts)
    }

    pub fn camera_translation(
        &mut self,
        asp: f32,
        cursor_dx: f32,
        cursor_dy: f32,
    ) {
        let mp = self.mat4_col_major(asp);
        let sx = (mp[3 + 4 * 3] - mp[0 + 4 * 3]) / mp[0 + 4 * 0];
        let sy = (mp[3 + 4 * 3] - mp[1 + 4 * 3]) / mp[1 + 4 * 1];
        self.cam_pos[0] -= sx * cursor_dx as f32;
        self.cam_pos[1] -= sy * cursor_dy as f32;
    }
}