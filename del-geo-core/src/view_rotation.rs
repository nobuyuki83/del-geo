pub struct Trackball {
    pub quaternion: [f32; 4],
}

impl Trackball {
    pub const fn new() -> Self {
        Self {
            quaternion: [0., 0., 0., 1.],
        }
    }

    pub fn mat4_col_major(&self) -> [f32; 16] {
        crate::quaternion::to_mat4_col_major(&self.quaternion)
    }

    pub fn camera_rotation(&mut self, cursor_dx: f64, cursor_dy: f64) {
        let dx = cursor_dx as f32;
        let dy = cursor_dy as f32;
        let a = (dx * dx + dy * dy).sqrt();
        if a == 0.0 {
            return;
        }
        let dq = crate::quaternion::normalized(&crate::quaternion::from_axisangle(&[-dy, dx, 0.]));
        self.quaternion = crate::quaternion::mult_quaternion(&dq, &self.quaternion);
    }
}

impl Default for Trackball {
    fn default() -> Self {
        Self::new()
    }
}
