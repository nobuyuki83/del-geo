//! View rotation by trackball UI

use crate::quaternion::Quaternion;
#[derive(Debug, Clone, Copy)]
pub struct Trackball<Real>
where
    Real: num_traits::Float,
{
    pub quaternion: [Real; 4],
}

impl<Real> Trackball<Real>
where
    Real: num_traits::Float,
{
    pub fn new() -> Self {
        let zero = Real::zero();
        let one = Real::one();
        Self {
            quaternion: [zero, zero, zero, one],
        }
    }
    pub fn mat4_col_major(&self) -> [Real; 16] {
        self.quaternion.to_mat4_col_major()
    }

    pub fn camera_rotation(&mut self, cursor_dx: Real, cursor_dy: Real) {
        let dx = cursor_dx;
        let dy = cursor_dy;
        let a = (dx * dx + dy * dy).sqrt();
        let zero = Real::zero();
        if a.is_zero() {
            return;
        }
        let dq = crate::quaternion::from_axisangle(&[-dy, dx, zero]).normalized();
        self.quaternion = dq.mult_quaternion(&self.quaternion);
    }
}

impl<Real> Default for Trackball<Real>
where
    Real: num_traits::Float,
{
    fn default() -> Self {
        Self::new()
    }
}
