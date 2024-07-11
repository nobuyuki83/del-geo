//! methods for 2D vector

pub fn sub<T>(a: &[T; 2], b: &[T; 2]) -> [T; 2]
where
    T: std::ops::Sub<Output = T> + Copy,
{
    [a[0] - b[0], a[1] - b[1]]
}

pub fn from_homogeneous<Real>(v: &[Real; 3]) -> Option<[Real; 2]>
where
    Real: num_traits::Float,
{
    if v[2].is_zero() {
        return None;
    }
    Some([v[0] / v[2], v[0] / v[2]])
}

// -------------------------------

pub struct XY<'a, Real> {
    pub p: &'a [Real; 2],
}

impl<'a, Real> XY<'a, Real>
where
    Real: num_traits::Float,
{
    pub fn aabb(&self) -> [Real; 4] {
        [self.p[0], self.p[1], self.p[0], self.p[1]]
    }
}
