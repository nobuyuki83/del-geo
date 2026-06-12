//! methods for 2D or 3D edge (line segment)

pub const EDGE_FACE2IDX: [usize; 3] = [0, 1, 2];
pub const EDGE_IDX2NODE: [usize; 2] = [1, 0];

pub fn length<T, const N: usize>(p0: &[T; N], p1: &[T; N]) -> T
where
    T: num_traits::Float,
{
    let mut x = T::zero();
    for i in 0..N {
        x = x + (p0[i] - p1[i]) * (p0[i] - p1[i]);
    }
    x.sqrt()
}
