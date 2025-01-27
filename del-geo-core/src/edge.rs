//! methods for 2D or 3D edge (line segment)

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
