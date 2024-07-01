//! methods for 2D or 3D edge (line segment)

pub fn length_<T, const N: usize>(p0: &[T], p1: &[T]) -> T
where
    T: num_traits::Float + std::ops::AddAssign,
{
    let mut x = T::zero();
    for i in 0..N {
        x += (p0[i] - p1[i]) * (p0[i] - p1[i]);
    }
    x.sqrt()
}
