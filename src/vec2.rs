//! methods for 2D vector

use num_traits::AsPrimitive;

pub fn sub_<T>(a: &[T; 2], b: &[T; 2]) -> [T; 2]
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

pub fn to_array_from_vtx2xy<T, Index>(vtx2xyz: &[T], i_vtx: Index) -> [T; 2]
    where
        T: Copy,
        Index: AsPrimitive<usize>
{
    let i_vtx: usize = i_vtx.as_();
    [
        vtx2xyz[i_vtx * 2],
        vtx2xyz[i_vtx * 2 + 1],
    ]
}

// -------------------------------
pub fn rotate90<T>(v: &nalgebra::Vector2<T>) -> nalgebra::Vector2<T>
where
    T: nalgebra::RealField + Copy,
{
    nalgebra::Vector2::<T>::new(-v[1], v[0])
}

pub fn to_na<T>(vtx2xyz: &[T], i_vtx: usize) -> nalgebra::Vector2<T>
where
    T: Copy + nalgebra::RealField,
{
    nalgebra::Vector2::<T>::from_row_slice(&vtx2xyz[i_vtx * 2..(i_vtx + 1) * 2])
}

pub fn norm_squared<T>(v: &nalgebra::Vector2<T>) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    v[0] * v[0] + v[1] * v[1]
}

pub fn basis<T>(i_dim: usize, eps: T) -> nalgebra::Vector2<T>
where
    T: nalgebra::RealField,
{
    let mut b = nalgebra::Vector2::<T>::zeros();
    b[i_dim] = eps;
    b
}
