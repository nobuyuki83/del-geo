//! methods for 2D vector

pub fn sub_<T>(
    a: &[T;2],
    b: &[T;2]) -> [T; 2]
    where T: std::ops::Sub<Output=T> + Copy
{
    [a[0] - b[0], a[1] - b[1]]
}


pub fn to_na<T>(vtx2xyz: &[T], i_vtx: usize) -> nalgebra::Vector2::<T>
    where T: Copy + nalgebra::RealField
{
    nalgebra::Vector2::<T>::from_row_slice(&vtx2xyz[i_vtx *2..(i_vtx +1)*2])
}

pub fn norm_squared<T>(v: &nalgebra::Vector2::<T>) -> T
where T: std::ops::Mul<Output=T> + std::ops::Add<Output = T> + Copy
{
    v[0]*v[0] + v[1]*v[1]
}