//! methods for 2x2 matrix
pub fn polar_decomposition<T>(
    m: &nalgebra::Matrix2<T>,
) -> (nalgebra::Matrix2<T>, nalgebra::Matrix2<T>)
where
    T: nalgebra::RealField + Copy,
{
    let x = m[(0, 0)] + m[(1, 1)];
    let y = m[(1, 0)] - m[(0, 1)];
    let scale = T::one() / (x * x + y * y).sqrt();
    let c = x * scale;
    let s = y * scale;
    let u_mat = nalgebra::Matrix2::<T>::new(c, -s, s, c);
    let p_mat = u_mat * m;
    (u_mat, p_mat)
}
