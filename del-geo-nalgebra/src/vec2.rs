pub fn rotate90<T>(v: &nalgebra::Vector2<T>) -> nalgebra::Vector2<T>
where
    T: nalgebra::RealField + Copy,
{
    nalgebra::Vector2::<T>::new(-v[1], v[0])
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
