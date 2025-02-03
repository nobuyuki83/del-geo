pub trait Mat2ColMajor<T: num_traits::Float>
where
    Self: Sized,
{
    fn scale(&self, s: T) -> Self;
    fn mult_mat_col_major(&self, other: &Self) -> Self;
    fn mult_vec(&self, other: &[T; 2]) -> [T; 2];
    fn add(&self, other: &Self) -> Self;
    fn transpose(&self) -> Self;
}

impl<Real> Mat2ColMajor<Real> for [Real; 4]
where
    Real: num_traits::Float,
{
    fn scale(&self, s: Real) -> Self {
        std::array::from_fn(|i| self[i] * s)
    }

    fn add(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i] + other[i])
    }

    fn mult_mat_col_major(&self, other: &Self) -> Self {
        [
            self[0] * other[0] + self[2] * other[1],
            self[1] * other[0] + self[3] * other[1],
            self[0] * other[2] + self[2] * other[3],
            self[1] * other[2] + self[3] * other[3],
        ]
    }

    fn mult_vec(&self, other: &[Real; 2]) -> [Real; 2] {
        mult_vec(self, other)
    }

    fn transpose(&self) -> Self {
        [self[0], self[2], self[1], self[3]]
    }
}

pub fn from_outer_product<T>(a: &[T; 2], b: &[T; 2]) -> [T; 4]
where
    T: num_traits::Float,
{
    [a[0] * b[0], a[1] * b[0], a[0] * b[1], a[1] * b[1]]
}

pub fn mult_vec<T>(a: &[T; 4], b: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    [a[0] * b[0] + a[2] * b[1], a[1] * b[0] + a[3] * b[1]]
}

pub fn from_identity<T>() -> [T; 4]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, one]
}

pub fn add_four<T>(a: &[T; 4], b: &[T; 4], c: &[T; 4], d: &[T; 4]) -> [T; 4]
where
    T: num_traits::Float,
{
    std::array::from_fn(|i| a[i] + b[i] + c[i] + d[i])
}
