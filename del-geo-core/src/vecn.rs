pub trait Arr<T, const N: usize> {
    fn add(self, other: &[T; N]) -> Self;
    fn sub(self, other: &[T; N]) -> Self;
    fn norm(self) -> T;
}

impl<T, const N: usize> Arr<T, N> for [T; N]
where
    T: num_traits::Float + Copy + std::iter::Sum,
{
    fn add(self, other: &[T; N]) -> Self {
        std::array::from_fn(|i| self[i] + other[i])
    }
    fn sub(self, other: &[T; N]) -> Self {
        std::array::from_fn(|i| self[i] - other[i])
    }
    fn norm(self) -> T {
        self.iter().map(|&v| v * v).sum::<T>().sqrt()
    }
}

#[test]
fn test_add() {
    // Test with different array sizes
    assert_eq!([1., 2.].add(&[2., 3.]), [3., 5.]);
    assert_eq!([1., 2.].sub(&[2., 3.]), [-1., -1.]);
    assert_eq!([4., 3.].norm(), 5.);
    // assert_eq!([1, 2, 3].add(&[2, 3, 4]), [3, 5, 7]);
    // assert_eq!([1, 2, 3, 4].add(&[2, 3, 4, 5]), [3, 5, 7, 9]);
}
