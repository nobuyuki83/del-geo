pub trait Arr<T, const N: usize> {
    fn add(self, other: T) -> Self;
}

impl<T, const N: usize> Arr<[T; N], N> for [T; N]
where
    T: std::ops::Add<Output = T> + Copy,
{
    fn add(self, other: [T; N]) -> Self {
        std::array::from_fn(|i| self[i] + other[i])
    }
}

#[test]
fn test_add() {
    // Test with different array sizes
    assert_eq!([1, 2].add([2, 3]), [3, 5]);
    assert_eq!([1, 2, 3].add([2, 3, 4]), [3, 5, 7]);
    assert_eq!([1, 2, 3, 4].add([2, 3, 4, 5]), [3, 5, 7, 9]);
}
