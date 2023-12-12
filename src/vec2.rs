//! methods for 2D vector

pub fn sub_<T>(
    a: &[T],
    b: &[T]) -> [T; 2]
    where T: std::ops::Sub<Output=T> + Copy
{
    [a[0] - b[0], a[1] - b[1]]
}