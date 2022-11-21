pub fn squared_norm<T>(p: &[T]) -> T
    where T: std::ops::Mul<Output=T> + std::ops::Add<Output=T> + Copy
{
    assert!(p.len() == 3);
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
}

pub fn norm<T>(
    v: &[T]) -> T
    where T: num_traits::Float
{
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn normalize<T>(
    v: &mut [T]) -> T
where T: num_traits::Float + std::ops::MulAssign
{
    let l = norm(v);
    let linv = T::one() / l;
    v[0] *= linv;
    v[1] *= linv;
    v[2] *= linv;
    l
}

pub fn cross_mut<T>(
    vo: &mut [T],
    v1: &[T],
    v2: &[T])
    where T: std::ops::Mul<Output=T> + std::ops::Sub<Output=T> + Copy
{
    vo[0] = v1[1] * v2[2] - v2[1] * v1[2];
    vo[1] = v1[2] * v2[0] - v2[2] * v1[0];
    vo[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

pub fn cross<T>(
    v1: &[T],
    v2: &[T]) -> [T; 3]
    where T: std::ops::Mul<Output=T> + std::ops::Sub<Output=T> + Copy
{
    [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]]
}

pub fn dot<T>(
    a: &[T],
    b: &[T]) -> T
    where T: std::ops::Mul<Output=T> + std::ops::Add<Output=T> + Copy
{
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn sub<T>(
    a: &[T],
    b: &[T]) -> [T; 3]
    where T: std::ops::Sub<Output=T> + Copy
{
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale<T>(
    a: &[T],
    s: T) -> [T; 3]
    where T: Copy + std::ops::Mul<Output=T> {
    [s * a[0], s * a[1], s * a[2]]
}

pub fn distance<T>(
    p0: &[T],
    p1: &[T]) -> T
    where T: num_traits::Float
{
    let v0 = p1[0] - p0[0];
    let v1 = p1[1] - p0[1];
    let v2 = p1[2] - p0[2];
    (v0 * v0 + v1 * v1 + v2 * v2).sqrt()
}