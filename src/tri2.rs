use num_traits::AsPrimitive;

pub fn area_<T>(p0: &[T], p1: &[T], p2: &[T]) -> T
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + 'static + Copy,
          f64: AsPrimitive<T>
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    0.5_f64.as_() * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
}

pub fn dldx_(
    p0: &[f32],
    p1: &[f32],
    p2: &[f32]) -> [[f32; 3]; 3] {
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let a0 = area_(p0, p1, p2);
    let tmp1 = 0.5 / a0;
    [
        [
            tmp1 * (p1[1] - p2[1]),
            tmp1 * (p2[1] - p0[1]),
            tmp1 * (p0[1] - p1[1]),
        ],
        [
            tmp1 * (p2[0] - p1[0]),
            tmp1 * (p0[0] - p2[0]),
            tmp1 * (p1[0] - p0[0]),
        ],
        [
            tmp1 * (p1[0] * p2[1] - p2[0] * p1[1]),
            tmp1 * (p2[0] * p0[1] - p0[0] * p2[1]),
            tmp1 * (p0[0] * p1[1] - p1[0] * p0[1]),
        ]
    ]
}