/// Find the distance of two ranges, return None if they are overlapped
///
pub fn distance_to_range<Real>(a: (Real, Real), b: (Real, Real)) -> Option<Real>
where
    Real: num_traits::Float,
{
    debug_assert!(a.0 <= a.1);
    debug_assert!(b.0 <= b.1);
    if a.0 > b.1 {
        return Some(a.0 - b.1);
    }
    if b.0 > a.1 {
        return Some(b.0 - a.1);
    }
    None
}

pub fn intersection_length(r1: &[f32; 2], r2: &[f32; 2]) -> Option<f32> {
    debug_assert!(r1[0] <= r1[1]);
    debug_assert!(r2[0] <= r2[1]);
    // separated
    if r1[1] <= r2[0] || r1[0] >= r2[1] {
        return None;
    };
    let vmin = r1[0].max(r2[0]);
    let vmax = r1[1].min(r2[1]);
    Some(vmax - vmin)
}
