/// Find the distance of two ranges, return None if they are overlapped
pub fn distance_to_range<Real>(a: (Real, Real), b: (Real, Real)) -> Option<Real>
where
    Real: num_traits::Float,
{
    if a.0 > b.1 {
        return Some(a.0 - b.1);
    }
    if b.0 > a.1 {
        return Some(b.0 - a.1);
    }
    None
}
