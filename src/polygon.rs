use num_traits::AsPrimitive;

pub fn arclength_polygon<T, const X: usize>(
    vtxs: &Vec<nalgebra::base::SVector<T, X>>) -> T
    where T: nalgebra::RealField + Copy,
          f64: num_traits::AsPrimitive<T>
{
    if vtxs.len() < 2 { return 0_f64.as_(); }
    let np = vtxs.len();
    let mut len: T = 0_f64.as_();
    for ip0 in 0..np {
        let ip1 = (ip0 + 1) % np;
        len += (vtxs[ip0] - vtxs[ip1]).norm();
    }
    return len;
}