use num_traits::AsPrimitive;

pub struct ControlPoints<'a, Real, const N: usize> {
    pub p0: &'a [Real; N],
    pub p1: &'a [Real; N],
    pub p2: &'a [Real; N],
    pub p3: &'a [Real; N],
}

pub fn eval<Real, const N: usize>(
    p0: &[Real; N],
    p1: &[Real; N],
    p2: &[Real; N],
    p3: &[Real; N],
    t0: Real,
) -> [Real; N]
where
    Real: num_traits::Float + Copy + std::iter::Sum,
{
    let one = Real::one();
    let three = one + one + one;
    let t1 = one - t0;
    use crate::vecn::VecN;
    crate::vecn::add_four_vectors(
        &p0.scale(t1 * t1 * t1),
        &p1.scale(three * t0 * t1 * t1),
        &p2.scale(three * t0 * t0 * t1),
        &p3.scale(t0 * t0 * t0),
    )
}

pub fn sample_uniform_param<Real, const N: usize>(
    ndiv: usize,
    p0: &[Real; N],
    p1: &[Real; N],
    p2: &[Real; N],
    p3: &[Real; N],
    is_include_endpoint_start: bool,
    is_include_endpoint_end: bool,
) -> Vec<[Real; N]>
where
    Real: Copy + 'static + num_traits::Float + std::iter::Sum,
    usize: AsPrimitive<Real>,
{
    let mut ret: Vec<[Real; N]> = vec![];
    if is_include_endpoint_start {
        ret.push(*p0);
    }
    for idiv in 1..ndiv {
        let t0: Real = idiv.as_() / ndiv.as_();
        let p0 = eval(p0, p1, p2, p3, t0);
        ret.push(p0);
    }
    if is_include_endpoint_end {
        ret.push(*p3);
    }
    ret
}

pub fn arclength_from_vtx2vecn<T, const N: usize>(vtxs: &[[T; N]]) -> T
where
    T: num_traits::Float + Copy + std::iter::Sum,
{
    use crate::vecn::VecN;
    if vtxs.len() < 2 {
        return T::zero();
    }
    let mut len: T = T::zero();
    for ip0 in 0..vtxs.len() - 1 {
        len = len + vtxs[ip0].sub(&vtxs[ip0 + 1]).norm();
    }
    len
}
pub fn sample_uniform_length<Real, const N: usize>(
    cps: ControlPoints<Real, N>,
    target_edge_length: Real,
    is_include_endpoint_start: bool,
    is_include_endpoint_end: bool,
    ndiv_sample: usize,
) -> Vec<[Real; N]>
where
    Real: Copy
        + 'static
        + num_traits::Float
        + AsPrimitive<usize>
        + std::iter::Sum
        + std::fmt::Debug
        + std::fmt::Display,
    usize: AsPrimitive<Real>,
    f64: AsPrimitive<Real>,
{
    use crate::vecn::VecN;
    let mut ret: Vec<[Real; N]> = vec![];
    if is_include_endpoint_start {
        ret.push(*cps.p0);
    }
    let ps = sample_uniform_param(ndiv_sample, cps.p0, cps.p1, cps.p2, cps.p3, true, true);
    assert_eq!(ps.len(), ndiv_sample + 1);
    let len = arclength_from_vtx2vecn(&ps);
    let ndiv_out: usize = (len / target_edge_length).ceil().as_();
    let elen: Real = len / ndiv_out.as_();
    // dbg!(len, ndiv_out, elen, target_edge_length);
    let mut len_to_go = elen;
    let mut traveled_len_in_edge = Real::zero();
    let mut i_div = 0;
    loop {
        if is_include_endpoint_start {
            if ret.len() == ndiv_out {
                break;
            }
        } else if ret.len() == ndiv_out - 1 {
            break;
        }
        if i_div == ps.len() - 1 {
            break;
        }
        let len_edge = ps[i_div + 1].sub(&ps[i_div]).norm();
        // println!("{} {} {} {} {}", i_div, len_to_go, traveled_len_in_edge, elen, len_edge);
        assert!(len_edge > traveled_len_in_edge);
        if len_edge - traveled_len_in_edge >= len_to_go {
            // there is a sampled point in this edge
            (traveled_len_in_edge, len_to_go) = (traveled_len_in_edge + len_to_go, elen);
            {
                // output point
                let r0 = traveled_len_in_edge / len_edge;
                assert!(r0 >= Real::zero() && r0 <= Real::one(), "{}", r0);
                let t = (i_div.as_() + r0) / ndiv_sample.as_();
                assert!(
                    t < Real::one(),
                    "t={} idiv={} ndiv_sample={} r0={}",
                    t,
                    i_div,
                    ndiv_sample,
                    r0
                );
                assert!(t > Real::zero(), "t={}", t);
                let q = eval(cps.p0, cps.p1, cps.p2, cps.p3, t);
                ret.push(q);
            }
        } else {
            // move next edge
            len_to_go = len_to_go - (len_edge - traveled_len_in_edge);
            traveled_len_in_edge = Real::zero();
            i_div += 1;
        }
    }
    if is_include_endpoint_end {
        ret.push(*cps.p3);
    }
    ret
}

#[test]
fn test() {
    use crate::vec2::Vec2;
    let p0 = [0.1f32, 0.2];
    let p1 = [0.4, 0.3];
    let p2 = [1.1, 1.3];
    let p3 = [1.3, 0.8];
    let elen_trg = 0.1;
    let ps = sample_uniform_length(
        ControlPoints {
            p0: &p0,
            p1: &p1,
            p2: &p2,
            p3: &p3,
        },
        elen_trg,
        true,
        true,
        30,
    );
    for ip in 0..ps.len() - 1 {
        let q0 = ps[ip];
        let q1 = ps[ip + 1];
        let elen = q0.sub(&q1).norm();
        let dev = (elen - elen_trg).abs();
        // println!("{} {}", ip, dev);
        assert!(dev < 0.007, "{}", dev);
    }
}
