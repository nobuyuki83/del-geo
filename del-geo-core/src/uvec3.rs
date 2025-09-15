//! methods for unit 3D vector

use num_traits::AsPrimitive;
pub trait UVec3<Real> {
    fn map_to_unit2_octahedron(&self) -> [Real; 2];
    fn map_to_unit2_equal_area(&self) -> [Real; 2];
}
impl<Real> UVec3<Real> for [Real; 3]
where
    Real: num_traits::Float + 'static + AsPrimitive<f64>,
    f64: AsPrimitive<Real>,
{
    fn map_to_unit2_octahedron(&self) -> [Real; 2] {
        map_to_unit2_octahedron(self)
    }
    fn map_to_unit2_equal_area(&self) -> [Real; 2] {
        map_to_unit2_equal_area(self)
    }
}

pub fn map_to_unit2_octahedron<Real>(dir: &[Real; 3]) -> [Real; 2]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let half = one / (one + one);
    let zero = Real::zero();
    let n = dir[0].abs() + dir[1].abs() + dir[2].abs();
    let oct = [dir[0] / n, dir[1] / n];
    let oct = if dir[2] < zero {
        [
            (one - oct[1].abs()) * oct[0].signum(),
            (one - oct[0].abs()) * oct[1].signum(),
        ]
    } else {
        oct
    };
    [oct[0] * half + half, oct[1] * half + half]
}

/// mapping the surface of sphere (uvec3) to square unit2 (\[0,1\]^2).
/// (0,0,1) is map to the center (0.5, 0.5),
/// the output is texture coordinate i.e., (0.0, 0.0) is the left-bottom
///
/// <https://github.com/mmp/pbrt-v4/blob/1ae72cfa7344e79a7815a21ed3da746cdccee59b/src/pbrt/util/math.cpp#L317>
pub fn map_to_unit2_equal_area<Real>(d: &[Real; 3]) -> [Real; 2]
where
    Real: num_traits::Float + 'static + AsPrimitive<f64>,
    f64: AsPrimitive<Real>,
{
    let one = Real::one();
    let zero = Real::zero();
    let half = one / (one + one);
    let x = d[0].abs();
    let y = d[1].abs();
    let z = d[2].abs();
    let r = (one - z).sqrt();
    /*
    let phi: Real = {
        // fast approximation of atan2
        let a = x.max(y);
        let b = x.min(y);
        let b = if a == zero { zero } else { b / a };
        let b: f64 = b.as_();
        let t1 = 0.406758566246788489601959989e-5;
        let t2 = 0.636226545274016134946890922156;
        let t3 = 0.61572017898280213493197203466e-2;
        let t4 = -0.247333733281268944196501420480;
        let t5 = 0.881770664775316294736387951347e-1;
        let t6 = 0.419038818029165735901852432784e-1;
        let t7 = -0.251390972343483509333252996350e-1;
        let phi = t1 + b * (t2 + b * (t3 + b * (t4 + b * (t5 + b * (t6 + b * t7))))); // horner
        let phi = phi.as_();
        if x < y { one - phi } else { phi }
    };
    */
    let phi = y.atan2(x);
    let phi = phi * std::f64::consts::FRAC_2_PI.as_();
    let v = phi * r;
    let u = r - v;
    let (u, v) = if d[2] < zero {
        (one - v, one - u)
    } else {
        (u, v)
    };
    let u = u.copysign(d[0]);
    let v = v.copysign(d[1]);
    [u * half + half, v * half + half]
}
