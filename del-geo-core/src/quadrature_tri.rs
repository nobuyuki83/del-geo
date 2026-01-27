pub const DEGREE2IDX: [usize; 4] = [0, 1, 4, 11];

#[allow(clippy::excessive_precision)]
pub const IDX2POSWEIGHT_F32: [[f32; 3]; 11] = [
    // linear
    [0.3333333333, 0.3333333333, 1.0],
    // quadratic
    [0.1666666667, 0.1666666667, 0.3333333333],
    [0.6666666667, 0.1666666667, 0.3333333333],
    [0.1666666667, 0.6666666667, 0.3333333333],
    // cubic
    [0.1012865073, 0.1012865073, 0.1259391805],
    [0.7974269854, 0.1012865073, 0.1259391805],
    [0.1012865073, 0.7974269854, 0.1259391805],
    [0.4701420641, 0.0597158718, 0.1323941527],
    [0.4701420641, 0.4701420641, 0.1323941527],
    [0.0597158718, 0.4701420641, 0.1323941527],
    [0.3333333333, 0.3333333333, 0.225],
];

#[allow(clippy::excessive_precision)]
pub const IDX2POSWEIGHT_F64: [[f64; 3]; 11] = [
    // linear
    [0.3333333333, 0.3333333333, 1.0],
    // quadratic
    [0.1666666667, 0.1666666667, 0.3333333333],
    [0.6666666667, 0.1666666667, 0.3333333333],
    [0.1666666667, 0.6666666667, 0.3333333333],
    // cubic
    [0.1012865073, 0.1012865073, 0.1259391805],
    [0.7974269854, 0.1012865073, 0.1259391805],
    [0.1012865073, 0.7974269854, 0.1259391805],
    [0.4701420641, 0.0597158718, 0.1323941527],
    [0.4701420641, 0.4701420641, 0.1323941527],
    [0.0597158718, 0.4701420641, 0.1323941527],
    [0.3333333333, 0.3333333333, 0.225],
];

pub struct Quad<Real> {
    _marker: std::marker::PhantomData<fn() -> Real>,
}

pub trait QuadratureTri<Real> {
    fn hoge(i_gauss_degree: usize) -> &'static [[Real; 3]];
}

impl<Real> QuadratureTri<f64> for Quad<Real> {
    fn hoge(i_gauss_degree: usize) -> &'static [[f64; 3]] {
        let i0 = DEGREE2IDX[i_gauss_degree];
        let i1 = DEGREE2IDX[i_gauss_degree + 1];
        &IDX2POSWEIGHT_F64[i0..i1]
    }
}

impl<Real> QuadratureTri<f32> for Quad<Real> {
    fn hoge(i_gauss_degree: usize) -> &'static [[f32; 3]] {
        let i0 = DEGREE2IDX[i_gauss_degree];
        let i1 = DEGREE2IDX[i_gauss_degree + 1];
        &IDX2POSWEIGHT_F32[i0..i1]
    }
}
