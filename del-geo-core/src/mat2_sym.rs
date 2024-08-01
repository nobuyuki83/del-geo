//! symmetric matrix [[a,b],[b,c]] parameterized as [a,b,c]

/// ax^2 + 2bxy + cy^2 = 1
pub fn prinsipal_directions<Real>(coeff: &[Real; 3])
                                  -> ([Real; 2], [[Real; 2]; 2])
    where Real: num_traits::Float + std::fmt::Debug
{
    let zero = Real::zero();
    let one = Real::one();
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    if b.is_zero() {
        return ([a, c], [[one,zero],[zero,one]])
    }
    let two = one + one;
    let four = two + two;
    let half = one / two;
    let tmp = ((a - c) * (a - c) + four * b * b).sqrt();
    let lam0 = half * (a + c - tmp);
    let lam1 = half * (a + c + tmp);
    let det0 = a - c + tmp;
    let det1 = a - c - tmp;
    if det0.abs() > det1.abs() {
        let evec0 = [-two * b, det0];
        let evec1 = [det0, two * b];
        ([lam0, lam1], [evec0, evec1])
    } else {
        let evec0 = [det1, two * b];
        let evec1 = [-two * b, det1];
        ([lam0, lam1], [evec0, evec1])
    }
}


pub fn matvec<Real>(coeff: &[Real; 3], v: &[Real; 2]) -> [Real; 2]
    where Real: num_traits::Float
{
    [
        coeff[0] * v[0] + coeff[1] * v[1],
        coeff[1] * v[0] + coeff[2] * v[1]
    ]
}

pub fn aabb2<Real>(coeff: &[Real; 3]) -> [Real; 4]
    where Real: num_traits::Float + std::fmt::Debug
{
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    let det = a * c - b * b;
    let minx = (c / det).sqrt();
    let miny = (a / det).sqrt();
    [-minx, -miny, minx, miny]
}


#[test]
fn test_prinsipal_directions() {
    let mats = [[4., 0., 1.], [4., 0.2, 1.], [1., 0., 4.], [1., -0.1, 4.],[4.,0.,4.],[4.0,1.99,1.0]];
    for mat in mats {
        let ([lam0, lam1], [evec0, evec1]) = prinsipal_directions::<f32>(&mat);
        let evec0 = nalgebra::Vector2::<f32>::from_row_slice(&evec0);
        let evec1 = nalgebra::Vector2::<f32>::from_row_slice(&evec1);
        assert!(evec0.norm() > 1.0e-5);
        assert!(evec1.norm() > 1.0e-5);
        let evec0 = evec0.normalize();
        let evec1 = evec1.normalize();
        assert!(evec0.dot(&evec1).abs()<1.0e-5, "{}", evec0.dot(&evec1));
        {
            let tmp0 = nalgebra::Vector2::<f32>::from_row_slice(
                &matvec(&mat, arrayref::array_ref![evec0.as_slice(),0,2]));
            let diff = (tmp0 - lam0 * evec0).norm();
            assert!(diff < 1.0e-6, "{}", diff);
        }
        {
            let tmp0 = nalgebra::Vector2::<f32>::from_row_slice(
                &matvec(&mat, arrayref::array_ref![evec1.as_slice(),0,2]));
            assert!((tmp0 - lam1 * evec1).norm() < 1.0e-6);
        }
        let evec0 = evec0 / lam0.sqrt();
        let evec1 = evec1 / lam1.sqrt();
        let sig = nalgebra::Matrix2::<f32>::new(mat[0], mat[1], mat[1], mat[2]);
        let aabb = aabb2(&mat);
        let ndiv = 128;
        let mut sdf_max = f32::MIN;
        for i in 0..ndiv {
            let theta = (i as f32) * 2.0 * std::f32::consts::PI / (ndiv as f32);
            let v = evec0 * theta.cos() + evec1 * theta.sin();
            let radrad = (v.transpose() * sig * v)[0];
            assert!((radrad - 1.).abs() < 1.0e-3,"{}",radrad);
            assert!(crate::aabb2::is_inlcude_point(&aabb, &[v[0], v[1]]));
            let sdf = crate::aabb2::sdf(&aabb, &[v[0], v[1]]);
            assert!(sdf<=0f32);
            if sdf > sdf_max { sdf_max = sdf; }
        }
        assert!(sdf_max > -1.0e-3,"{}", sdf_max);
    }
}


