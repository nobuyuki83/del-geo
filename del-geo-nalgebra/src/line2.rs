//! methods for 2D line (parameterized by origin and direction vector)

pub fn intersection<T>(
    ps: &nalgebra::Vector2<T>, // point
    pd: &nalgebra::Vector2<T>, // point
    qs: &nalgebra::Vector2<T>, // point
    qd: &nalgebra::Vector2<T>,
) -> nalgebra::Vector2<T>
where
    T: nalgebra::RealField + Copy,
{
    let qn = crate::vec2::rotate90(qd);
    let t = (qs - ps).dot(&qn) / pd.dot(&qn);
    ps + pd.scale(t)
}

pub fn dw_intersection<T>(
    ps: &nalgebra::Vector2<T>, // point
    pd: &nalgebra::Vector2<T>, // point
    qs: &nalgebra::Vector2<T>, // point
    qd: &nalgebra::Vector2<T>,
) -> (
    nalgebra::Vector2<T>,
    nalgebra::Matrix2<T>,
    nalgebra::Matrix2<T>,
)
where
    T: nalgebra::RealField + Copy,
{
    let qn = crate::vec2::rotate90(qd);
    let a = qs - ps;
    let b = T::one() / pd.dot(&qn);
    let t = a.dot(&qn) * b;
    let dt_dqn = a.scale(b) - pd.scale(a.dot(&qn) * b * b);
    let dt_dqd = nalgebra::Vector2::<T>::new(dt_dqn.y, -dt_dqn.x);
    let dt_dqs = qn.scale(b);
    (
        ps + pd.scale(t),
        pd * dt_dqs.transpose(),
        pd * dt_dqd.transpose(),
    )
}

#[test]
fn test_dw_intersection() {
    type Vec2 = nalgebra::Vector2<f64>;
    let ps = Vec2::new(2.5, 0.3);
    let pd = Vec2::new(-2.3, 0.1);
    let qs0 = Vec2::new(1.4, 1.2);
    let qd0 = Vec2::new(-1.5, -2.1);
    let (t0, dtds0, dtdd0) = dw_intersection(&ps, &pd, &qs0, &qd0);
    let eps = 1.0e-5;
    for i_dim in 0..2 {
        let dqs = crate::vec2::basis(i_dim, 1.);
        let qs1 = qs0 + dqs.scale(eps);
        let (t1, _, _) = dw_intersection(&ps, &pd, &qs1, &qd0);
        let v0 = (t1 - t0) / eps;
        let v1 = dtds0 * dqs;
        assert!((v0 - v1).norm() < 1.0e-5);
    }
    for i_dim in 0..2 {
        let dqd = crate::vec2::basis(i_dim, 1.);
        let qd1 = qd0 + crate::vec2::basis(i_dim, eps);
        let (t1, _, _) = dw_intersection(&ps, &pd, &qs0, &qd1);
        let v0 = (t1 - t0) / eps;
        let v1 = dtdd0 * dqd;
        assert!((v0 - v1).norm() < 1.0e-5);
    }
}

pub fn dw_intersection_against_bisector<Real>(
    ls: &nalgebra::Vector2<Real>, // source of line
    ld: &nalgebra::Vector2<Real>, // direction of line
    p0: &nalgebra::Vector2<Real>, // point 0
    p1: &nalgebra::Vector2<Real>,
) -> (
    nalgebra::Vector2<Real>,
    nalgebra::Matrix2<Real>,
    nalgebra::Matrix2<Real>,
)
where
    Real: nalgebra::RealField + Copy,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = Real::one() + Real::one();
    let half = Real::one() / two;
    let a = (p0 + p1) * half;
    let b = crate::vec2::rotate90(&(p1 - p0));
    let (r, drda, drdb) = dw_intersection(ls, ld, &a, &b);
    let dbdp0 = nalgebra::Matrix2::<Real>::new(zero, one, -one, zero);
    let dbdp1 = nalgebra::Matrix2::<Real>::new(zero, -one, one, zero);
    let drdp0 = drda.scale(half) + drdb * dbdp0;
    let drdp1 = drda.scale(half) + drdb * dbdp1;
    (r, drdp0, drdp1)
}

#[test]
fn test_dw_intersection_against_bisector() {
    type Vec2 = nalgebra::Vector2<f64>;
    let ps = Vec2::new(2.5, 0.3);
    let pd = Vec2::new(-2.3, 0.1);
    let qs0 = Vec2::new(1.4, 1.2);
    let qd0 = Vec2::new(-0.5, -2.1);
    let (t0, dtds0, dtdd0) = dw_intersection_against_bisector(&ps, &pd, &qs0, &qd0);
    let eps = 1.0e-5;
    for i_dim in 0..2 {
        let dqs = crate::vec2::basis(i_dim, 1.);
        let qs1 = qs0 + dqs.scale(eps);
        let (t1, _, _) = dw_intersection_against_bisector(&ps, &pd, &qs1, &qd0);
        let v0 = (t1 - t0) / eps;
        let v1 = dtds0 * dqs;
        assert!((v0 - v1).norm() < 1.0e-5);
    }
    for i_dim in 0..2 {
        let dqd = crate::vec2::basis(i_dim, 1.);
        let qd1 = qd0 + crate::vec2::basis(i_dim, eps);
        let (t1, _, _) = dw_intersection_against_bisector(&ps, &pd, &qs0, &qd1);
        let v0 = (t1 - t0) / eps;
        let v1 = dtdd0 * dqd;
        assert!((v0 - v1).norm() < 1.0e-5);
    }
}
