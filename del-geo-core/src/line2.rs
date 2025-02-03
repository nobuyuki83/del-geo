//! methods for 2D line (parameterized by origin and direction vector)

pub fn intersection<T>(
    ps: &[T; 2], // point
    pd: &[T; 2], // point
    qs: &[T; 2], // point
    qd: &[T; 2],
) -> [T; 2]
where
    T: num_traits::Float + Copy,
{
    use crate::vec2::Vec2;
    let qn = crate::vec2::rotate90(qd);
    let t = qs.sub(ps).dot(&qn) / pd.dot(&qn);
    ps.add(&pd.scale(t))
}

pub fn dw_intersection<T>(
    ps: &[T; 2], // point
    pd: &[T; 2], // point
    qs: &[T; 2], // point
    qd: &[T; 2],
) -> ([T; 2], [T; 4], [T; 4])
where
    T: num_traits::Float + Copy,
{
    use crate::vec2::Vec2;
    let qn = crate::vec2::rotate90(qd);
    let a = qs.sub(ps);
    let b = T::one() / pd.dot(&qn);
    let t = a.dot(&qn) * b;
    let dt_dqn = a.scale(b).sub(&pd.scale(a.dot(&qn) * b * b));
    let dt_dqd = [dt_dqn[1], -dt_dqn[0]];
    let dt_dqs = qn.scale(b);
    (
        ps.add(&pd.scale(t)),
        crate::mat2_col_major::from_outer_product(pd, &dt_dqs),
        crate::mat2_col_major::from_outer_product(pd, &dt_dqd),
    )
}

#[test]
fn test_dw_intersection() {
    use crate::vec2::Vec2;
    let ps = [2.5, 0.3];
    let pd = [-2.3, 0.1];
    let qs0 = [1.4, 1.2];
    let qd0 = [-1.5, -2.1];
    let (t0, dtds0, dtdd0) = dw_intersection(&ps, &pd, &qs0, &qd0);
    let eps = 1.0e-5;
    for i_dim in 0..2 {
        let dqs = crate::vec2::basis(i_dim, 1.);
        let qs1 = qs0.add(&dqs.scale(eps));
        let (t1, _, _) = dw_intersection(&ps, &pd, &qs1, &qd0);
        let v0 = t1.sub(&t0).scale(1. / eps);
        let v1 = crate::mat2_col_major::mult_vec(&dtds0, &dqs);
        assert!(v0.sub(&v1).norm() < 1.0e-5);
    }
    for i_dim in 0..2 {
        let dqd = crate::vec2::basis(i_dim, 1.);
        let qd1 = qd0.add(&crate::vec2::basis(i_dim, eps));
        let (t1, _, _) = dw_intersection(&ps, &pd, &qs0, &qd1);
        let v0 = t1.sub(&t0).scale(1. / eps);
        let v1 = crate::mat2_col_major::mult_vec(&dtdd0, &dqd);
        assert!(v0.sub(&v1).norm() < 1.0e-5);
    }
}

pub fn dw_intersection_against_bisector<Real>(
    ls: &[Real; 2], // source of line
    ld: &[Real; 2], // direction of line
    p0: &[Real; 2], // point 0
    p1: &[Real; 2],
) -> ([Real; 2], [Real; 4], [Real; 4])
where
    Real: num_traits::Float + Copy,
{
    use crate::mat2_col_major::Mat2ColMajor;
    use crate::vec2::Vec2;
    let zero = Real::zero();
    let one = Real::one();
    let two = Real::one() + Real::one();
    let half = Real::one() / two;
    let a = p0.add(p1).scale(half);
    let b = crate::vec2::rotate90(&p1.sub(p0));
    let (r, drda, drdb) = dw_intersection(ls, ld, &a, &b);
    let dbdp1 = [zero, one, -one, zero];
    let dbdp0 = [zero, -one, one, zero];
    let drdp0 = drda.scale(half).add(&drdb.mult_mat_col_major(&dbdp0));
    let drdp1 = drda.scale(half).add(&drdb.mult_mat_col_major(&dbdp1));
    (r, drdp0, drdp1)
}

#[test]
fn test_dw_intersection_against_bisector() {
    use crate::vec2::Vec2;
    let ps = [2.5, 0.3];
    let pd = [-2.3, 0.1];
    let qs0 = [1.4, 1.2];
    let qd0 = [-0.5, -2.1];
    let (t0, dtds0, dtdd0) = dw_intersection_against_bisector(&ps, &pd, &qs0, &qd0);
    let eps = 1.0e-5;
    for i_dim in 0..2 {
        let dqs = crate::vec2::basis(i_dim, 1.);
        let qs1 = qs0.add(&dqs.scale(eps));
        let (t1, _, _) = dw_intersection_against_bisector(&ps, &pd, &qs1, &qd0);
        let v0 = t1.sub(&t0).scale(1. / eps);
        let v1 = crate::mat2_col_major::mult_vec(&dtds0, &dqs);
        assert!(v0.sub(&v1).norm() < 1.0e-5);
    }
    for i_dim in 0..2 {
        let dqd = crate::vec2::basis(i_dim, 1.);
        let qd1 = qd0.add(&crate::vec2::basis(i_dim, eps));
        let (t1, _, _) = dw_intersection_against_bisector(&ps, &pd, &qs0, &qd1);
        let v0 = t1.sub(&t0).scale(1. / eps);
        let v1 = crate::mat2_col_major::mult_vec(&dtdd0, &dqd);
        assert!(v0.sub(&v1).norm() < 1.0e-5);
    }
}
