

fn squared_norm<Real>(sm: &[Real; 6]) -> Real
    where Real: num_traits::Float
{
    let two = Real::one() + Real::one();
    sm[0] * sm[0] + sm[1] * sm[1] + sm[2] * sm[2] + two * (sm[3] * sm[3] + sm[4] * sm[4] + sm[5] * sm[5])
}

pub fn eigen_decomp(
    sm: [f64; 6],
    nitr: usize) -> Option<([f64; 9], [f64; 3])>
{
    let mut u = [0f64; 9];
    let mut l = [0f64; 3];
    // initialize u as identity matrix
    u[0] = 1.0;
    u[4] = 1.0;
    u[8] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    u[5] = 0.0;
    u[6] = 0.0;
    u[7] = 0.0;
    l[0] = 0.0;
    l[1] = 0.0;
    l[2] = 0.0;
    let dnrm = squared_norm(&sm);
    if dnrm < 1.0e-30 { return None; } // this matrix is too small
    let scale = dnrm.sqrt();
    let invscl = 1.0 / scale;
    let mut sms = [sm[0] * invscl, sm[1] * invscl, sm[2] * invscl, sm[3] * invscl, sm[4] * invscl, sm[5] * invscl];
    for itr in 0..nitr {
        let m = [sms[0], sms[1], sms[2], sms[3], sms[4], sms[5]];
        let v = [u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]];
        let a12 = sms[3].abs();
        let a20 = sms[4].abs();
        let a01 = sms[5].abs();
        if a12 >= a20 && a12 >= a01 {
            // when a12 sms[3] is the biggest
            let t = 0.5 * (2. * m[3]).atan2(m[2] - m[1]);
            let ct = t.cos();
            let st = t.sin();
            sms[1] = ct * ct * m[1] + st * st * m[2] - 2. * st * ct * m[3];
            sms[2] = ct * ct * m[2] + st * st * m[1] + 2. * st * ct * m[3];
            sms[3] = 0.; // (ct*ct-st*st)*m[3]+st*ct*(m[1]-m[2]);
            sms[4] = st * m[5] + ct * m[4];
            sms[5] = ct * m[5] - st * m[4];
            //
            u[1] = ct * v[1] - st * v[2];
            u[2] = st * v[1] + ct * v[2];
            u[4] = ct * v[4] - st * v[5];
            u[5] = st * v[4] + ct * v[5];
            u[7] = ct * v[7] - st * v[8];
            u[8] = st * v[7] + ct * v[8];
        } else if a20 >= a01 && a20 >= a12 {
            // when a20 sms[4] is the biggest
            // the above condition statement shoud pass exactly once for each iteration.
            let t = 0.5 * (2. * m[4]).atan2(m[2] - m[0]);
            let ct = t.cos();
            let st = t.sin();
            sms[0] = ct * ct * m[0] + st * st * m[2] - 2. * st * ct * m[4];
            sms[2] = ct * ct * m[2] + st * st * m[0] + 2. * st * ct * m[4];
            sms[3] = st * m[5] + ct * m[3];
            sms[4] = 0.; // (ct*ct-st*st)*m[4]+st*ct*(m[0]-m[2]);
            sms[5] = ct * m[5] - st * m[3];
            //
            u[0] = ct * v[0] - st * v[2];
            u[2] = st * v[0] + ct * v[2];
            u[3] = ct * v[3] - st * v[5];
            u[5] = st * v[3] + ct * v[5];
            u[6] = ct * v[6] - st * v[8];
            u[8] = st * v[6] + ct * v[8];
        } else {
            // when a01 sms[5] is the biggest
            // the condition statement shoud pass exactly once for each iteration.
            let t = 0.5 * (2. * m[5]).atan2(m[1] - m[0]);
            let ct = t.cos();
            let st = t.sin();
            sms[0] = ct * ct * m[0] + st * st * m[1] - 2. * st * ct * m[5];
            sms[1] = ct * ct * m[1] + st * st * m[0] + 2. * st * ct * m[5];
            sms[3] = st * m[4] + ct * m[3];
            sms[4] = ct * m[4] - st * m[3];
            sms[5] = 0.; // (ct*ct-st*st)*m[5]+st*ct*(m[0]-m[1]);
            //
            u[0] = ct * v[0] - st * v[1];
            u[1] = st * v[0] + ct * v[1];
            u[3] = ct * v[3] - st * v[4];
            u[4] = st * v[3] + ct * v[4];
            u[6] = ct * v[6] - st * v[7];
            u[7] = st * v[6] + ct * v[7];
        }
    }
    l[0] = scale * sms[0];
    l[1] = scale * sms[1];
    l[2] = scale * sms[2];
    Some((u, l))
}

#[test]
fn hoge() {
    use rand::SeedableRng;
    use rand::Rng;
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    // std::uniform_real_distribution < double > dist(-50.0, 50.0);
    for itr in 0..1000 {
        let sm = {
            let mut sm = [0f64; 6];
            for i in 0..6 {
                sm[i] = rng.gen::<f64>() * 50.;
            }
            sm
        };
        let Some((u,l)) = eigen_decomp(
            sm, 20) else { todo!() };
        {
            let ut = crate::mat3_row_major::transpose(&u);
            let utu = crate::mat3_row_major::mult_mat_row_major(&u, &ut);
            let id = crate::mat3_row_major::from_identity();
            let diff = crate::mat3_row_major::sub(&id, &utu);
            let diffnorm = crate::mat3_row_major::squared_norm(&diff);
            assert!( diffnorm < 1.0e-20 );
        }
    }
}