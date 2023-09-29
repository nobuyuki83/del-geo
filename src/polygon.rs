pub fn arclength_polygon<T, const X: usize>(
    vtxs: &Vec<nalgebra::base::SVector<T, X>>) -> T
    where T: nalgebra::RealField + Copy,
          f64: num_traits::AsPrimitive<T>
{
    if vtxs.len() < 2 { return T::zero(); }
    let np = vtxs.len();
    let mut len: T = T::zero();
    for ip0 in 0..np {
        let ip1 = (ip0 + 1) % np;
        len += (vtxs[ip0] - vtxs[ip1]).norm();
    }
    return len;
}

fn match_frames_of_loop_ends(
    vtx2xyz: &nalgebra::Matrix3xX::<f32>,
    vtx2bin0: &nalgebra::Matrix3xX::<f32>) -> nalgebra::Matrix3xX::<f32>
{
    let num_vtx = vtx2xyz.shape().1;
    let theta = {
        let x0 = vtx2bin0.column(0);
        let xn = vtx2bin0.column(num_vtx - 1);
        let vn0 = (vtx2xyz.column(0) - vtx2xyz.column(num_vtx - 1)).normalize();
        let v01 = (vtx2xyz.column(1) - vtx2xyz.column(0)).normalize();
        assert!(x0.dot(&v01).abs() < 1.0e-6);
        let rot = crate::mat3::minimum_rotation_matrix(vn0, v01);
        let x1a = rot * xn;
        let y0 = v01.cross(&x0);
        assert!(x1a.dot(&v01).abs() < 1.0e-6);
        assert!((y0.norm()-1.).abs() < 1.0e-6);
        let c0 = x1a.dot(&x0);
        let s0 = x1a.dot(&y0);
        f32::atan2(s0, c0)
    };
    let theta_step = theta / num_vtx as f32;
    let mut vtx2bin1 = nalgebra::Matrix3xX::<f32>::zeros(num_vtx);
    for iseg in 0..num_vtx {
        let dtheta = theta_step * iseg as f32;
        let x0 = vtx2bin0.column(iseg);
        let ivtx0 = iseg;
        let ivtx1 = (iseg + 1) % num_vtx;
        let v01 = (vtx2xyz.column(ivtx1) - vtx2xyz.column(ivtx0)).normalize();
        let y0 = v01.cross(&x0);
        assert!((x0.cross(&y0).dot(&v01)-1.).abs()<1.0e-5);
        let x1 = x0.scale(dtheta.sin()) + y0.scale(dtheta.cos());
        vtx2bin1.column_mut(iseg).copy_from(&x1);
    }
    vtx2bin1
}

pub fn smooth_frame(vtx2xyz: &nalgebra::Matrix3xX::<f32>) -> nalgebra::Matrix3xX::<f32> {
    let vtx2bin0 = crate::polyline::parallel_transport_polyline(&vtx2xyz);
    match_frames_of_loop_ends(&vtx2xyz, &vtx2bin0)
}

pub fn tube_mesh(
    vtx2xyz: &nalgebra::Matrix3xX::<f32>,
    vtx2bin: &nalgebra::Matrix3xX::<f32>,
    rad: f32) -> (Vec<usize>, Vec<f32>) {
    let n = 8;
    let dtheta = std::f32::consts::PI * 2. / n as f32;
    let num_vtx = vtx2xyz.ncols();
    let mut pnt2xyz = Vec::<f32>::new();
    for ipnt in 0..num_vtx {
        let p0 = vtx2xyz.column(ipnt).into_owned();
        let p1 = vtx2xyz.column((ipnt+1)%num_vtx).into_owned();
        let z0 = (p1-p0).normalize();
        let x0 = vtx2bin.column(ipnt);
        let y0 = z0.cross(&x0);
        for i in 0..n {
            let theta = dtheta*i as f32;
            let v0 =  x0.scale(theta.cos()) + y0.scale(theta.sin());
            let q0 = p0 + v0.scale(rad);
            q0.iter().for_each(|&v| pnt2xyz.push(v));
        }
    }

    let mut tri2pnt = Vec::<usize>::new();
    for iseg in 0..num_vtx {
        let ipnt0 = iseg;
        let ipnt1 = (ipnt0 + 1) % num_vtx;
        for i in 0..n {
            tri2pnt.push(ipnt0 * n + i);
            tri2pnt.push(ipnt0 * n + (i+1)%n);
            tri2pnt.push(ipnt1 * n + i);
            //
            tri2pnt.push(ipnt1 * n + (i+1)%n);
            tri2pnt.push(ipnt1 * n + i);
            tri2pnt.push(ipnt0 * n + (i+1)%n);
        }
    }
    (tri2pnt, pnt2xyz)
}