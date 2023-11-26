pub fn from_list_of_vertices(
    idx2vtx: &[usize],
    vtx2xyz: &[f32],
    eps: f32) -> (nalgebra::Vector3::<f32>, nalgebra::Vector3::<f32>)
{
    assert!(!idx2vtx.is_empty());
    let mut min = nalgebra::Vector3::<f32>::zeros();
    let mut max = nalgebra::Vector3::<f32>::zeros();
    {
        let i_vtx = idx2vtx[0];
        {
            let cgx = vtx2xyz[i_vtx * 3 + 0];
            min.x = cgx - eps;
            max.x = cgx + eps;
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            min.y = cgy - eps;
            max.y = cgy + eps;
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            min.z = cgz - eps;
            max.z = cgz + eps;
        }
    }
    for &i_vtx in idx2vtx.iter().skip(1) {
        {
            let cgx = vtx2xyz[i_vtx * 3 + 0];
            min.x = if cgx - eps < min.x { cgx - eps } else { min.x };
            max.x = if cgx + eps > max.x { cgx + eps } else { max.x };
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            min.y = if cgy - eps < min.y { cgy - eps } else { min.y };
            max.y = if cgy + eps > max.y { cgy + eps } else { max.y };
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            min.z = if cgz - eps < min.z { cgz - eps } else { min.z };
            max.z = if cgz + eps > max.z { cgz + eps } else { max.z };
        }
    }
    assert!(min.x <= max.x);
    assert!(min.y <= max.y);
    assert!(min.z <= max.z);
    (min, max)
}

pub fn from_two_aabbs_slice6(
    i0: &[f32],
    i1: &[f32]) -> [f32; 6]
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    let mut o = [0_f32; 6];
    for i in 0..3 {
        o[i + 0] = if i0[i+0] < i1[i+0] { i0[i+0] } else { i1[i+0] };
        o[i + 3] = if i0[i+3] > i1[i+3] { i0[i+3] } else { i1[i+3] };
    }
    o
}