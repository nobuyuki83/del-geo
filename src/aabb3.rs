//! methods for 3D Axis-aligned Bounding Box (AABB)

#[allow(clippy::identity_op)]
pub fn from_list_of_vertices(
    idx2vtx: &[usize],
    vtx2xyz: &[f32],
    eps: f32) -> [f32;6]
{
    assert!(!idx2vtx.is_empty());
    let mut aabb = [0_f32;6];
    {
        let i_vtx = idx2vtx[0];
        {
            let cgx = vtx2xyz[i_vtx * 3 + 0];
            aabb[0] = cgx - eps;
            aabb[3] = cgx + eps;
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            aabb[1] = cgy - eps;
            aabb[4] = cgy + eps;
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            aabb[2] = cgz - eps;
            aabb[5] = cgz + eps;
        }
    }
    for &i_vtx in idx2vtx.iter().skip(1) {
        {
            let cgx = vtx2xyz[i_vtx * 3 + 0];
            aabb[0] = if cgx - eps < aabb[0] { cgx - eps } else { aabb[0] };
            aabb[3] = if cgx + eps > aabb[3] { cgx + eps } else { aabb[3]};
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            aabb[1] = if cgy - eps < aabb[1] { cgy - eps } else { aabb[1] };
            aabb[4] = if cgy + eps > aabb[4] { cgy + eps } else { aabb[4] };
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            aabb[2] = if cgz - eps < aabb[2] { cgz - eps } else { aabb[2] };
            aabb[5] = if cgz + eps > aabb[5] { cgz + eps } else { aabb[5] };
        }
    }
    assert!(aabb[0] <= aabb[3]);
    assert!(aabb[1] <= aabb[4]);
    assert!(aabb[2] <= aabb[5]);
    aabb
}

#[allow(clippy::identity_op)]
pub fn from_two_aabbs_slice6(
    i0: &[f32],
    i1: &[f32]) -> [f32; 6]
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    let mut o = [0_f32; 6];
    for i in 0..3 {
        o[i + 0] = if i0[i + 0] < i1[i + 0] { i0[i + 0] } else { i1[i + 0] };
        o[i + 3] = if i0[i + 3] > i1[i + 3] { i0[i + 3] } else { i1[i + 3] };
    }
    o
}

pub fn is_active(i0: &[f32]) -> bool
{
    i0[0] <= i0[3]
}

pub fn is_intersect(
    i0: &[f32],
    i1: &[f32]) -> bool
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    if !is_active(i0) { return false; }
    if !is_active(i1) { return false; }
    if i0[0] > i1[3] { return false; }
    if i0[1] > i1[4] { return false; }
    if i0[2] > i1[5] { return false; }
    if i0[3] < i1[0] { return false; }
    if i0[4] < i1[1] { return false; }
    if i0[5] < i1[2] { return false; }
    true
}