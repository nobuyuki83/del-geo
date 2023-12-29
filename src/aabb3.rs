//! methods for 3D Axis-aligned Bounding Box (AABB)


use num_traits::AsPrimitive;

#[allow(clippy::identity_op)]
pub fn from_vtx2xyz<T>(
    vtx2xyz: &[T],
    eps: T) -> [T;6]
    where T: num_traits::Float
{
    assert!(!vtx2xyz.is_empty());
    let mut aabb = [T::zero();6];
    {
        {
            let cgx = vtx2xyz[0];
            aabb[0] = cgx - eps;
            aabb[3] = cgx + eps;
        }
        {
            let cgy = vtx2xyz[1];
            aabb[1] = cgy - eps;
            aabb[4] = cgy + eps;
        }
        {
            let cgz = vtx2xyz[2];
            aabb[2] = cgz - eps;
            aabb[5] = cgz + eps;
        }
    }
    for i_vtx in 1..vtx2xyz.len()/3 {
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
pub fn from_list_of_vertices<T>(
    idx2vtx: &[usize],
    vtx2xyz: &[T],
    eps: T) -> [T;6]
where T: num_traits::Float
{
    assert!(!idx2vtx.is_empty());
    let mut aabb = [T::zero();6];
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

pub fn center<T>(aabb: &[T;6]) -> [T;3]
where T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>
{
    [
        (aabb[0]+aabb[3])*0.5f64.as_(),
        (aabb[1]+aabb[4])*0.5f64.as_(),
        (aabb[2]+aabb[5])*0.5f64.as_()]
}

pub fn max_edge_size<T>(aabb: &[T;6]) -> T
where T: num_traits::Float
{
    let lx = aabb[3] - aabb[0];
    let ly = aabb[4] - aabb[1];
    let lz = aabb[5] - aabb[2];
    if lx > ly {
        if lx > lz { return lx; }
        else { return lz; }
    }
    if ly > lz {return ly; }
    return lz;
}

#[allow(clippy::identity_op)]
pub fn from_two_aabbs_slice6<T>(
    i0: &[T],
    i1: &[T]) -> [T; 6]
where T: num_traits::Float
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    let mut o = [T::zero(); 6];
    for i in 0..3 {
        o[i + 0] = if i0[i + 0] < i1[i + 0] { i0[i + 0] } else { i1[i + 0] };
        o[i + 3] = if i0[i + 3] > i1[i + 3] { i0[i + 3] } else { i1[i + 3] };
    }
    o
}

pub fn is_active<T>(i0: &[T]) -> bool
where T: std::cmp::PartialOrd
{
    i0[0] <= i0[3]
}

pub fn is_intersect<T>(
    i0: &[T],
    i1: &[T]) -> bool
where T: std::cmp::PartialOrd
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