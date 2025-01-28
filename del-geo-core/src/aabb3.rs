//! methods for 3D Axis-aligned Bounding Box (AABB)

/// trait for 3D Axis-aligned Bounding Box (AABB)
pub trait AABB3Trait<T> {
    fn sample<Reng: rand::Rng>(&self, reng: &mut Reng) -> [T; 3]
    where
        rand::distr::StandardUniform: rand::distr::Distribution<T>;

    fn center(&self) -> [T; 3];
    fn size(&self) -> [T; 3];
    fn scale(&self, s: T) -> Self;
    fn volume(&self) -> T;
    fn xyz_from_hex_index(&self, i_vtx: usize) -> [T; 3];
    fn max_edge_size(&self) -> T;
    fn is_active(&self) -> bool;
    fn is_intersect(&self, i1: &Self) -> bool;
    fn set_as_cube(&mut self, xyz: &[T; 3], eps: T);
    fn add_point(&mut self, xyz: &[T; 3], eps: T);
    fn is_possible_distance_to_aabb2_smaller_than_threshold(
        &self,
        aabb1: &Self,
        threshold: T,
    ) -> bool;
}
impl<T> AABB3Trait<T> for [T; 6]
where
    T: num_traits::Float,
{
    fn sample<Reng>(&self, reng: &mut Reng) -> [T; 3]
    where
        Reng: rand::Rng,
        T: num_traits::Float,
        rand::distr::StandardUniform: rand::distr::Distribution<T>,
    {
        sample(self, reng)
    }
    fn center(&self) -> [T; 3] {
        center(self)
    }
    fn size(&self) -> [T; 3] {
        size(self)
    }
    fn scale(&self, s: T) -> Self {
        scale(self, s)
    }
    fn volume(&self) -> T {
        volume(self)
    }
    fn xyz_from_hex_index(&self, i_vtx: usize) -> [T; 3] {
        xyz_from_hex_index(self, i_vtx)
    }
    fn max_edge_size(&self) -> T {
        max_edge_size(self)
    }
    fn is_active(&self) -> bool {
        is_active(self)
    }
    fn is_intersect(&self, i1: &Self) -> bool {
        is_intersect(self, i1)
    }
    fn set_as_cube(&mut self, xyz: &[T; 3], eps: T) {
        set_as_cube(self, xyz, eps)
    }
    fn add_point(&mut self, xyz: &[T; 3], eps: T) {
        add_point(self, xyz, eps)
    }
    fn is_possible_distance_to_aabb2_smaller_than_threshold(
        &self,
        aabb1: &Self,
        threshold: T,
    ) -> bool {
        is_possible_distance_to_aabb2_smaller_than_threshold(self, aabb1, threshold)
    }
}

pub fn from_two_aabbs<T>(i0: &[T; 6], i1: &[T; 6]) -> [T; 6]
where
    T: num_traits::Float,
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    let mut o = [T::zero(); 6];
    for i in 0..3 {
        o[i] = if i0[i] < i1[i] { i0[i] } else { i1[i] };
        o[i + 3] = if i0[i + 3] > i1[i + 3] {
            i0[i + 3]
        } else {
            i1[i + 3]
        };
    }
    o
}

// Above: from method
// ----------------------------------
// Below: to method

// ----------------------------------

pub fn scale<T>(aabb: &[T; 6], s: T) -> [T; 6]
where
    T: num_traits::Float,
{
    let c = aabb.center();
    let size = aabb.size();
    let half = T::one() / (T::one() + T::one());
    [
        c[0] - size[0] * s * half,
        c[1] - size[1] * s * half,
        c[2] - size[2] * s * half,
        c[0] + size[0] * s * half,
        c[1] + size[1] * s * half,
        c[2] + size[2] * s * half,
    ]
}

pub fn center<T>(aabb: &[T; 6]) -> [T; 3]
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    [
        (aabb[0] + aabb[3]) * half,
        (aabb[1] + aabb[4]) * half,
        (aabb[2] + aabb[5]) * half,
    ]
}

pub fn size<T>(aabb: &[T; 6]) -> [T; 3]
where
    T: num_traits::Float,
{
    [aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2]]
}

pub fn volume<T>(aabb: &[T; 6]) -> T
where
    T: num_traits::Float,
{
    (aabb[3] - aabb[0]) * (aabb[4] - aabb[1]) * (aabb[5] - aabb[2])
}

pub fn xyz_from_hex_index<Real>(aabb: &[Real; 6], i_vtx: usize) -> [Real; 3]
where
    Real: num_traits::Float,
{
    match i_vtx {
        0 => [aabb[0], aabb[1], aabb[2]],
        1 => [aabb[3], aabb[1], aabb[2]],
        2 => [aabb[3], aabb[4], aabb[2]],
        3 => [aabb[0], aabb[4], aabb[2]],
        4 => [aabb[0], aabb[1], aabb[5]],
        5 => [aabb[3], aabb[1], aabb[5]],
        6 => [aabb[3], aabb[4], aabb[5]],
        7 => [aabb[0], aabb[4], aabb[5]],
        _ => panic!(),
    }
}

pub fn max_edge_size<T>(aabb: &[T; 6]) -> T
where
    T: num_traits::Float,
{
    let lx = aabb[3] - aabb[0];
    let ly = aabb[4] - aabb[1];
    let lz = aabb[5] - aabb[2];
    if lx > ly {
        return if lx > lz { lx } else { lz };
    }
    if ly > lz {
        return ly;
    }
    lz
}

pub fn is_active<T>(i0: &[T; 6]) -> bool
where
    T: PartialOrd,
{
    i0[0] <= i0[3]
}

pub fn is_intersect<T>(i0: &[T; 6], i1: &[T; 6]) -> bool
where
    T: PartialOrd,
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    if !is_active(i0) {
        return false;
    }
    if !is_active(i1) {
        return false;
    }
    if i0[0] > i1[3] {
        return false;
    }
    if i0[1] > i1[4] {
        return false;
    }
    if i0[2] > i1[5] {
        return false;
    }
    if i0[3] < i1[0] {
        return false;
    }
    if i0[4] < i1[1] {
        return false;
    }
    if i0[5] < i1[2] {
        return false;
    }
    true
}

/// return a vec3 sampled inside a aabb
pub fn sample<Reng, T>(aabb: &[T; 6], reng: &mut Reng) -> [T; 3]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<T>,
{
    let r0 = reng.random::<T>();
    let r1 = reng.random::<T>();
    let r2 = reng.random::<T>();
    [
        aabb[0] + r0 * (aabb[3] - aabb[0]),
        aabb[1] + r1 * (aabb[4] - aabb[1]),
        aabb[2] + r2 * (aabb[5] - aabb[2]),
    ]
}

pub fn set_as_cube<Real>(aabb: &mut [Real; 6], xyz: &[Real; 3], eps: Real)
where
    Real: num_traits::Float,
{
    {
        let cgx = xyz[0];
        aabb[0] = cgx - eps;
        aabb[3] = cgx + eps;
    }
    {
        let cgy = xyz[1];
        aabb[1] = cgy - eps;
        aabb[4] = cgy + eps;
    }
    {
        let cgz = xyz[2];
        aabb[2] = cgz - eps;
        aabb[5] = cgz + eps;
    }
}

pub fn add_point<Real>(aabb: &mut [Real; 6], xyz: &[Real; 3], eps: Real)
where
    Real: num_traits::Float,
{
    aabb[0] = aabb[0].min(xyz[0] - eps);
    aabb[3] = aabb[3].max(xyz[0] + eps);
    aabb[1] = aabb[1].min(xyz[1] - eps);
    aabb[4] = aabb[4].max(xyz[1] + eps);
    aabb[2] = aabb[2].min(xyz[2] - eps);
    aabb[5] = aabb[5].max(xyz[2] + eps);
}

pub fn is_possible_distance_to_aabb2_smaller_than_threshold<Real>(
    aabb0: &[Real; 6],
    aabb1: &[Real; 6],
    threshold: Real,
) -> bool
where
    Real: num_traits::Float,
{
    for i_dim in 0..3 {
        let range0 = (aabb0[i_dim], aabb0[i_dim + 3]);
        let range1 = (aabb1[i_dim], aabb1[i_dim + 3]);
        if let Some(dist) = crate::range::distance_to_range(range0, range1) {
            if dist > threshold {
                return false;
            }
        }
    }
    true
}

// --------------------------

pub type AABB3<Real> = crate::aabb::AABB<Real, 3, 6>;

pub fn from_slice<Real>(s: &[Real]) -> AABB3<Real>
where
    Real: num_traits::Float,
{
    AABB3 {
        aabb: s[..6].try_into().unwrap(),
    }
}

pub fn from_aabbs<Real>(aabbs: &[Real], i_aabb: usize) -> AABB3<Real>
where
    Real: num_traits::Float,
{
    AABB3 {
        aabb: aabbs[(i_aabb * 6)..((i_aabb * 6) + 6)].try_into().unwrap(),
    }
}
