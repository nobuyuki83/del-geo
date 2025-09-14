/// trait for 4D vector
pub trait Vec4<Real>
where
    Self: Sized,
{
    fn add_in_place(&mut self, other: &Self);
}

impl<Real> Vec4<Real> for [Real; 4]
where
    Real: num_traits::Float,
{
    fn add_in_place(&mut self, other: &Self) {
        self[0] = self[0] + other[0];
        self[1] = self[1] + other[1];
        self[2] = self[2] + other[2];
        self[3] = self[3] + other[3];
    }
}

pub fn distance<Real>(a: &[Real; 4], b: &[Real; 4]) -> Real
where
    Real: num_traits::Float,
{
    a.iter()
        .zip(b.iter())
        .fold(Real::zero(), |acc, (&a, &b)| acc + (a - b) * (a - b))
        .sqrt()
}

pub fn length<Real>(a: &[Real; 4]) -> Real
where
    Real: num_traits::Float,
{
    a.iter().fold(Real::zero(), |acc, &a| acc + a * a).sqrt()
}
