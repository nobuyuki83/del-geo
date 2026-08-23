pub fn to_quaternion<Real>(a: &[[Real; 3]; 3]) -> [Real; 4]
where
    Real: num_traits::Float + std::fmt::Debug,
{
    use slice_of_array::SliceFlatExt;
    crate::mat3_col_major::to_quaternion(a.flat().try_into().unwrap())
}
