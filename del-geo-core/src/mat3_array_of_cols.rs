pub fn to_quaternion(a: &[[f32; 3]; 3]) -> [f32; 4] {
    use slice_of_array::SliceFlatExt;
    crate::mat3_col_major::to_quaternion(a.flat().try_into().unwrap())
}
