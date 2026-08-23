/// Subtract the mean of the columns (centroid) from each column in place.
pub fn centerize<T>(x: &mut [[T; 3]; 4])
where
    T: num_traits::Float,
{
    let scale = T::one() / (T::one() + T::one() + T::one() + T::one());
    let mov = [
        (x[0][0] + x[1][0] + x[2][0] + x[3][0]) * scale,
        (x[0][1] + x[1][1] + x[2][1] + x[3][1]) * scale,
        (x[0][2] + x[1][2] + x[2][2] + x[3][2]) * scale,
    ];
    for col in x.iter_mut() {
        col[0] = col[0] - mov[0];
        col[1] = col[1] - mov[1];
        col[2] = col[2] - mov[2];
    }
}

pub fn scale_in_place<T>(x: &mut [[T; 3]; 4], s: T)
where
    T: num_traits::Float,
{
    for i in 0..4 {
        for j in 0..3 {
            x[i][j] = x[i][j] * s;
        }
    }
}

pub fn add_scaled<T>(x0: &[[T; 3]; 4], dx: &[[T; 3]; 4], t: T) -> [[T; 3]; 4]
where
    T: num_traits::Float,
{
    let mut x = *x0;
    for i in 0..4 {
        for j in 0..3 {
            x[i][j] = x0[i][j] + t * dx[i][j];
        }
    }
    x
}
