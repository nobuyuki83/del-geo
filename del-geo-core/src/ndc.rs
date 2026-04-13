//! functions for Normalized Device Coordinate (NDC), that is `[-1, 1]^3`.

pub fn sample_inside_uniformly<Reng, T>(rng: &mut Reng) -> [T; 3]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<T>,
{
    use rand::RngExt;
    let one = T::one();
    let two = one + one;
    std::array::from_fn(|_i| rng.random() * two - one)
}

pub fn to_image_coordinate(ndc: &[f32; 3], (img_width, img_height): (usize, usize)) -> [f32; 2] {
    let x0 = (ndc[0] + 1.0) * 0.5;
    let y0 = (-ndc[1] + 1.0) * 0.5;
    [img_width as f32 * x0, img_height as f32 * y0]
}
