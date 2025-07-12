//! functions for Normalized Device Coordinate (NDC), that is `[-1, 1]^3`.

pub fn sample_inside_uniformly<Reng, T>(rng: &mut Reng) -> [T; 3]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<T>,
{
    let one = T::one();
    let two = one + one;
    std::array::from_fn(|_i| rng.random() * two - one)
}
