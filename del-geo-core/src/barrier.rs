pub fn w_energy(mut g: f32, ghat: f32, offset: f32) -> f32 {
    g -= offset;

    if g <= 0.0 {
        f32::INFINITY
    } else if g >= ghat {
        0.0
    } else {
        -(g - ghat) * (g - ghat) * (g / ghat).ln()
    }
}

pub fn dw_energy(mut g: f32, ghat: f32, offset: f32) -> f32 {
    g -= offset;

    if g <= 0.0 {
        f32::NEG_INFINITY
    } else if g >= ghat {
        0.0
    } else {
        (ghat - g) * (2.0 * g * (g / ghat).ln() + g - ghat) / g
    }
}

pub fn ddw_energy(mut g: f32, ghat: f32, offset: f32) -> f32 {
    g -= offset;

    if g <= 0.0 {
        f32::INFINITY
    } else if g >= ghat {
        0.0
    } else {
        -2.0 * (g / ghat).ln() + ghat * (ghat + 2.0 * g) / (g * g) - 3.0
    }
}
