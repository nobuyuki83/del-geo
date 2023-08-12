use num_traits::AsPrimitive;

pub fn resample<T, const X: usize>(
    stroke0: &Vec<nalgebra::base::SVector<T, X>>,
    l: T) -> Vec<nalgebra::base::SVector<T, X>>
    where T: nalgebra::RealField + Copy,
          f64: num_traits::AsPrimitive<T>
{
    if stroke0.len() == 0 {
        return vec!();
    }
    let mut stroke = Vec::<nalgebra::base::SVector<T, X>>::new();
    stroke.push(stroke0[0]);
    let mut jcur = 0;
    let mut rcur: T = 0_f64.as_();
    let mut lcur = l;
    loop {
        if jcur >= stroke0.len() - 1 { break; }
        let lenj = (stroke0[jcur + 1] - stroke0[jcur]).norm();
        let lenjr = lenj * (1_f64.as_() - rcur);
        if lenjr > lcur { // put point in this segment
            rcur += lcur / lenj;
            stroke.push(stroke0[jcur].scale(1_f64.as_() - rcur) + stroke0[jcur + 1].scale(rcur));
            lcur = l;
        } else { // next segment
            lcur -= lenjr;
            rcur = 0_f64.as_();
            jcur += 1;
        }
    }
    stroke
}