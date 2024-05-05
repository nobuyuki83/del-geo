pub fn intersection<T>(
    ps: &nalgebra::Vector2<T>, // point
    pd: &nalgebra::Vector2<T>, // point
    qs: &nalgebra::Vector2<T>, // point
    qd: &nalgebra::Vector2<T>) -> nalgebra::Vector2<T>
    where T: nalgebra::RealField + Copy
{
    let qn = crate::vec2::rotate90(qd);
    let t = (qs-ps).dot(&qn)/(pd.dot(&qn));
    ps + pd.scale(t)
}