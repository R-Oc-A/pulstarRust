#[derive(Debug)]
pub enum MathErrors{
    DivisionByZero,
    CosineBiggerThanOne,
    VectorLengthZero,
    DifferentVectorBase,
    OutOfBounds,
    NotAdequateNumberOfElements,
}

pub const MACHINE_PRECISION:f64 = 1.0e-8;