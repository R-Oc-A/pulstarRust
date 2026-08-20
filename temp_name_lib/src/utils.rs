/// This enum helps to encapsulate the errors that may arrise in the profile and pulstar codes.
#[derive(Debug)]
pub enum MathErrors{
    DivisionByZero,
    CosineBiggerThanOne,
    VectorLengthZero,
    DifferentVectorBase,
    OutOfBounds,
    NotAdequateNumberOfElements,
    NotDefinedForRadialOrSectorialPulsations,
    FunctionNotFound,
    CorruptHeaderOfDataFrame,
    OrderOfExpansionNotSupported,
    RequestUnrelatedRotationRegime,
}

pub const MACHINE_PRECISION:f64 = 1.0e-8;