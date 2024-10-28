#![no_std]

extern crate alloc;

pub mod arith;
mod fields;
mod groups;

use crate::fields::FieldElement;
use crate::groups::{G1Params, G2Params, GroupElement, GroupParams};

use alloc::vec::Vec;
use core::fmt::Display;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_bigint::BigUint;
use rand::Rng;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fr(fields::Fr);

impl Fr {
    pub fn zero() -> Self {
        Fr(fields::Fr::zero())
    }
    pub fn one() -> Self {
        Fr(fields::Fr::one())
    }
    pub fn random<R: Rng>(rng: &mut R) -> Self {
        Fr(fields::Fr::random(rng))
    }
    pub fn pow(&self, exp: Fr) -> Self {
        Fr(self.0.pow(exp.0))
    }
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        fields::Fr::from_str(s).map(Fr)
    }
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Fr)
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn interpret(buf: &[u8; 64]) -> Fr {
        Fr(fields::Fr::interpret(buf))
    }
    pub fn from_slice(slice: &[u8]) -> Result<Self, FieldError> {
        arith::U256::from_slice(slice)
            .map_err(|_| FieldError::InvalidSliceLength) // todo: maybe more sensful error handling
            .map(Fr::new_mul_factor)
    }
    pub fn from_bytes_be_mod_order(slice: &[u8]) -> Result<Self, FieldError> {
        let mut modulus_bytes = [0u8; 32];
        Fr::modulus().to_big_endian(&mut modulus_bytes).unwrap();
        let modulus = BigUint::from_bytes_be(&modulus_bytes);

        let num = BigUint::from_bytes_be(slice) % modulus;

        Fr::from_slice(&num.to_bytes_be())
    }
    pub fn to_big_endian(&self, slice: &mut [u8]) -> Result<(), FieldError> {
        // NOTE: serialized in Montgomery form (as in the original bn crate)
        self.0
            .to_mont()
            .to_big_endian(slice)
            .map_err(|_| FieldError::InvalidSliceLength)
    }
    pub fn new(val: arith::U256) -> Option<Self> {
        fields::Fr::new(val).map(Fr)
    }
    pub fn new_mul_factor(val: arith::U256) -> Self {
        Fr(fields::Fr::new_mul_factor(val))
    }
    pub fn modulus() -> arith::U256 {
        fields::Fr::modulus()
    }
    pub fn into_u256(self) -> arith::U256 {
        (self.0).into()
    }
    pub fn set_bit(&mut self, bit: usize, to: bool) {
        self.0.set_bit(bit, to);
    }
}

impl Add<Fr> for Fr {
    type Output = Fr;

    fn add(self, other: Fr) -> Fr {
        Fr(self.0 + other.0)
    }
}

impl Sub<Fr> for Fr {
    type Output = Fr;

    fn sub(self, other: Fr) -> Fr {
        Fr(self.0 - other.0)
    }
}

impl Neg for Fr {
    type Output = Fr;

    fn neg(self) -> Fr {
        Fr(-self.0)
    }
}

impl Mul for Fr {
    type Output = Fr;

    fn mul(self, other: Fr) -> Fr {
        Fr(self.0 * other.0)
    }
}

impl Div for Fr {
    type Output = Fr;

    fn div(self, other: Fr) -> Fr {
        Fr(self.0 / other.0)
    }
}

impl AddAssign<Fr> for Fr {
    fn add_assign(&mut self, other: Fr) {
        *self = *self + other;
    }
}

impl SubAssign<Fr> for Fr {
    fn sub_assign(&mut self, other: Fr) {
        *self = *self - other;
    }
}

impl MulAssign<Fr> for Fr {
    fn mul_assign(&mut self, other: Fr) {
        *self = *self * other;
    }
}

impl DivAssign<Fr> for Fr {
    fn div_assign(&mut self, other: Fr) {
        *self = *self / other;
    }
}

#[derive(Debug)]
pub enum FieldError {
    InvalidSliceLength,
    InvalidU512Encoding,
    NotMember,
}

impl Display for FieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            FieldError::InvalidSliceLength => write!(f, "Invalid slice length"),
            FieldError::InvalidU512Encoding => write!(f, "Invalid U512 encoding"),
            FieldError::NotMember => write!(f, "Not a member of the field"),
        }
    }
}

#[derive(Debug)]
pub enum CurveError {
    InvalidEncoding,
    NotMember,
    Field(FieldError),
    ToAffineConversion,
}

impl Display for CurveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CurveError::InvalidEncoding => write!(f, "Invalid encoding"),
            CurveError::NotMember => write!(f, "Not a member of the curve"),
            CurveError::Field(fe) => write!(f, "Field error: {:?}", fe),
            CurveError::ToAffineConversion => write!(f, "Failed to convert to affine coordinates"),
        }
    }
}

impl From<FieldError> for CurveError {
    fn from(fe: FieldError) -> Self {
        CurveError::Field(fe)
    }
}

pub use crate::groups::Error as GroupError;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd)]
#[repr(C)]
pub struct Fq(pub fields::Fq);

impl Fq {
    pub fn zero() -> Self {
        Fq(fields::Fq::zero())
    }
    pub fn one() -> Self {
        Fq(fields::Fq::one())
    }
    pub fn random<R: Rng>(rng: &mut R) -> Self {
        Fq(fields::Fq::random(rng))
    }
    pub fn pow(&self, exp: Fq) -> Self {
        Fq(self.0.pow(exp.0))
    }
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        fields::Fq::from_str(s).map(Fq)
    }
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Fq)
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn interpret(buf: &[u8; 64]) -> Fq {
        Fq(fields::Fq::interpret(buf))
    }
    pub fn from_slice(slice: &[u8]) -> Result<Self, FieldError> {
        arith::U256::from_slice(slice)
            .map_err(|_| FieldError::InvalidSliceLength) // todo: maybe more sensful error handling
            .and_then(|x| fields::Fq::new(x).ok_or(FieldError::NotMember))
            .map(Fq)
    }
    pub fn from_be_bytes_mod_order(bytes: &[u8]) -> Result<Self, FieldError> {
        let mut modulus_bytes = [0u8; 32];
        Fq::modulus().to_big_endian(&mut modulus_bytes).unwrap();
        let modulus = BigUint::from_bytes_be(&modulus_bytes);

        let num = BigUint::from_bytes_be(bytes) % modulus;

        Fq::from_slice(&num.to_bytes_be())
    }

    pub fn to_mont_big_endian(&self, slice: &mut [u8]) -> Result<(), FieldError> {
        let a: arith::U256 = self.0.to_mont().into();
        a.to_big_endian(slice)
            .map_err(|_| FieldError::InvalidSliceLength)
    }

    pub fn to_big_endian(&self, slice: &mut [u8]) -> Result<(), FieldError> {
        let a: arith::U256 = self.0.into();
        a.to_big_endian(slice)
            .map_err(|_| FieldError::InvalidSliceLength)
    }

    pub fn from_u256(u256: arith::U256) -> Result<Self, FieldError> {
        Ok(Fq(fields::Fq::new(u256).ok_or(FieldError::NotMember)?))
    }
    pub fn into_u256(self) -> arith::U256 {
        (self.0).into()
    }
    pub fn modulus() -> arith::U256 {
        fields::Fq::modulus()
    }

    pub fn sqrt(&self) -> Option<Self> {
        self.0.sqrt().map(Fq)
    }
}

impl Add<Fq> for Fq {
    type Output = Fq;

    fn add(self, other: Fq) -> Fq {
        Fq(self.0 + other.0)
    }
}

impl Sub<Fq> for Fq {
    type Output = Fq;

    fn sub(self, other: Fq) -> Fq {
        Fq(self.0 - other.0)
    }
}

impl Neg for Fq {
    type Output = Fq;

    fn neg(self) -> Fq {
        Fq(-self.0)
    }
}

impl Mul for Fq {
    type Output = Fq;

    fn mul(self, other: Fq) -> Fq {
        Fq(self.0 * other.0)
    }
}

impl Ord for Fq {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq2(pub fields::Fq2);

impl Fq2 {
    pub fn one() -> Fq2 {
        Fq2(fields::Fq2::one())
    }

    pub fn i() -> Fq2 {
        Fq2(fields::Fq2::i())
    }

    pub fn zero() -> Fq2 {
        Fq2(fields::Fq2::zero())
    }

    /// Initalizes new F_q2(a + bi, a is real coeff, b is imaginary)
    pub fn new(a: Fq, b: Fq) -> Fq2 {
        Fq2(fields::Fq2::new(a.0, b.0))
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    pub fn pow(&self, exp: arith::U256) -> Self {
        Fq2(self.0.pow(exp))
    }

    pub fn real(&self) -> Fq {
        Fq(*self.0.real())
    }

    pub fn imaginary(&self) -> Fq {
        Fq(*self.0.imaginary())
    }

    pub fn sqrt(&self) -> Option<Self> {
        self.0.sqrt().map(Fq2)
    }

    pub fn from_slice(bytes: &[u8]) -> Result<Self, FieldError> {
        let u512 = arith::U512::from_slice(bytes).map_err(|_| FieldError::InvalidU512Encoding)?;
        let (res, c0) = u512.divrem(&Fq::modulus());
        Ok(Fq2::new(
            Fq::from_u256(c0).map_err(|_| FieldError::NotMember)?,
            Fq::from_u256(res.ok_or(FieldError::NotMember)?).map_err(|_| FieldError::NotMember)?,
        ))
    }
}

impl Add<Fq2> for Fq2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Fq2(self.0 + other.0)
    }
}

impl Sub<Fq2> for Fq2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Fq2(self.0 - other.0)
    }
}

impl Neg for Fq2 {
    type Output = Self;

    fn neg(self) -> Self {
        Fq2(-self.0)
    }
}

impl Mul for Fq2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Fq2(self.0 * other.0)
    }
}

impl Div for Fq2 {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Fq2(self.0 / other.0)
    }
}

pub trait Group:
    Send
    + Sync
    + Copy
    + Clone
    + PartialEq
    + Eq
    + Sized
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Fr, Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
    fn is_zero(&self) -> bool;
    fn normalize(&mut self);
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct G1(groups::G1);

impl G1 {
    pub fn new(x: Fq, y: Fq, z: Fq) -> Self {
        G1(groups::G1::new(x.0, y.0, z.0))
    }

    pub fn zero() -> Self {
        G1(groups::G1::zero())
    }

    pub fn x(&self) -> Fq {
        Fq(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq {
        Fq(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq) {
        *self.0.y_mut() = y.0
    }

    pub fn z(&self) -> Fq {
        Fq(*self.0.z())
    }

    pub fn set_z(&mut self, z: Fq) {
        *self.0.z_mut() = z.0
    }

    pub fn b() -> Fq {
        Fq(G1Params::coeff_b())
    }

    pub fn from_compressed(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 33 {
            return Err(CurveError::InvalidEncoding);
        }

        let sign = bytes[0];
        let fq = Fq::from_slice(&bytes[1..])?;
        let x = fq;
        let y_squared = (fq * fq * fq) + Self::b();

        let mut y = y_squared.sqrt().ok_or(CurveError::NotMember)?;

        if (sign == 2 && y.into_u256().get_bit(0).expect("bit 0 always exist; qed"))
            || (sign == 3 && !y.into_u256().get_bit(0).expect("bit 0 always exist; qed"))
        {
            y = y.neg();
        } else if sign != 3 && sign != 2 {
            return Err(CurveError::InvalidEncoding);
        }
        AffineG1::new(x, y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }

    pub fn msm(points: &[Self], scalars: &[Fr]) -> Self {
        G1(groups::G1::msm_variable_base(
            &points.iter().map(|p| p.0).collect::<Vec<_>>(),
            &scalars.iter().map(|x| x.0).collect::<Vec<_>>(),
        ))
    }

    pub fn double(&self) -> Self {
        G1(self.0.double())
    }
}

impl Group for G1 {
    fn zero() -> Self {
        G1(groups::G1::zero())
    }
    fn one() -> Self {
        G1(groups::G1::one())
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        G1(groups::G1::random(rng))
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn normalize(&mut self) {
        let new = match self.0.to_affine() {
            Some(a) => a,
            None => return,
        };

        self.0 = new.to_jacobian();
    }
}

impl Add<G1> for G1 {
    type Output = G1;

    fn add(self, other: G1) -> G1 {
        G1(self.0 + other.0)
    }
}

impl Sub<G1> for G1 {
    type Output = G1;

    fn sub(self, other: G1) -> G1 {
        G1(self.0 - other.0)
    }
}

impl Neg for G1 {
    type Output = G1;

    fn neg(self) -> G1 {
        G1(-self.0)
    }
}

impl Mul<Fr> for G1 {
    type Output = G1;

    fn mul(self, other: Fr) -> G1 {
        G1(self.0 * other.0)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct AffineG1(groups::AffineG1);

impl Default for AffineG1 {
    fn default() -> Self {
        AffineG1(groups::AffineG::one())
    }
}

impl AffineG1 {
    pub fn new(x: Fq, y: Fq) -> Result<Self, GroupError> {
        Ok(AffineG1(groups::AffineG1::new(x.0, y.0)?))
    }

    pub fn new_unchecked(x: Fq, y: Fq) -> Self {
        AffineG1(groups::AffineG1::new_unchecked(x.0, y.0))
    }

    pub fn zero() -> Self {
        AffineG1(groups::AffineG1::zero())
    }

    pub fn one() -> Self {
        AffineG1(groups::AffineG1::one())
    }

    pub fn x(&self) -> Fq {
        Fq(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq {
        Fq(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq) {
        *self.0.y_mut() = y.0
    }

    pub fn from_jacobian(g1: G1) -> Option<Self> {
        g1.0.to_affine().map(AffineG1)
    }

    pub fn get_ys_from_x_unchecked(x: Fq) -> Option<(Fq, Fq)> {
        groups::AffineG1::get_ys_from_x_unchecked(x.0).map(|(neq_y, y)| (Fq(neq_y), Fq(y)))
    }

    pub fn msm(points: &[Self], scalars: &[Fr]) -> Self {
        AffineG1(groups::AffineG1::msm_variable_base(
            &points.iter().map(|p| p.0).collect::<Vec<_>>(),
            &scalars.iter().map(|x| x.0).collect::<Vec<_>>(),
        ))
    }
}

impl Neg for AffineG1 {
    type Output = AffineG1;

    fn neg(self) -> AffineG1 {
        AffineG1(-self.0)
    }
}

impl Into<G1> for AffineG1 {
    fn into(self) -> G1 {
        G1(self.0.to_jacobian())
    }
}

impl Into<AffineG1> for G1 {
    fn into(self) -> AffineG1 {
        AffineG1(
            self.0
                .to_affine()
                .expect("Unable to convert G1 to AffineG1"),
        )
    }
}

impl Add<AffineG1> for AffineG1 {
    type Output = AffineG1;

    fn add(self, other: AffineG1) -> AffineG1 {
        AffineG1(self.0 + other.0)
    }
}

impl Sub<AffineG1> for AffineG1 {
    type Output = AffineG1;

    fn sub(self, other: AffineG1) -> AffineG1 {
        AffineG1(self.0 - other.0)
    }
}

impl Mul<Fr> for AffineG1 {
    type Output = AffineG1;

    fn mul(self, other: Fr) -> AffineG1 {
        AffineG1(self.0 * other.0)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct G2(groups::G2);

impl G2 {
    pub fn new(x: Fq2, y: Fq2, z: Fq2) -> Self {
        G2(groups::G2::new(x.0, y.0, z.0))
    }

    pub fn x(&self) -> Fq2 {
        Fq2(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq2) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq2 {
        Fq2(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq2) {
        *self.0.y_mut() = y.0
    }

    pub fn z(&self) -> Fq2 {
        Fq2(*self.0.z())
    }

    pub fn set_z(&mut self, z: Fq2) {
        *self.0.z_mut() = z.0
    }

    pub fn b() -> Fq2 {
        Fq2(G2Params::coeff_b())
    }

    pub fn from_compressed(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 65 {
            return Err(CurveError::InvalidEncoding);
        }

        let sign = bytes[0];
        let x = Fq2::from_slice(&bytes[1..])?;

        let y_squared = (x * x * x) + G2::b();
        let y = y_squared.sqrt().ok_or(CurveError::NotMember)?;
        let y_neg = -y;
        let y_gt = y.0.to_u512() > y_neg.0.to_u512();

        let e_y = if sign == 10 {
            if y_gt {
                y_neg
            } else {
                y
            }
        } else if sign == 11 {
            if y_gt {
                y
            } else {
                y_neg
            }
        } else {
            return Err(CurveError::InvalidEncoding);
        };

        AffineG2::new(x, e_y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }
}

impl Group for G2 {
    fn zero() -> Self {
        G2(groups::G2::zero())
    }
    fn one() -> Self {
        G2(groups::G2::one())
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        G2(groups::G2::random(rng))
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn normalize(&mut self) {
        let new = match self.0.to_affine() {
            Some(a) => a,
            None => return,
        };

        self.0 = new.to_jacobian();
    }
}

impl Add<G2> for G2 {
    type Output = G2;

    fn add(self, other: G2) -> G2 {
        G2(self.0 + other.0)
    }
}

impl Sub<G2> for G2 {
    type Output = G2;

    fn sub(self, other: G2) -> G2 {
        G2(self.0 - other.0)
    }
}

impl Neg for G2 {
    type Output = G2;

    fn neg(self) -> G2 {
        G2(-self.0)
    }
}

impl Mul<Fr> for G2 {
    type Output = G2;

    fn mul(self, other: Fr) -> G2 {
        G2(self.0 * other.0)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Gt(fields::Fq12);

impl Gt {
    pub fn zero() -> Self {
        Gt(fields::Fq12::zero())
    }
    pub fn one() -> Self {
        Gt(fields::Fq12::one())
    }
    pub fn pow(&self, exp: Fr) -> Self {
        Gt(self.0.pow(exp.0))
    }
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Gt)
    }
    pub fn final_exponentiation(&self) -> Option<Self> {
        self.0.final_exponentiation().map(Gt)
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn is_one(&self) -> bool {
        self == &Gt::one()
    }
}

impl Mul<Gt> for Gt {
    type Output = Gt;

    fn mul(self, other: Gt) -> Gt {
        Gt(self.0 * other.0)
    }
}

pub fn pairing(p: G1, q: G2) -> Gt {
    Gt(groups::pairing(&p.0, &q.0))
}

pub fn pairing_batch(pairs: &[(G1, G2)]) -> Gt {
    let mut ps: Vec<groups::G1> = Vec::new();
    let mut qs: Vec<groups::G2> = Vec::new();
    for (p, q) in pairs {
        ps.push(p.0);
        qs.push(q.0);
    }
    Gt(groups::pairing_batch(&ps, &qs))
}

pub fn miller_loop_batch(pairs: &[(G2, G1)]) -> Result<Gt, CurveError> {
    let mut ps: Vec<groups::G2Precomp> = Vec::new();
    let mut qs: Vec<groups::AffineG<groups::G1Params>> = Vec::new();
    for (p, q) in pairs {
        ps.push(
            p.0.to_affine()
                .ok_or(CurveError::ToAffineConversion)?
                .precompute(),
        );
        qs.push(q.0.to_affine().ok_or(CurveError::ToAffineConversion)?);
    }
    Ok(Gt(groups::miller_loop_batch(&ps, &qs)))
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
#[repr(C)]
pub struct AffineG2(groups::AffineG2);

impl AffineG2 {
    pub fn zero() -> Self {
        AffineG2(groups::AffineG2::zero())
    }

    pub fn one() -> Self {
        AffineG2(groups::AffineG2::one())
    }

    pub fn new(x: Fq2, y: Fq2) -> Result<Self, GroupError> {
        Ok(AffineG2(groups::AffineG2::new(x.0, y.0)?))
    }

    pub fn new_unchecked(x: Fq2, y: Fq2) -> Self {
        AffineG2(groups::AffineG2::new_unchecked(x.0, y.0))
    }

    pub fn x(&self) -> Fq2 {
        Fq2(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq2) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq2 {
        Fq2(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq2) {
        *self.0.y_mut() = y.0
    }

    pub fn from_jacobian(g2: G2) -> Option<Self> {
        g2.0.to_affine().map(AffineG2)
    }

    pub fn get_ys_from_x_unchecked(x: Fq2) -> Option<(Fq2, Fq2)> {
        groups::AffineG2::get_ys_from_x_unchecked(x.0).map(|(neq_y, y)| (Fq2(neq_y), Fq2(y)))
    }
}

impl Neg for AffineG2 {
    type Output = AffineG2;

    fn neg(self) -> AffineG2 {
        AffineG2(-self.0)
    }
}

impl From<AffineG2> for G2 {
    fn from(affine: AffineG2) -> Self {
        G2(affine.0.to_jacobian())
    }
}

impl From<G2> for AffineG2 {
    fn from(g2: G2) -> Self {
        AffineG2::new(g2.x() / g2.z(), g2.y() / g2.z()).expect("Unable to convert G2 to AffineG2")
    }
}

#[cfg(test)]
mod tests {
    use super::{Fq, Fq2, G1, G2};
    use alloc::vec::Vec;

    fn hex(s: &'static str) -> Vec<u8> {
        use rustc_hex::FromHex;
        s.from_hex().unwrap()
    }

    #[test]
    fn g1_from_compressed() {
        let g1 = G1::from_compressed(&hex(
            "0230644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46",
        ))
        .expect("Invalid g1 decompress result");
        assert_eq!(
            g1.x(),
            Fq::from_str(
                "21888242871839275222246405745257275088696311157297823662689037894645226208582"
            )
            .unwrap()
        );
        assert_eq!(
            g1.y(),
            Fq::from_str(
                "3969792565221544645472939191694882283483352126195956956354061729942568608776"
            )
            .unwrap()
        );
        assert_eq!(g1.z(), Fq::one());
    }

    #[test]
    fn g2_from_compressed() {
        let g2 = G2::from_compressed(
            &hex("0a023aed31b5a9e486366ea9988b05dba469c6206e58361d9c065bbea7d928204a761efc6e4fa08ed227650134b52c7f7dd0463963e8a4bf21f4899fe5da7f984a")
        ).expect("Valid g2 point hex encoding");

        assert_eq!(
            g2.x(),
            Fq2::new(
                Fq::from_str(
                    "5923585509243758863255447226263146374209884951848029582715967108651637186684"
                )
                .unwrap(),
                Fq::from_str(
                    "5336385337059958111259504403491065820971993066694750945459110579338490853570"
                )
                .unwrap(),
            )
        );

        assert_eq!(
            g2.y(),
            Fq2::new(
                Fq::from_str(
                    "10374495865873200088116930399159835104695426846400310764827677226300185211748"
                )
                .unwrap(),
                Fq::from_str(
                    "5256529835065685814318509161957442385362539991735248614869838648137856366932"
                )
                .unwrap(),
            )
        );

        // 0b prefix is point reflection on the curve
        let g2 = -G2::from_compressed(
            &hex("0b023aed31b5a9e486366ea9988b05dba469c6206e58361d9c065bbea7d928204a761efc6e4fa08ed227650134b52c7f7dd0463963e8a4bf21f4899fe5da7f984a")
        ).expect("Valid g2 point hex encoding");

        assert_eq!(
            g2.x(),
            Fq2::new(
                Fq::from_str(
                    "5923585509243758863255447226263146374209884951848029582715967108651637186684"
                )
                .unwrap(),
                Fq::from_str(
                    "5336385337059958111259504403491065820971993066694750945459110579338490853570"
                )
                .unwrap(),
            )
        );

        assert_eq!(
            g2.y(),
            Fq2::new(
                Fq::from_str(
                    "10374495865873200088116930399159835104695426846400310764827677226300185211748"
                )
                .unwrap(),
                Fq::from_str(
                    "5256529835065685814318509161957442385362539991735248614869838648137856366932"
                )
                .unwrap(),
            )
        );

        // valid point but invalid sign prefix
        assert!(
            G2::from_compressed(
                &hex("0c023aed31b5a9e486366ea9988b05dba469c6206e58361d9c065bbea7d928204a761efc6e4fa08ed227650134b52c7f7dd0463963e8a4bf21f4899fe5da7f984a")
            ).is_err()
        );
    }
}
