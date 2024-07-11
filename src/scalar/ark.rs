use ark_bls12_381::Fr;
use ark_ff::{
    BigInteger, BigInteger256, FftField, Field as ArkField, LegendreSymbol, One as ArkOne,
    PrimeField as ArkPrimeField, SqrtPrecomputation, UniformRand, Zero as ArkZero,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Compress, Flags, Read, SerializationError, Valid, Validate, Write,
};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
use core::str::FromStr;
use ff::{Field, PrimeField};
use num_bigint::BigUint;
use num_traits::Num;
use core::iter;
use subtle::ConstantTimeEq;

use crate::Scalar;

impl<'a> Add<&'a mut Scalar> for Scalar {
    type Output = Scalar;

    fn add(self, rhs: &'a mut Scalar) -> Self::Output {
        self + *rhs
    }
}

impl<'a> Sub<&'a mut Scalar> for Scalar {
    type Output = Scalar;

    fn sub(self, rhs: &'a mut Scalar) -> Self::Output {
        self - *rhs
    }
}

impl<'a> Mul<&'a mut Scalar> for Scalar {
    type Output = Scalar;

    fn mul(self, rhs: &'a mut Scalar) -> Self::Output {
        self * *rhs
    }
}

impl<'a> Div<&'a mut Scalar> for Scalar {
    type Output = Scalar;

    fn div(self, rhs: &'a mut Scalar) -> Self::Output {
        self / *rhs
    }
}

impl<'a> AddAssign<&'a mut Scalar> for Scalar {
    fn add_assign(&mut self, rhs: &'a mut Scalar) {
        *self += *rhs;
    }
}

impl<'a> SubAssign<&'a mut Scalar> for Scalar {
    fn sub_assign(&mut self, rhs: &'a mut Scalar) {
        *self -= *rhs;
    }
}

impl<'a> MulAssign<&'a mut Scalar> for Scalar {
    fn mul_assign(&mut self, rhs: &'a mut Scalar) {
        *self *= *rhs;
    }
}

impl<'a> DivAssign<&'a mut Scalar> for Scalar {
    fn div_assign(&mut self, rhs: &'a mut Scalar) {
        *self /= *rhs;
    }
}

unsafe impl Send for Scalar {}

unsafe impl Sync for Scalar {}

impl CanonicalSerialize for Scalar {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress,
    ) -> Result<(), SerializationError> {
        let bytes = self.to_le_bytes();
        Ok(writer.write_all(&bytes)?)
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        32
    }
}

impl CanonicalSerializeWithFlags for Scalar {
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: W,
        flags: F,
    ) -> Result<(), SerializationError> {
        let bytes = self.to_le_bytes();
        let fr = ark_bls12_381::Fr::deserialize_compressed(&bytes[..])
            .expect("Failed to serialize scalar");
        fr.serialize_with_flags(writer, flags)
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        let bytes = self.to_le_bytes();
        let fr = ark_bls12_381::Fr::deserialize_compressed(&bytes[..])
            .expect("Failed to serialize scalar");
        <ark_bls12_381::Fr as CanonicalSerializeWithFlags>::serialized_size_with_flags::<F>(&fr)
    }
}

impl CanonicalDeserialize for Scalar {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; 32];
        reader.read_exact(&mut bytes)?;
        let s = Self::from_le_bytes(&bytes);
        match validate {
            Validate::No => Ok(s.expect("a valid scalar")),
            Validate::Yes => Option::<Self>::from(s).ok_or(SerializationError::InvalidData),
        }
    }
}

impl CanonicalDeserializeWithFlags for Scalar {
    fn deserialize_with_flags<R: Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        let mut bytes = [0u8; 32];
        reader.read_exact(&mut bytes)?;
        let s = Self::from_le_bytes(&bytes);
        Option::<Self>::from(s)
            .map(|s| (s, F::default()))
            .ok_or(SerializationError::InvalidData)
    }
}

impl Valid for Scalar {
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl ArkOne for Scalar {
    fn one() -> Self {
        Self::ONE
    }

    fn set_one(&mut self) {
        *self = Self::ONE;
    }

    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        *self == Self::ONE
    }
}

impl ArkZero for Scalar {
    fn zero() -> Self {
        Self::ZERO
    }

    fn set_zero(&mut self) {
        *self = Self::ZERO;
    }

    fn is_zero(&self) -> bool
    where
        Self: PartialEq,
    {
        *self == Self::ZERO
    }
}

impl UniformRand for Scalar {
    fn rand<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
        Self::random(rng)
    }
}

impl FftField for Scalar {
    const GENERATOR: Self = <Self as PrimeField>::MULTIPLICATIVE_GENERATOR;
    const TWO_ADICITY: u32 = <Self as PrimeField>::S;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self::from_raw_unchecked([
        0x3829971f439f0d2bu64,
        0xb63683508c2280b9u64,
        0xd09b681922c813b4u64,
        0x16a2a19edfe81f20u64,
    ]);
    const SMALL_SUBGROUP_BASE: Option<u32> = Some(3);
    const SMALL_SUBGROUP_BASE_ADICITY: Option<u32> = Some(1);
    const LARGE_SUBGROUP_ROOT_OF_UNITY: Option<Self> = Some(Self::from_raw_unchecked([
        0x02b93785357e7917,
        0x85aedb297ca15150,
        0xea45ce5f9f533109,
        0x3c1d00c4965f33c8,
    ]));
}

impl FromStr for Scalar {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Ark expects base 10
        let num = BigUint::from_str_radix(s, 10).map_err(|_| ())?;
        let modulus = BigUint::from_bytes_be(&<Self as ArkPrimeField>::MODULUS.to_bytes_be());
        if num > modulus {
            return Err(());
        }
        Ok(Self::from(num))
    }
}

impl From<BigInteger256> for Scalar {
    fn from(repr: BigInteger256) -> Self {
        Self::from_bigint(repr).expect("Failed to convert from BigInt")
    }
}

impl From<Scalar> for BigInteger256 {
    fn from(scalar: Scalar) -> Self {
        scalar.into_bigint()
    }
}

impl From<BigUint> for Scalar {
    fn from(value: BigUint) -> Self {
        let modulus = BigUint::from_bytes_be(&<Self as ArkPrimeField>::MODULUS.to_bytes_be());
        let value = value % modulus;
        let bytes = value.to_bytes_be();
        let mut be_bytes = [0u8; 32];
        be_bytes[32 - bytes.len()..].copy_from_slice(&bytes);
        Self::from_be_bytes(&be_bytes).expect("Failed to convert from BigUint")
    }
}

impl From<Scalar> for BigUint {
    fn from(scalar: Scalar) -> Self {
        BigUint::from_bytes_be(&scalar.to_be_bytes())
    }
}

impl ArkPrimeField for Scalar {
    type BigInt = BigInteger256;
    const MODULUS: Self::BigInt = BigInteger256::new([
        0xffff_ffff_0000_0001,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ]);
    const MODULUS_MINUS_ONE_DIV_TWO: Self::BigInt = BigInteger256::new([
        0x7fff_2dff_7fff_ffff,
        0x04d0_ec02_a9de_d201,
        0x94ce_bea4_199c_ec04,
        0x0000_0000_39f6_d3a9,
    ]);
    const MODULUS_BIT_SIZE: u32 = super::MODULUS_BITS;
    const TRACE: Self::BigInt = BigInteger256::new([
        0xfffe5bfeffffffff,
        0x09a1d80553bda402,
        0x299d7d483339d808,
        0x0000000073eda753,
    ]);
    const TRACE_MINUS_ONE_DIV_TWO: Self::BigInt = BigInteger256::new([
        0x7fff2dff7fffffff,
        0x04d0ec02a9ded201,
        0x94cebea4199cec04,
        0x0000000039f6d3a9,
    ]);

    fn from_bigint(repr: Self::BigInt) -> Option<Self> {
        let modulus = <Self as ArkPrimeField>::MODULUS;
        if repr >= modulus {
            return None;
        }
        Option::<Self>::from(Self::from_raw(repr.0))
    }

    fn into_bigint(self) -> Self::BigInt {
        let tmp = Self::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);
        BigInteger256::new(tmp.0)
    }
}

impl ArkField for Scalar {
    type BasePrimeField = Self;
    type BasePrimeFieldIter = iter::Once<Self::BasePrimeField>;
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = Some(SqrtPrecomputation::Case3Mod4 {
        modulus_plus_one_div_four: &[
            0xbfffffffc0000000,
            0x54ef6900bfff96ff,
            0x0cce760202687601,
            0x1cfb69d4ca675f52,
        ],
    });
    const ZERO: Self = <Self as Field>::ZERO;
    const ONE: Self = <Self as Field>::ONE;

    fn extension_degree() -> u64 {
        1
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        iter::once(*self)
    }

    fn from_base_prime_field_elems(elems: &[Self::BasePrimeField]) -> Option<Self> {
        if elems.len() != (Self::extension_degree() as usize) {
            return None;
        }
        Some(elems[0])
    }

    fn from_base_prime_field(elem: Self::BasePrimeField) -> Self {
        elem
    }

    fn double(&self) -> Self {
        self.double()
    }

    fn double_in_place(&mut self) -> &mut Self {
        *self = self.double();
        self
    }

    fn neg_in_place(&mut self) -> &mut Self {
        *self = -*self;
        self
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        ark_bls12_381::Fr::from_random_bytes_with_flags(bytes)
            .map(|(fr, flags)| (Scalar::from(fr.0), flags))
    }

    fn legendre(&self) -> LegendreSymbol {
        // s = self^((MODULUS - 1) // 2)
        let s = self.pow(&[
            0x7fff_2dff_7fff_ffff,
            0x04d0_ec02_a9de_d201,
            0x94ce_bea4_199c_ec04,
            0x0000_0000_39f6_d3a9,
        ]);
        if s.ct_eq(&Self::ZERO).into() {
            LegendreSymbol::Zero
        } else if s.ct_eq(&Self::ONE).into() {
            LegendreSymbol::QuadraticResidue
        } else {
            LegendreSymbol::QuadraticNonResidue
        }
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn square_in_place(&mut self) -> &mut Self {
        *self = self.square();
        self
    }

    fn inverse(&self) -> Option<Self> {
        Option::<Self>::from(self.invert())
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        let inv = Option::<Self>::from(self.invert());
        if let Some(inv) = inv {
            *self = inv;
            Some(self)
        } else {
            None
        }
    }

    fn frobenius_map_in_place(&mut self, _: usize) {}
}

impl From<Fr> for Scalar {
    fn from(value: Fr) -> Self {
        let mut bytes = std::vec::Vec::with_capacity(32);
        value
            .0
            .serialize_compressed(&mut bytes)
            .expect("Failed to serialize Fr");
        debug_assert_eq!(bytes.len(), 32);
        Self::from_le_bytes(&<[u8; 32]>::try_from(bytes.as_slice()).unwrap())
            .expect("Failed to convert bytes to Scalar")
    }
}

impl From<Scalar> for Fr {
    fn from(value: Scalar) -> Self {
        let bytes = value.to_le_bytes();
        Fr::deserialize_compressed(&bytes[..]).expect("Failed to deserialize Fr")
    }
}