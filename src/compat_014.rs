use crate::{Bls12381G1, Bls12381G2, Scalar};
use elliptic_curve_014::{
    bigint::{Odd, U384},
    consts::U48,
    rand_core::TryRng,
    Curve, PrimeCurve,
};
use subtle::{Choice, ConstantTimeEq, CtOption};

impl ff_014::Field for Scalar {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    fn try_random<R: TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let mut buf = [0u8; 64];
        rng.try_fill_bytes(&mut buf)?;
        Ok(Self::from_bytes_wide(&buf))
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn double(&self) -> Self {
        self.double()
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff_014::helpers::sqrt_ratio_generic(num, div)
    }

    fn sqrt(&self) -> CtOption<Self> {
        ff_014::helpers::sqrt_tonelli_shanks(
            self,
            [
                0x7fff_2dff_7fff_ffff,
                0x04d0_ec02_a9de_d201,
                0x94ce_bea4_199c_ec04,
                0x0000_0000_39f6_d3a9,
            ],
        )
    }

    fn is_zero_vartime(&self) -> bool {
        self.0 == Self::ZERO.0
    }
}

impl ff_014::PrimeField for Scalar {
    type Repr = [u8; 32];

    fn from_repr(mut repr: Self::Repr) -> CtOption<Self> {
        Self::from_le_bytes(&repr).or_else(|| {
            repr.reverse();
            let mut bytes = [0u8; 48];
            bytes[16..].copy_from_slice(&repr);
            let scalar = Self::from_okm(&bytes);
            CtOption::new(scalar, !scalar.ct_eq(&Self::ZERO))
        })
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_le_bytes()
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_le_bytes()[0] & 1)
    }

    const MODULUS: &'static str =
        "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
    const NUM_BITS: u32 = <Scalar as ff::PrimeField>::NUM_BITS;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const TWO_INV: Self = <Scalar as ff::PrimeField>::TWO_INV;
    const MULTIPLICATIVE_GENERATOR: Self = <Scalar as ff::PrimeField>::MULTIPLICATIVE_GENERATOR;
    const S: u32 = <Scalar as ff::PrimeField>::S;
    const ROOT_OF_UNITY: Self = <Scalar as ff::PrimeField>::ROOT_OF_UNITY;
    const ROOT_OF_UNITY_INV: Self = <Scalar as ff::PrimeField>::ROOT_OF_UNITY_INV;
    const DELTA: Self = <Scalar as ff::PrimeField>::DELTA;
}

#[cfg(all(feature = "bits", not(target_pointer_width = "64")))]
type ReprBits014 = [u32; 8];

#[cfg(all(feature = "bits", target_pointer_width = "64"))]
type ReprBits014 = [u64; 4];

#[cfg(feature = "bits")]
impl ff_014::PrimeFieldBits for Scalar {
    type ReprBits = ReprBits014;

    fn to_le_bits(&self) -> ff_014::FieldBits<Self::ReprBits> {
        let bytes = self.to_le_bytes();

        #[cfg(not(target_pointer_width = "64"))]
        let limbs = [
            u32::from_le_bytes(bytes[0..4].try_into().unwrap()),
            u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
            u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
            u32::from_le_bytes(bytes[12..16].try_into().unwrap()),
            u32::from_le_bytes(bytes[16..20].try_into().unwrap()),
            u32::from_le_bytes(bytes[20..24].try_into().unwrap()),
            u32::from_le_bytes(bytes[24..28].try_into().unwrap()),
            u32::from_le_bytes(bytes[28..32].try_into().unwrap()),
        ];

        #[cfg(target_pointer_width = "64")]
        let limbs = [
            u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
            u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
            u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
            u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
        ];

        ff_014::FieldBits::new(limbs)
    }

    fn char_le_bits() -> ff_014::FieldBits<Self::ReprBits> {
        #[cfg(not(target_pointer_width = "64"))]
        {
            ff_014::FieldBits::new([
                0x0000_0001,
                0xffff_ffff,
                0xfffe_5bfe,
                0x53bd_a402,
                0x09a1_d805,
                0x3339_d808,
                0x299d_7d48,
                0x73ed_a753,
            ])
        }

        #[cfg(target_pointer_width = "64")]
        ff_014::FieldBits::new([
            0xffff_ffff_0000_0001,
            0x53bd_a402_fffe_5bfe,
            0x3339_d808_09a1_d805,
            0x73ed_a753_299d_7d48,
        ])
    }
}

impl Curve for Bls12381G1 {
    type FieldBytesSize = U48;
    type Uint = U384;

    const ORDER: Odd<U384> = Odd::<U384>::from_be_hex(
        "0000000000000000000000000000000073eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    );
}

impl PrimeCurve for Bls12381G1 {}

impl elliptic_curve_014::point::PointCompression for Bls12381G1 {
    const COMPRESS_POINTS: bool = true;
}

impl Curve for Bls12381G2 {
    type FieldBytesSize = U48;
    type Uint = U384;

    const ORDER: Odd<U384> = Odd::<U384>::from_be_hex(
        "0000000000000000000000000000000073eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    );
}

impl PrimeCurve for Bls12381G2 {}

impl elliptic_curve_014::point::PointCompression for Bls12381G2 {
    const COMPRESS_POINTS: bool = true;
}

#[cfg(feature = "groups")]
#[allow(deprecated)]
mod groups {
    use crate::{G1Affine, G1Compressed, G1Projective, G1Uncompressed, G2Affine, G2Compressed};
    use crate::{G2Projective, G2Uncompressed, Scalar};
    use elliptic_curve::consts::{U48, U96};
    use elliptic_curve::generic_array::GenericArray;
    use elliptic_curve::hash2curve::ExpandMsgXmd;
    use elliptic_curve::hash2curve::Sgn0;
    use elliptic_curve_014::rand_core::TryRng;
    use group_014::{Curve, CurveAffine, Group, GroupEncoding, UncompressedEncoding};
    use subtle::{Choice, CtOption};

    impl elliptic_curve_014::point::AffineCoordinates for G1Affine {
        type FieldRepr = GenericArray<u8, U48>;

        fn from_coordinates(x: &Self::FieldRepr, y: &Self::FieldRepr) -> CtOption<Self> {
            let mut bytes = [0u8; G1Affine::UNCOMPRESSED_BYTES];
            bytes[..48].copy_from_slice(x.as_ref());
            bytes[48..].copy_from_slice(y.as_ref());
            G1Affine::from_uncompressed(&bytes)
        }

        fn x(&self) -> Self::FieldRepr {
            GenericArray::<u8, U48>::clone_from_slice(&self.x.to_bytes())
        }

        fn y(&self) -> Self::FieldRepr {
            GenericArray::<u8, U48>::clone_from_slice(&self.y.to_bytes())
        }

        fn x_is_odd(&self) -> Choice {
            self.x.sgn0()
        }

        fn y_is_odd(&self) -> Choice {
            self.y.sgn0()
        }
    }

    impl Group for G1Projective {
        type Scalar = Scalar;

        fn try_random<R: TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
            let mut ikm = [0u8; 32];
            rng.try_fill_bytes(&mut ikm)?;
            Ok(Self::hash::<ExpandMsgXmd<sha2::Sha256>>(
                &ikm,
                b"BLS12381G1_XMD:SHA-256_SSWU_RO_",
            ))
        }

        fn identity() -> Self {
            Self::IDENTITY
        }

        fn generator() -> Self {
            Self::GENERATOR
        }

        fn is_identity(&self) -> Choice {
            self.is_identity()
        }

        fn double(&self) -> Self {
            self.double()
        }
    }

    #[cfg(feature = "alloc")]
    impl group_014::WnafGroup for G1Projective {
        fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
            const RECOMMENDATIONS: [usize; 12] =
                [1, 3, 7, 20, 43, 120, 273, 563, 1630, 3128, 7933, 62569];

            let mut ret = 4;
            for r in &RECOMMENDATIONS {
                if num_scalars > *r {
                    ret += 1;
                } else {
                    break;
                }
            }

            ret
        }
    }

    impl Curve for G1Projective {
        type Affine = G1Affine;

        fn batch_normalize(points: &[Self], out: &mut [Self::Affine]) {
            Self::batch_normalize(points, out);
        }

        fn to_affine(&self) -> Self::Affine {
            self.into()
        }
    }

    impl group_014::prime::PrimeGroup for G1Projective {}
    impl group_014::prime::PrimeCurve for G1Projective {}

    impl CurveAffine for G1Affine {
        type Curve = G1Projective;
        type Scalar = Scalar;

        fn identity() -> Self {
            Self::identity()
        }

        fn generator() -> Self {
            Self::generator()
        }

        fn is_identity(&self) -> Choice {
            self.is_identity()
        }

        fn to_curve(&self) -> Self::Curve {
            self.into()
        }
    }

    impl GroupEncoding for G1Projective {
        type Repr = G1Compressed;

        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            <G1Affine as GroupEncoding>::from_bytes(bytes).map(Self::from)
        }

        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            <G1Affine as GroupEncoding>::from_bytes_unchecked(bytes).map(Self::from)
        }

        fn to_bytes(&self) -> Self::Repr {
            <G1Affine as GroupEncoding>::to_bytes(&G1Affine::from(self))
        }
    }

    impl GroupEncoding for G1Affine {
        type Repr = G1Compressed;

        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            <Self as group::GroupEncoding>::from_bytes(bytes)
        }

        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            <Self as group::GroupEncoding>::from_bytes_unchecked(bytes)
        }

        fn to_bytes(&self) -> Self::Repr {
            <Self as group::GroupEncoding>::to_bytes(self)
        }
    }

    impl UncompressedEncoding for G1Affine {
        type Uncompressed = G1Uncompressed;

        fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
            <Self as group::UncompressedEncoding>::from_uncompressed(bytes)
        }

        fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
            <Self as group::UncompressedEncoding>::from_uncompressed_unchecked(bytes)
        }

        fn to_uncompressed(&self) -> Self::Uncompressed {
            <Self as group::UncompressedEncoding>::to_uncompressed(self)
        }
    }

    impl elliptic_curve_014::point::AffineCoordinates for G2Affine {
        type FieldRepr = GenericArray<u8, U96>;

        fn from_coordinates(x: &Self::FieldRepr, y: &Self::FieldRepr) -> CtOption<Self> {
            let mut bytes = [0u8; G2Affine::UNCOMPRESSED_BYTES];
            bytes[..96].copy_from_slice(x.as_ref());
            bytes[96..].copy_from_slice(y.as_ref());
            G2Affine::from_uncompressed(&bytes)
        }

        fn x(&self) -> Self::FieldRepr {
            let mut res = GenericArray::<u8, U96>::default();
            res[0..48].copy_from_slice(&self.x.c1.to_bytes()[..]);
            res[48..96].copy_from_slice(&self.x.c0.to_bytes()[..]);
            res
        }

        fn y(&self) -> Self::FieldRepr {
            let mut res = GenericArray::<u8, U96>::default();
            res[0..48].copy_from_slice(&self.y.c1.to_bytes()[..]);
            res[48..96].copy_from_slice(&self.y.c0.to_bytes()[..]);
            res
        }

        fn x_is_odd(&self) -> Choice {
            self.x.sgn0()
        }

        fn y_is_odd(&self) -> Choice {
            self.y.sgn0()
        }
    }

    impl Group for G2Projective {
        type Scalar = Scalar;

        fn try_random<R: TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
            let mut ikm = [0u8; 32];
            rng.try_fill_bytes(&mut ikm)?;
            Ok(Self::hash::<ExpandMsgXmd<sha2::Sha256>>(
                &ikm,
                b"BLS12381G2_XMD:SHA-256_SSWU_RO_",
            ))
        }

        fn identity() -> Self {
            Self::IDENTITY
        }

        fn generator() -> Self {
            Self::GENERATOR
        }

        fn is_identity(&self) -> Choice {
            self.is_identity()
        }

        fn double(&self) -> Self {
            self.double()
        }
    }

    #[cfg(feature = "alloc")]
    impl group_014::WnafGroup for G2Projective {
        fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
            const RECOMMENDATIONS: [usize; 11] =
                [1, 3, 8, 20, 47, 126, 260, 826, 1501, 4555, 84071];

            let mut ret = 4;
            for r in &RECOMMENDATIONS {
                if num_scalars > *r {
                    ret += 1;
                } else {
                    break;
                }
            }

            ret
        }
    }

    impl Curve for G2Projective {
        type Affine = G2Affine;

        fn batch_normalize(points: &[Self], out: &mut [Self::Affine]) {
            Self::batch_normalize(points, out);
        }

        fn to_affine(&self) -> Self::Affine {
            self.into()
        }
    }

    impl group_014::prime::PrimeGroup for G2Projective {}
    impl group_014::prime::PrimeCurve for G2Projective {}

    impl CurveAffine for G2Affine {
        type Curve = G2Projective;
        type Scalar = Scalar;

        fn identity() -> Self {
            Self::identity()
        }

        fn generator() -> Self {
            Self::generator()
        }

        fn is_identity(&self) -> Choice {
            self.is_identity()
        }

        fn to_curve(&self) -> Self::Curve {
            self.into()
        }
    }

    impl GroupEncoding for G2Projective {
        type Repr = G2Compressed;

        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            <G2Affine as GroupEncoding>::from_bytes(bytes).map(Self::from)
        }

        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            <G2Affine as GroupEncoding>::from_bytes_unchecked(bytes).map(Self::from)
        }

        fn to_bytes(&self) -> Self::Repr {
            <G2Affine as GroupEncoding>::to_bytes(&G2Affine::from(self))
        }
    }

    impl GroupEncoding for G2Affine {
        type Repr = G2Compressed;

        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            <Self as group::GroupEncoding>::from_bytes(bytes)
        }

        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            <Self as group::GroupEncoding>::from_bytes_unchecked(bytes)
        }

        fn to_bytes(&self) -> Self::Repr {
            <Self as group::GroupEncoding>::to_bytes(self)
        }
    }

    impl UncompressedEncoding for G2Affine {
        type Uncompressed = G2Uncompressed;

        fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
            <Self as group::UncompressedEncoding>::from_uncompressed(bytes)
        }

        fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
            <Self as group::UncompressedEncoding>::from_uncompressed_unchecked(bytes)
        }

        fn to_uncompressed(&self) -> Self::Uncompressed {
            <Self as group::UncompressedEncoding>::to_uncompressed(self)
        }
    }

    #[cfg(feature = "pairings")]
    mod gt {
        use crate::{pairings::GtRepr, Gt, Scalar};
        use ff_014::Field;
        use group_014::{Group, GroupEncoding};
        use subtle::{Choice, ConstantTimeEq, CtOption};

        impl Group for Gt {
            type Scalar = Scalar;

            fn try_random<R: elliptic_curve_014::rand_core::TryRng + ?Sized>(
                rng: &mut R,
            ) -> Result<Self, R::Error> {
                loop {
                    let scalar = <Scalar as Field>::try_random(rng)?;
                    if !bool::from(scalar.is_zero()) {
                        return Ok(Self::generator() * scalar);
                    }
                }
            }

            fn identity() -> Self {
                Self::IDENTITY
            }

            fn generator() -> Self {
                <Self as group::Group>::generator()
            }

            fn is_identity(&self) -> Choice {
                self.ct_eq(&Self::IDENTITY)
            }

            fn double(&self) -> Self {
                self.double()
            }
        }

        impl GroupEncoding for Gt {
            type Repr = GtRepr;

            fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
                <Self as group::GroupEncoding>::from_bytes(bytes)
            }

            fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
                <Self as group::GroupEncoding>::from_bytes_unchecked(bytes)
            }

            fn to_bytes(&self) -> Self::Repr {
                <Self as group::GroupEncoding>::to_bytes(self)
            }
        }
    }
}
