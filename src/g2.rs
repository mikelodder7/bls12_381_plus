//! This module provides an implementation of the $\mathbb{G}_2$ group of BLS12-381.

use core::borrow::Borrow;
use core::fmt;
use core::iter::Sum;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use group::{
    prime::{PrimeCurve, PrimeCurveAffine, PrimeGroup},
    Curve, Group, GroupEncoding, UncompressedEncoding,
};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "alloc")]
use group::WnafGroup;

use crate::fp::Fp;
use crate::fp2::Fp2;
#[cfg(feature = "hashing")]
use crate::hash_to_field::ExpandMsg;
use crate::Scalar;

/// This is an element of $\mathbb{G}_2$ represented in the affine coordinate space.
/// It is ideal to keep elements in this representation to reduce memory usage and
/// improve performance through the use of mixed curve model arithmetic.
///
/// Values of `G2Affine` are guaranteed to be in the $q$-order subgroup unless an
/// "unchecked" API was misused.
#[cfg_attr(docsrs, doc(cfg(feature = "groups")))]
#[derive(Copy, Clone, Debug)]
pub struct G2Affine {
    pub(crate) x: Fp2,
    pub(crate) y: Fp2,
    infinity: Choice,
}

impl Default for G2Affine {
    fn default() -> G2Affine {
        G2Affine::identity()
    }
}

impl fmt::Display for G2Affine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a G2Projective> for G2Affine {
    fn from(p: &'a G2Projective) -> G2Affine {
        let zinv = p.z.invert().unwrap_or(Fp2::zero());
        let x = p.x * zinv;
        let y = p.y * zinv;

        let tmp = G2Affine {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        G2Affine::conditional_select(&tmp, &G2Affine::identity(), zinv.is_zero())
    }
}

impl From<G2Projective> for G2Affine {
    fn from(p: G2Projective) -> G2Affine {
        G2Affine::from(&p)
    }
}

impl ConstantTimeEq for G2Affine {
    fn ct_eq(&self, other: &Self) -> Choice {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        (self.infinity & other.infinity)
            | ((!self.infinity)
                & (!other.infinity)
                & self.x.ct_eq(&other.x)
                & self.y.ct_eq(&other.y))
    }
}

impl ConditionallySelectable for G2Affine {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        G2Affine {
            x: Fp2::conditional_select(&a.x, &b.x, choice),
            y: Fp2::conditional_select(&a.y, &b.y, choice),
            infinity: Choice::conditional_select(&a.infinity, &b.infinity, choice),
        }
    }
}

impl Eq for G2Affine {}
impl PartialEq for G2Affine {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a G2Affine {
    type Output = G2Affine;

    #[inline]
    fn neg(self) -> G2Affine {
        G2Affine {
            x: self.x,
            y: Fp2::conditional_select(&-self.y, &Fp2::one(), self.infinity),
            infinity: self.infinity,
        }
    }
}

impl Neg for G2Affine {
    type Output = G2Affine;

    #[inline]
    fn neg(self) -> G2Affine {
        -&self
    }
}

impl<'a, 'b> Add<&'b G2Projective> for &'a G2Affine {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Projective) -> G2Projective {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b G2Affine> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Affine) -> G2Projective {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Projective> for &'a G2Affine {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Projective) -> G2Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Affine> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Affine) -> G2Projective {
        self + (-rhs)
    }
}

impl<T> Sum<T> for G2Projective
where
    T: Borrow<G2Projective>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl_binops_additive!(G2Projective, G2Affine);
impl_binops_additive_specify_output!(G2Affine, G2Projective, G2Projective);

const B: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0xaa27_0000_000c_fff3,
        0x53cc_0032_fc34_000a,
        0x478f_e97a_6b0a_807f,
        0xb1d3_7ebe_e6ba_24d7,
        0x8ec9_733b_bf78_ab2f,
        0x09d6_4551_3d83_de7e,
    ]),
    c1: Fp::from_raw_unchecked([
        0xaa27_0000_000c_fff3,
        0x53cc_0032_fc34_000a,
        0x478f_e97a_6b0a_807f,
        0xb1d3_7ebe_e6ba_24d7,
        0x8ec9_733b_bf78_ab2f,
        0x09d6_4551_3d83_de7e,
    ]),
};

const B3: Fp2 = Fp2::add(&Fp2::add(&B, &B), &B);

impl G2Affine {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> G2Affine {
        G2Affine {
            x: Fp2::zero(),
            y: Fp2::one(),
            infinity: Choice::from(1u8),
        }
    }

    /// Returns a fixed generator of the group. See [`notes::design`](notes/design/index.html#fixed-generators)
    /// for how this generator is chosen.
    pub fn generator() -> G2Affine {
        G2Affine {
            x: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xf5f2_8fa2_0294_0a10,
                    0xb3f5_fb26_87b4_961a,
                    0xa1a8_93b5_3e2a_e580,
                    0x9894_999d_1a3c_aee9,
                    0x6f67_b763_1863_366b,
                    0x0581_9192_4350_bcd7,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xa5a9_c075_9e23_f606,
                    0xaaa0_c59d_bccd_60c3,
                    0x3bb1_7e18_e286_7806,
                    0x1b1a_b6cc_8541_b367,
                    0xc2b6_ed0e_f215_8547,
                    0x1192_2a09_7360_edf3,
                ]),
            },
            y: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x4c73_0af8_6049_4c4a,
                    0x597c_fa1f_5e36_9c5a,
                    0xe7e6_856c_aa0a_635a,
                    0xbbef_b5e9_6e0d_495f,
                    0x07d3_a975_f0ef_25a2,
                    0x0083_fd8e_7e80_dae5,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xadc0_fc92_df64_b05d,
                    0x18aa_270a_2b14_61dc,
                    0x86ad_ac6a_3be4_eba0,
                    0x7949_5c4e_c93d_a33a,
                    0xe717_5850_a43c_caed,
                    0x0b2b_c2a1_63de_1bf2,
                ]),
            },
            infinity: Choice::from(0u8),
        }
    }

    /// Serializes this element into compressed form. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn to_compressed(&self) -> [u8; 96] {
        // Strictly speaking, self.x is zero already when self.infinity is true, but
        // to guard against implementation mistakes we do not assume this.
        let x = Fp2::conditional_select(&self.x, &Fp2::zero(), self.infinity);

        let mut res = [0; 96];

        (&mut res[0..48]).copy_from_slice(&x.c1.to_bytes()[..]);
        (&mut res[48..96]).copy_from_slice(&x.c0.to_bytes()[..]);

        // This point is in compressed form, so we set the most significant bit.
        res[0] |= 1u8 << 7;

        // Is this point at infinity? If so, set the second-most significant bit.
        res[0] |= u8::conditional_select(&0u8, &(1u8 << 6), self.infinity);

        // Is the y-coordinate the lexicographically largest of the two associated with the
        // x-coordinate? If so, set the third-most significant bit so long as this is not
        // the point at infinity.
        res[0] |= u8::conditional_select(
            &0u8,
            &(1u8 << 5),
            (!self.infinity) & self.y.lexicographically_largest(),
        );

        res
    }

    /// Serializes this element into uncompressed form. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn to_uncompressed(&self) -> [u8; 192] {
        let mut res = [0; 192];

        let x = Fp2::conditional_select(&self.x, &Fp2::zero(), self.infinity);
        let y = Fp2::conditional_select(&self.y, &Fp2::zero(), self.infinity);

        res[0..48].copy_from_slice(&x.c1.to_bytes()[..]);
        res[48..96].copy_from_slice(&x.c0.to_bytes()[..]);
        res[96..144].copy_from_slice(&y.c1.to_bytes()[..]);
        res[144..192].copy_from_slice(&y.c0.to_bytes()[..]);

        // Is this point at infinity? If so, set the second-most significant bit.
        res[0] |= u8::conditional_select(&0u8, &(1u8 << 6), self.infinity);

        res
    }

    /// Attempts to deserialize an uncompressed element. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn from_uncompressed(bytes: &[u8; 192]) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(bytes)
            .and_then(|p| CtOption::new(p, p.is_on_curve() & p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 192]) -> CtOption<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[0] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[0] >> 6) & 1);
        let sort_flag_set = Choice::from((bytes[0] >> 5) & 1);

        // Attempt to obtain the x-coordinate
        let xc1 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            Fp::from_bytes(&tmp)
        };
        let xc0 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

            Fp::from_bytes(&tmp)
        };

        // Attempt to obtain the y-coordinate
        let yc1 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[96..144]);

            Fp::from_bytes(&tmp)
        };
        let yc0 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[144..192]);

            Fp::from_bytes(&tmp)
        };

        xc1.and_then(|xc1| {
            xc0.and_then(|xc0| {
                yc1.and_then(|yc1| {
                    yc0.and_then(|yc0| {
                        let x = Fp2 {
                            c0: xc0,
                            c1: xc1
                        };
                        let y = Fp2 {
                            c0: yc0,
                            c1: yc1
                        };

                        // Create a point representing this value
                        let p = G2Affine::conditional_select(
                            &G2Affine {
                                x,
                                y,
                                infinity: infinity_flag_set,
                            },
                            &G2Affine::identity(),
                            infinity_flag_set,
                        );

                        CtOption::new(
                            p,
                            // If the infinity flag is set, the x and y coordinates should have been zero.
                            ((!infinity_flag_set) | (infinity_flag_set & x.is_zero() & y.is_zero())) &
                            // The compression flag should not have been set, as this is an uncompressed element
                            (!compression_flag_set) &
                            // The sort flag should not have been set, as this is an uncompressed element
                            (!sort_flag_set),
                        )
                    })
                })
            })
        })
    }

    /// Attempts to deserialize a compressed element. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn from_compressed(bytes: &[u8; 96]) -> CtOption<Self> {
        // We already know the point is on the curve because this is established
        // by the y-coordinate recovery procedure in from_compressed_unchecked().

        Self::from_compressed_unchecked(bytes).and_then(|p| CtOption::new(p, p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 96]) -> CtOption<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[0] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[0] >> 6) & 1);
        let sort_flag_set = Choice::from((bytes[0] >> 5) & 1);

        // Attempt to obtain the x-coordinate
        let xc1 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            Fp::from_bytes(&tmp)
        };
        let xc0 = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

            Fp::from_bytes(&tmp)
        };

        xc1.and_then(|xc1| {
            xc0.and_then(|xc0| {
                let x = Fp2 { c0: xc0, c1: xc1 };

                // If the infinity flag is set, return the value assuming
                // the x-coordinate is zero and the sort bit is not set.
                //
                // Otherwise, return a recovered point (assuming the correct
                // y-coordinate can be found) so long as the infinity flag
                // was not set.
                CtOption::new(
                    G2Affine::identity(),
                    infinity_flag_set & // Infinity flag should be set
                    compression_flag_set & // Compression flag should be set
                    (!sort_flag_set) & // Sort flag should not be set
                    x.is_zero(), // The x-coordinate should be zero
                )
                .or_else(|| {
                    // Recover a y-coordinate given x by y = sqrt(x^3 + 4)
                    ((x.square() * x) + B).sqrt().and_then(|y| {
                        // Switch to the correct y-coordinate if necessary.
                        let y = Fp2::conditional_select(
                            &y,
                            &-y,
                            y.lexicographically_largest() ^ sort_flag_set,
                        );

                        CtOption::new(
                            G2Affine {
                                x,
                                y,
                                infinity: infinity_flag_set,
                            },
                            (!infinity_flag_set) & // Infinity flag should not be set
                            compression_flag_set, // Compression flag should be set
                        )
                    })
                })
            })
        })
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.infinity
    }

    /// Returns true if this point is free of an $h$-torsion component, and so it
    /// exists within the $q$-order subgroup $\mathbb{G}_2$. This should always return true
    /// unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> Choice {
        const FQ_MODULUS_BYTES: [u8; 32] = [
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115,
        ];

        // Clear the r-torsion from the point and check if it is the identity
        G2Projective::from(*self)
            .multiply(&FQ_MODULUS_BYTES)
            .is_identity()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // y^2 - x^3 ?= 4(u + 1)
        (self.y.square() - (self.x.square() * self.x)).ct_eq(&B) | self.infinity
    }
}

impl_serde!(
    G2Affine,
    |p: &G2Affine| p.to_compressed(),
    |arr: &[u8; 96]| G2Affine::from_compressed(arr),
    96
);

/// This is an element of $\mathbb{G}_2$ represented in the projective coordinate space.
#[cfg_attr(docsrs, doc(cfg(feature = "groups")))]
#[derive(Copy, Clone, Debug)]
pub struct G2Projective {
    pub(crate) x: Fp2,
    pub(crate) y: Fp2,
    pub(crate) z: Fp2,
}

impl Default for G2Projective {
    fn default() -> G2Projective {
        G2Projective::identity()
    }
}

impl fmt::Display for G2Projective {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a G2Affine> for G2Projective {
    fn from(p: &'a G2Affine) -> G2Projective {
        G2Projective {
            x: p.x,
            y: p.y,
            z: Fp2::conditional_select(&Fp2::one(), &Fp2::zero(), p.infinity),
        }
    }
}

impl From<G2Affine> for G2Projective {
    fn from(p: G2Affine) -> G2Projective {
        G2Projective::from(&p)
    }
}

impl ConstantTimeEq for G2Projective {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz, yz, z) equal to (x'z', y'z', z') when converted to affine?

        let x1 = self.x * other.z;
        let x2 = other.x * self.z;

        let y1 = self.y * other.z;
        let y2 = other.y * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for G2Projective {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        G2Projective {
            x: Fp2::conditional_select(&a.x, &b.x, choice),
            y: Fp2::conditional_select(&a.y, &b.y, choice),
            z: Fp2::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl Eq for G2Projective {}
impl PartialEq for G2Projective {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn neg(self) -> G2Projective {
        G2Projective {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}

impl Neg for G2Projective {
    type Output = G2Projective;

    #[inline]
    fn neg(self) -> G2Projective {
        -&self
    }
}

impl<'a, 'b> Add<&'b G2Projective> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Projective) -> G2Projective {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Projective> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Projective) -> G2Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G2Projective {
    type Output = G2Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other.to_bytes())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G2Affine {
    type Output = G2Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        G2Projective::from(self).multiply(&other.to_bytes())
    }
}

impl_binops_additive!(G2Projective, G2Projective);
impl_binops_multiplicative!(G2Projective, Scalar);
impl_binops_multiplicative_mixed!(G2Affine, Scalar, G2Projective);

#[inline(always)]
fn mul_by_3b(x: Fp2) -> Fp2 {
    x * B3
}

impl G2Projective {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> G2Projective {
        G2Projective {
            x: Fp2::zero(),
            y: Fp2::one(),
            z: Fp2::zero(),
        }
    }

    /// Returns a fixed generator of the group. See [`notes::design`](notes/design/index.html#fixed-generators)
    /// for how this generator is chosen.
    pub fn generator() -> G2Projective {
        G2Projective {
            x: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xf5f2_8fa2_0294_0a10,
                    0xb3f5_fb26_87b4_961a,
                    0xa1a8_93b5_3e2a_e580,
                    0x9894_999d_1a3c_aee9,
                    0x6f67_b763_1863_366b,
                    0x0581_9192_4350_bcd7,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xa5a9_c075_9e23_f606,
                    0xaaa0_c59d_bccd_60c3,
                    0x3bb1_7e18_e286_7806,
                    0x1b1a_b6cc_8541_b367,
                    0xc2b6_ed0e_f215_8547,
                    0x1192_2a09_7360_edf3,
                ]),
            },
            y: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x4c73_0af8_6049_4c4a,
                    0x597c_fa1f_5e36_9c5a,
                    0xe7e6_856c_aa0a_635a,
                    0xbbef_b5e9_6e0d_495f,
                    0x07d3_a975_f0ef_25a2,
                    0x0083_fd8e_7e80_dae5,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xadc0_fc92_df64_b05d,
                    0x18aa_270a_2b14_61dc,
                    0x86ad_ac6a_3be4_eba0,
                    0x7949_5c4e_c93d_a33a,
                    0xe717_5850_a43c_caed,
                    0x0b2b_c2a1_63de_1bf2,
                ]),
            },
            z: Fp2::one(),
        }
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> G2Projective {
        // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.y.square();
        let z3 = t0 + t0;
        let z3 = z3 + z3;
        let z3 = z3 + z3;
        let t1 = self.y * self.z;
        let t2 = self.z.square();
        let t2 = mul_by_3b(t2);
        let x3 = t2 * z3;
        let y3 = t0 + t2;
        let z3 = t1 * z3;
        let t1 = t2 + t2;
        let t2 = t1 + t2;
        let t0 = t0 - t2;
        let y3 = t0 * y3;
        let y3 = x3 + y3;
        let t1 = self.x * self.y;
        let x3 = t0 * t1;
        let x3 = x3 + x3;

        let tmp = G2Projective {
            x: x3,
            y: y3,
            z: z3,
        };

        G2Projective::conditional_select(&tmp, &G2Projective::identity(), self.is_identity())
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &G2Projective) -> G2Projective {
        // Algorithm 7, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t2 = self.z * rhs.z;
        let t3 = self.x + self.y;
        let t4 = rhs.x + rhs.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = self.y + self.z;
        let x3 = rhs.y + rhs.z;
        let t4 = t4 * x3;
        let x3 = t1 + t2;
        let t4 = t4 - x3;
        let x3 = self.x + self.z;
        let y3 = rhs.x + rhs.z;
        let x3 = x3 * y3;
        let y3 = t0 + t2;
        let y3 = x3 - y3;
        let x3 = t0 + t0;
        let t0 = x3 + t0;
        let t2 = mul_by_3b(t2);
        let z3 = t1 + t2;
        let t1 = t1 - t2;
        let y3 = mul_by_3b(y3);
        let x3 = t4 * y3;
        let t2 = t3 * t1;
        let x3 = t2 - x3;
        let y3 = y3 * t0;
        let t1 = t1 * z3;
        let y3 = t1 + y3;
        let t0 = t0 * t3;
        let z3 = z3 * t4;
        let z3 = z3 + t0;

        G2Projective {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &G2Affine) -> G2Projective {
        // Algorithm 8, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t3 = rhs.x + rhs.y;
        let t4 = self.x + self.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = rhs.y * self.z;
        let t4 = t4 + self.y;
        let y3 = rhs.x * self.z;
        let y3 = y3 + self.x;
        let x3 = t0 + t0;
        let t0 = x3 + t0;
        let t2 = mul_by_3b(self.z);
        let z3 = t1 + t2;
        let t1 = t1 - t2;
        let y3 = mul_by_3b(y3);
        let x3 = t4 * y3;
        let t2 = t3 * t1;
        let x3 = t2 - x3;
        let y3 = y3 * t0;
        let t1 = t1 * z3;
        let y3 = t1 + y3;
        let t0 = t0 * t3;
        let z3 = z3 * t4;
        let z3 = z3 + t0;

        let tmp = G2Projective {
            x: x3,
            y: y3,
            z: z3,
        };

        G2Projective::conditional_select(&tmp, &self, rhs.is_identity())
    }

    fn multiply(&self, by: &[u8]) -> G2Projective {
        let mut acc = G2Projective::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading bit because it's always unset for Fq
        // elements.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G2Projective::conditional_select(&acc, &(acc + self), bit);
        }

        acc
    }

    #[cfg(feature = "endo")]
    fn psi(&self) -> G2Projective {
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2 {
            c0: Fp::zero(),
            c1: Fp::from_raw_unchecked([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
                0x14e4f04fe2db9068,
                0x14e56d3f1564853a,
            ]),
        };
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2 {
            c0: Fp::from_raw_unchecked([
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
                0x1d794e4fac7cf0b9,
                0x0bd592fc7d825ec8,
            ]),
            c1: Fp::from_raw_unchecked([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0x0e2b7eedbbfd87d2,
            ]),
        };

        G2Projective {
            // x = frobenius(x)/((u+1)^((p-1)/3))
            x: self.x.frobenius_map() * psi_coeff_x,
            // y = frobenius(y)/(u+1)^((p-1)/2)
            y: self.y.frobenius_map() * psi_coeff_y,
            // z = frobenius(z)
            z: self.z.frobenius_map(),
        }
    }

    #[cfg(feature = "endo")]
    fn psi2(&self) -> G2Projective {
        // 1 / 2 ^ ((q-1)/3)
        let psi2_coeff_x = Fp2 {
            c0: Fp::from_raw_unchecked([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x03f97d6e83d050d2,
                0x18f0206554638741,
            ]),
            c1: Fp::zero(),
        };

        G2Projective {
            // x = frobenius^2(x)/2^((p-1)/3); note that q^2 is the order of the field.
            x: self.x * psi2_coeff_x,
            // y = -frobenius^2(y); note that q^2 is the order of the field.
            y: self.y.neg(),
            // z = z
            z: self.z,
        }
    }

    /// Multiply `self` by `crate::BLS_X`, using double and add.
    #[cfg(feature = "endo")]
    fn mul_by_x(&self) -> G2Projective {
        let mut xself = G2Projective::identity();
        // NOTE: in BLS12-381 we can just skip the first bit.
        let mut x = crate::BLS_X >> 1;
        let mut acc = *self;
        while x != 0 {
            acc = acc.double();
            if x % 2 == 1 {
                xself += acc;
            }
            x >>= 1;
        }
        // finally, flip the sign
        if crate::BLS_X_IS_NEGATIVE {
            xself = -xself;
        }
        xself
    }

    /// Clears the cofactor, using [Budroni-Pintore](https://ia.cr/2017/419).
    /// This is equivalent to multiplying by $h\_\textrm{eff} = 3(z^2 - 1) \cdot
    /// h_2$, where $h_2$ is the cofactor of $\mathbb{G}\_2$ and $z$ is the
    /// parameter of BLS12-381.
    ///
    /// The endomorphism is only actually used if the crate feature `endo` is
    /// enabled, which it is by default.
    pub fn clear_cofactor(&self) -> G2Projective {
        #[cfg(feature = "endo")]
        fn clear_cofactor(this: &G2Projective) -> G2Projective {
            let t1 = this.mul_by_x(); // [x] P
            let t2 = this.psi(); // psi(P)

            this.double().psi2() // psi^2(2P)
                + (t1 + t2).mul_by_x() // psi^2(2P) + [x^2] P + [x] psi(P)
                - t1 // psi^2(2P) + [x^2 - x] P + [x] psi(P)
                - t2 // psi^2(2P) + [x^2 - x] P + [x - 1] psi(P)
                - this // psi^2(2P) + [x^2 - x - 1] P + [x - 1] psi(P)
        }

        #[cfg(not(feature = "endo"))]
        fn clear_cofactor(this: &G2Projective) -> G2Projective {
            this.multiply(&[
                0x51, 0x55, 0xa9, 0xaa, 0x5, 0x0, 0x2, 0xe8, 0xb4, 0xf6, 0xbb, 0xde, 0xa, 0x4c,
                0x89, 0x59, 0xa3, 0xf6, 0x89, 0x66, 0xc0, 0xcb, 0x54, 0xe9, 0x1a, 0x7c, 0x47, 0xd7,
                0x69, 0xec, 0xc0, 0x2e, 0xb0, 0x12, 0x12, 0x5d, 0x1, 0xbf, 0x82, 0x6d, 0x95, 0xdb,
                0x31, 0x87, 0x17, 0x2f, 0x9c, 0x32, 0xe1, 0xff, 0x8, 0x15, 0x3, 0xff, 0x86, 0x99,
                0x68, 0xd7, 0x5a, 0x14, 0xe9, 0xa8, 0xe2, 0x88, 0x28, 0x35, 0x1b, 0xa9, 0xe, 0x6a,
                0x4c, 0x58, 0xb3, 0x75, 0xee, 0xf2, 0x8, 0x9f, 0xc6, 0xb,
            ])
        }

        clear_cofactor(self)
    }

    /// Converts a batch of `G2Projective` elements into `G2Affine` elements. This
    /// function will panic if `p.len() != q.len()`.
    pub fn batch_normalize(p: &[Self], q: &mut [G2Affine]) {
        assert_eq!(p.len(), q.len());

        let mut acc = Fp2::one();
        for (p, q) in p.iter().zip(q.iter_mut()) {
            // We use the `x` field of `G2Affine` to store the product
            // of previous z-coordinates seen.
            q.x = acc;

            // We will end up skipping all identities in p
            acc = Fp2::conditional_select(&(acc * p.z), &acc, p.is_identity());
        }

        // This is the inverse, as all z-coordinates are nonzero and the ones
        // that are not are skipped.
        acc = acc.invert().unwrap();

        for (p, q) in p.iter().rev().zip(q.iter_mut().rev()) {
            let skip = p.is_identity();

            // Compute tmp = 1/z
            let tmp = q.x * acc;

            // Cancel out z-coordinate in denominator of `acc`
            acc = Fp2::conditional_select(&(acc * p.z), &acc, skip);

            // Set the coordinates to the correct value
            q.x = p.x * tmp;
            q.y = p.y * tmp;
            q.infinity = Choice::from(0u8);

            *q = G2Affine::conditional_select(&q, &G2Affine::identity(), skip);
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // Y^2 Z = X^3 + b Z^3

        (self.y.square() * self.z).ct_eq(&(self.x.square() * self.x + self.z.square() * self.z * B))
            | self.z.is_zero()
    }

    #[cfg(feature = "hashing")]
    /// Use a random oracle to map a value to a curve point
    /// TODO: Make public once it works
    fn hash<X>(msg: &[u8], dst: &[u8]) -> Self
    where
        X: ExpandMsg,
    {
        {
            let u = Fp2::hash::<X>(msg, dst);
            let q0 = Self::osswu_map(&u[0]);
            let q1 = Self::osswu_map(&u[1]);
            q0 + q1
        }
        .clear_cofactor()
    }

    #[cfg(feature = "hashing")]
    /// Use injective encoding to map a value to a curve point
    pub fn encode<M, D, X>(msg: M, dst: D) -> Self
    where
        M: AsRef<[u8]>,
        D: AsRef<[u8]>,
        X: ExpandMsg,
    {
        let u = Fp2::encode::<X>(msg.as_ref(), dst.as_ref());
        Self::osswu_map(&u).isogeny_map().clear_cofactor()
    }

    #[cfg(feature = "hashing")]
    /// Optimized simplified swu map for q = 9 mod 16 where AB == 0
    /// TODO: Still WIP
    fn osswu_map(u: &Fp2) -> Self {
        // Taken from section 8.8.2 in
        // <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html>
        const A: Fp2 = Fp2 {
            c0: Fp([
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
            ]),
            c1: Fp([
                0xe53a000003135242u64,
                0x01080c0fdef80285u64,
                0xe7889edbe340f6bdu64,
                0x0b51375126310601u64,
                0x02d6985717c744abu64,
                0x1220b4e979ea5467u64,
            ]),
        };
        const B: Fp2 = Fp2 {
            c0: Fp([
                0x22ea00000cf89db2u64,
                0x6ec832df71380aa4u64,
                0x6e1b94403db5a66eu64,
                0x75bf3c53a79473bau64,
                0x3dd3a569412c0a34u64,
                0x125cdb5e74dc4fd1u64,
            ]),
            c1: Fp([
                0x22ea00000cf89db2u64,
                0x6ec832df71380aa4u64,
                0x6e1b94403db5a66eu64,
                0x75bf3c53a79473bau64,
                0x3dd3a569412c0a34u64,
                0x125cdb5e74dc4fd1u64,
            ]),
        };
        const Z: Fp2 = Fp2 {
            c0: Fp([
                0x87ebfffffff9555cu64,
                0x656fffe5da8ffffau64,
                0xfd0749345d33ad2u64,
                0xd951e663066576f4u64,
                0xde291a3d41e980d3u64,
                0x815664c7dfe040du64,
            ]),
            c1: Fp([
                0x43f5fffffffcaaaeu64,
                0x32b7fff2ed47fffdu64,
                0x7e83a49a2e99d69u64,
                0xeca8f3318332bb7au64,
                0xef148d1ea0f4c069u64,
                0x40ab3263eff0206u64,
            ]),
        };
        const XD1: Fp2 = Fp2::mul(&Z, &A);
        const C2: Fp2 = Fp2 {
            c0: Fp([0, 0, 0, 0, 0, 0]),
            c1: Fp([
                0x43f5fffffffcaaaeu64,
                0x32b7fff2ed47fffdu64,
                0x07e83a49a2e99d69u64,
                0xeca8f3318332bb7au64,
                0xef148d1ea0f4c069u64,
                0x040ab3263eff0206u64,
            ]),
        };
        const C3: Fp2 = Fp2 {
            c0: Fp([
                0x7bcfa7a25aa30fdau64,
                0xdc17dec12a927e7cu64,
                0x2f088dd86b4ebef1u64,
                0xd1ca2087da74d4a7u64,
                0x2da2596696cebc1du64,
                0x0e2b7eedbbfd87d2u64,
            ]),
            c1: Fp([
                0x7bcfa7a25aa30fdau64,
                0xdc17dec12a927e7cu64,
                0x2f088dd86b4ebef1u64,
                0xd1ca2087da74d4a7u64,
                0x2da2596696cebc1du64,
                0x0e2b7eedbbfd87d2u64,
            ]),
        };
        const C4: Fp2 = Fp2 {
            c0: Fp([
                0x486f252db11dd19cu64,
                0x791ffda2c3d18950u64,
                0x5af6c27debf95eb4u64,
                0x73b1fd8f2a929cdeu64,
                0xfc59602a1a90b871u64,
                0x08d7daafa8baddb3u64,
            ]),
            c1: Fp([
                0xb8640a067f5c429fu64,
                0xcfd425f04b4dc505u64,
                0x072d7e2ebb535cb1u64,
                0xd947b5f9d2b4754du64,
                0x46a7142740774afbu64,
                0x0c31864c32fb3b7eu64,
            ]),
        };
        const C5: Fp2 = Fp2 {
            c0: Fp([
                0xb8640a067f5c429fu64,
                0xcfd425f04b4dc505u64,
                0x072d7e2ebb535cb1u64,
                0xd947b5f9d2b4754du64,
                0x46a7142740774afbu64,
                0x0c31864c32fb3b7eu64,
            ]),
            c1: Fp([
                0x718fdad24ee1d90fu64,
                0xa58c025bed8276afu64,
                0x0c3a10230ab7976fu64,
                0xf0c54df5c8f275e1u64,
                0x4ec2478c28baf465u64,
                0x1129373a90c508e6u64,
            ]),
        };

        // tv1 = u^2
        let mut tv1 = u.square();
        // tv3 = Z * tv1
        let tv3 = Z * tv1;
        // tv5 = tv3^2
        let mut tv5 = tv3.square();
        // xd = tv5 + tv3
        let mut xd = tv5 + tv3;
        // x1n = xd + 1
        // x1n = x1n * B
        let x1n = (xd + Fp2::one()) * B;
        // xd = -A * xd
        xd *= A.neg();
        // e1 = xd == 0
        // xd = CMOV(xd, Z * A, e1)
        xd.conditional_assign(&XD1, xd.is_zero());
        // tv2 = xd^2
        let mut tv2 = xd.square();
        let gxd = tv2 * xd;
        // tv2 = A * tv2
        tv2 *= A;
        // gx1 = x1n^2
        // gx1 = gx1 + tv2
        // gx1 = gx1 * x1n
        let mut gx1 = (x1n.square() + tv2) * x1n;
        tv2 = B * gxd;
        gx1 += tv2;
        let mut tv4 = gxd.square();
        tv2 = tv4 * gxd;
        tv4 = tv4.square();
        tv2 *= tv4;
        tv2 *= gx1;
        tv4 = tv4.square();
        tv4 *= tv2;
        let mut y = tv4;
        Self::chain_p2m9div16(&mut y, tv4);
        y *= tv2;
        tv4 = y * C2;
        tv2 = tv4.square();
        tv2 *= gxd;
        y.conditional_assign(&tv4, tv2.ct_eq(&gx1));
        tv4 = y * C3;
        tv2 = tv4.square();
        tv2 *= gxd;
        y.conditional_assign(&tv4, tv2.ct_eq(&gx1));
        tv4 *= C2;
        tv2 = tv4.square();
        tv2 *= gxd;
        y.conditional_assign(&tv4, tv2.ct_eq(&gx1));
        let gx2 = gx1 * tv5 * tv3;
        tv5 = y * tv1 * u;
        tv1 = tv5 * C4;
        tv4 = tv1 * C2;
        tv2 = tv4.square();
        tv2 *= gxd;
        tv1.conditional_assign(&tv4, tv2.ct_eq(&gx2));
        tv4 = tv5 * C5;
        tv2 = tv4.square();
        tv2 *= gxd;
        tv1.conditional_assign(&tv4, tv2.ct_eq(&gx2));
        tv4 *= C2;
        tv2 = tv4.square();
        tv2 *= gxd;
        tv1.conditional_assign(&tv4, tv2.ct_eq(&gx2));
        tv2 = y.square();
        tv2 *= gxd;
        let e8 = tv2.ct_eq(&gx1);
        y = Fp2::conditional_select(&tv1, &y, e8);
        tv2 = tv3 * x1n;
        let xn = Fp2::conditional_select(&tv2, &x1n, e8);
        let e9 = u.sgn0() ^ y.sgn0();
        y.negate_if(e9);

        Self {
            x: xn * xd.invert().unwrap(),
            y,
            z: Fp2::one(),
        }
    }

    #[cfg(feature = "hashing")]
    /// Computes the isogeny map for this point
    fn isogeny_map(&self) -> Self {
        use crate::isogeny::g2::*;

        fn compute(xxs: &[Fp2], k: &[Fp2]) -> Fp2 {
            let mut xx = Fp2::zero();
            for i in 0..k.len() {
                xx += xxs[i] * k[i];
            }
            xx
        }

        let mut xs = [Fp2::one(); 4];
        xs[1] = self.x;
        xs[2] = self.x.square();
        xs[3] = xs[2] * self.x;

        let x_num = compute(&xs, &XNUM);
        let x_den = compute(&xs, &XDEN);
        let y_num = compute(&xs, &YNUM);
        let y_den = compute(&xs, &YDEN);

        let x = x_num * x_den.invert().unwrap();
        let y = self.y * y_num * y_den.invert().unwrap();
        Self {
            x,
            y,
            z: Fp2::one(),
        }
    }

    /// Addition chain implementing exponentiation by (p**2 - 9) // 16.
    fn chain_p2m9div16(tmpvar1: &mut Fp2, tmpvar0: Fp2) {
        *tmpvar1 = tmpvar1.square();
        //Self::sqr(tmpvar1, tmpvar0);                              /*    0 : 2 */
        let mut tmpvar2 = *tmpvar1;
        tmpvar2.mul_assign(tmpvar0);
        //Self::mul(&mut tmpvar2, tmpvar1, tmpvar0);                /*    1 : 3 */
        let mut tmpvar15 = tmpvar2;
        tmpvar15.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar15, &tmpvar2, tmpvar1);              /*    2 : 5 */
        let mut tmpvar3 = tmpvar15;
        tmpvar3.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar3, &tmpvar15, tmpvar1);              /*    3 : 7 */
        let mut tmpvar14 = tmpvar3;
        tmpvar14.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar14, &tmpvar3, tmpvar1);              /*    4 : 9 */
        let mut tmpvar13 = tmpvar14;
        tmpvar13.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar13, &tmpvar14, tmpvar1);             /*    5 : 11 */
        let mut tmpvar5 = tmpvar13;
        tmpvar5.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar5, &tmpvar13, tmpvar1);              /*    6 : 13 */
        let mut tmpvar10 = tmpvar5;
        tmpvar10.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar10, &tmpvar5, tmpvar1);              /*    7 : 15 */
        let mut tmpvar9 = tmpvar10;
        tmpvar9.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar9, &tmpvar10, tmpvar1);              /*    8 : 17 */
        let mut tmpvar16 = tmpvar9;
        tmpvar16.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar16, &tmpvar9, tmpvar1);              /*    9 : 19 */
        let mut tmpvar4 = tmpvar16;
        tmpvar4.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar4, &tmpvar16, tmpvar1);              /*   10 : 21 */
        let mut tmpvar7 = tmpvar4;
        tmpvar7.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar7, &tmpvar4, tmpvar1);               /*   11 : 23 */
        let mut tmpvar6 = tmpvar7;
        tmpvar6.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar6, &tmpvar7, tmpvar1);               /*   12 : 25 */
        let mut tmpvar12 = tmpvar6;
        tmpvar12.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar12, &tmpvar6, tmpvar1);              /*   13 : 27 */
        let mut tmpvar8 = tmpvar12;
        tmpvar8.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar8, &tmpvar12, tmpvar1);              /*   14 : 29 */
        let mut tmpvar11 = tmpvar8;
        tmpvar11.mul_assign(*tmpvar1);
        //Self::mul(&mut tmpvar11, &tmpvar8, tmpvar1);              /*   15 : 31 */
        *tmpvar1 = tmpvar4;
        *tmpvar1 = tmpvar1.square();
        //Self::sqr(tmpvar1, &tmpvar4);                             /*   16 : 42 */
        for _ in 0..2 {
            *tmpvar1 = tmpvar1.square();
        } /*   17 : 168 */
        tmpvar1.mul_assign(tmpvar0); /*   19 : 169 */
        for _ in 0..9 {
            *tmpvar1 = tmpvar1.square();
        } /*   20 : 86528 */
        tmpvar1.mul_assign(&tmpvar12); /*   29 : 86555 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*   30 : 1384880 */
        tmpvar1.mul_assign(&tmpvar5); /*   34 : 1384893 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*   35 : 88633152 */
        tmpvar1.mul_assign(&tmpvar14); /*   41 : 88633161 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*   42 : 1418130576 */
        tmpvar1.mul_assign(&tmpvar3); /*   46 : 1418130583 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*   47 : 45380178656 */
        tmpvar1.mul_assign(&tmpvar2); /*   52 : 45380178659 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*   53 : 11617325736704 */
        tmpvar1.mul_assign(&tmpvar5); /*   61 : 11617325736717 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*   62 : 185877211787472 */
        tmpvar1.mul_assign(&tmpvar3); /*   66 : 185877211787479 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*   67 : 2974035388599664 */
        tmpvar1.mul_assign(&tmpvar10); /*   71 : 2974035388599679 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*   72 : 761353059481517824 */
        tmpvar1.mul_assign(&tmpvar8); /*   80 : 761353059481517853 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*   81 : 48726595806817142592 */
        tmpvar1.mul_assign(&tmpvar13); /*   87 : 48726595806817142603 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*   88 : 779625532909074281648 */
        tmpvar1.mul_assign(&tmpvar5); /*   92 : 779625532909074281661 */
        for _ in 0..3 {
            *tmpvar1 = tmpvar1.square();
        } /*   93 : 6237004263272594253288 */
        tmpvar1.mul_assign(tmpvar0); /*   96 : 6237004263272594253289 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*   97 : 399168272849446032210496 */
        tmpvar1.mul_assign(&tmpvar10); /*  103 : 399168272849446032210511 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  104 : 102187077849458184245890816 */
        tmpvar1.mul_assign(&tmpvar8); /*  112 : 102187077849458184245890845 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  113 : 6539972982365323791737014080 */
        tmpvar1.mul_assign(&tmpvar4); /*  119 : 6539972982365323791737014101 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  120 : 1674233083485522890684675609856 */
        tmpvar1.mul_assign(&tmpvar9); /*  128 : 1674233083485522890684675609873 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  129 : 53575458671536732501909619515936 */
        tmpvar1.mul_assign(&tmpvar10); /*  134 : 53575458671536732501909619515951 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  135 : 3428829354978350880122215649020864 */
        tmpvar1.mul_assign(&tmpvar14); /*  141 : 3428829354978350880122215649020873 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  142 : 109722539359307228163910900768667936 */
        tmpvar1.mul_assign(&tmpvar10); /*  147 : 109722539359307228163910900768667951 */
        for _ in 0..2 {
            *tmpvar1 = tmpvar1.square();
        } /*  148 : 438890157437228912655643603074671804 */
        tmpvar1.mul_assign(tmpvar0); /*  150 : 438890157437228912655643603074671805 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  151 : 28088970075982650409961190596778995520 */
        tmpvar1.mul_assign(&tmpvar10); /*  157 : 28088970075982650409961190596778995535 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  158 : 3595388169725779252475032396387711428480 */
        tmpvar1.mul_assign(&tmpvar13); /*  165 : 3595388169725779252475032396387711428491 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  166 : 57526210715612468039600518342203382855856 */
        tmpvar1.mul_assign(&tmpvar3); /*  170 : 57526210715612468039600518342203382855863 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  171 : 3681677485799197954534433173901016502775232 */
        tmpvar1.mul_assign(&tmpvar14); /*  177 : 3681677485799197954534433173901016502775241 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  178 : 471254718182297338180407446259330112355230848 */
        tmpvar1.mul_assign(&tmpvar3); /*  185 : 471254718182297338180407446259330112355230855 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  186 : 15080150981833514821773038280298563595367387360 */
        tmpvar1.mul_assign(&tmpvar15); /*  191 : 15080150981833514821773038280298563595367387365 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  192 : 1930259325674689897186948899878216140207025582720 */
        tmpvar1.mul_assign(&tmpvar3); /*  199 : 1930259325674689897186948899878216140207025582727 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  200 : 61768298421590076709982364796102916486624818647264 */
        tmpvar1.mul_assign(&tmpvar3); /*  205 : 61768298421590076709982364796102916486624818647271 */
        for _ in 0..10 {
            *tmpvar1 = tmpvar1.square();
        } /*  206 : 63250737583708238551021941551209386482303814294805504 */
        tmpvar1.mul_assign(&tmpvar9); /*  216 : 63250737583708238551021941551209386482303814294805521 */
        for _ in 0..3 {
            *tmpvar1 = tmpvar1.square();
        } /*  217 : 506005900669665908408175532409675091858430514358444168 */
        tmpvar1.mul_assign(&tmpvar15); /*  220 : 506005900669665908408175532409675091858430514358444173 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  221 : 16192188821429309069061617037109602939469776459470213536 */
        tmpvar1.mul_assign(&tmpvar5); /*  226 : 16192188821429309069061617037109602939469776459470213549 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  227 : 4145200338285903121679773961500058352504262773624374668544 */
        tmpvar1.mul_assign(&tmpvar6); /*  235 : 4145200338285903121679773961500058352504262773624374668569 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  236 : 132646410825148899893752766768001867280136408755979989394208 */
        tmpvar1.mul_assign(&tmpvar7); /*  241 : 132646410825148899893752766768001867280136408755979989394231 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  242 : 8489370292809529593200177073152119505928730160382719321230784 */
        tmpvar1.mul_assign(&tmpvar13); /*  248 : 8489370292809529593200177073152119505928730160382719321230795 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  249 : 543319698739809893964811332681735648379438730264494036558770880 */
        tmpvar1.mul_assign(&tmpvar10); /*  255 : 543319698739809893964811332681735648379438730264494036558770895 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  256 : 34772460719347833213747925291631081496284078736927618339761337280 */
        tmpvar1.mul_assign(&tmpvar14); /*  262 : 34772460719347833213747925291631081496284078736927618339761337289 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  263 : 4450874972076522651359734437328778431524362078326735147489451172992 */
        tmpvar1.mul_assign(&tmpvar16); /*  270 : 4450874972076522651359734437328778431524362078326735147489451173011 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  271 : 142427999106448724843511501994520909808779586506455524719662437536352 */
        tmpvar1.mul_assign(&tmpvar14); /*  276 : 142427999106448724843511501994520909808779586506455524719662437536361 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  277 : 9115391942812718389984736127649338227761893536413153582058396002327104 */
        tmpvar1.mul_assign(&tmpvar10); /*  283 : 9115391942812718389984736127649338227761893536413153582058396002327119 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  284 : 583385084340013976959023112169557646576761186330441829251737344148935616 */
        tmpvar1.mul_assign(&tmpvar9); /*  290 : 583385084340013976959023112169557646576761186330441829251737344148935633 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  291 : 18668322698880447262688739589425844690456357962574138536055595012765940256 */
        tmpvar1.mul_assign(&tmpvar10); /*  296 : 18668322698880447262688739589425844690456357962574138536055595012765940271 */
        for _ in 0..2 {
            *tmpvar1 = tmpvar1.square();
        } /*  297 : 74673290795521789050754958357703378761825431850296554144222380051063761084 */
        tmpvar1.mul_assign(tmpvar0); /*  299 : 74673290795521789050754958357703378761825431850296554144222380051063761085 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  300 : 19116362443653577996993269339572064963027310553675917860920929293072322837760 */
        tmpvar1.mul_assign(&tmpvar15); /*  308 : 19116362443653577996993269339572064963027310553675917860920929293072322837765 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  309 : 2446894392787657983615138475465224315267495750870517486197878949513257323233920 */
        tmpvar1.mul_assign(&tmpvar15); /*  316 : 2446894392787657983615138475465224315267495750870517486197878949513257323233925 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  317 : 39150310284602527737842215607443589044279932013928279779166063192212117171742800 */
        tmpvar1.mul_assign(&tmpvar2); /*  321 : 39150310284602527737842215607443589044279932013928279779166063192212117171742803 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  322 : 5011239716429123550443803597752779397667831297782819811733256088603150997983078784 */
        tmpvar1.mul_assign(&tmpvar13); /*  329 : 5011239716429123550443803597752779397667831297782819811733256088603150997983078795 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  330 : 320719341851463907228403430256177881450741203058100467950928389670601663870917042880 */
        tmpvar1.mul_assign(&tmpvar10); /*  336 : 320719341851463907228403430256177881450741203058100467950928389670601663870917042895 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  337 : 5131509469623422515654454884098846103211859248929607487214854234729626621934672686320 */
        tmpvar1.mul_assign(&tmpvar5); /*  341 : 5131509469623422515654454884098846103211859248929607487214854234729626621934672686333 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  342 : 656833212111798082003770225164652301211117983862989758363501342045392207607638103850624 */
        tmpvar1.mul_assign(&tmpvar13); /*  349 : 656833212111798082003770225164652301211117983862989758363501342045392207607638103850635 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  350 : 42037325575155077248241294410537747277511550967231344535264085890905101286888838646440640 */
        tmpvar1.mul_assign(&tmpvar12); /*  356 : 42037325575155077248241294410537747277511550967231344535264085890905101286888838646440667 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  357 : 1345194418404962471943721421137207912880369630951403025128450748508963241180442836686101344 */
        tmpvar1.mul_assign(&tmpvar7); /*  362 : 1345194418404962471943721421137207912880369630951403025128450748508963241180442836686101367 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  363 : 43046221388958799102199085476390653212171828190444896804110423952286823717774170773955243744 */
        tmpvar1.mul_assign(&tmpvar15); /*  368 : 43046221388958799102199085476390653212171828190444896804110423952286823717774170773955243749 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  369 : 5509916337786726285081482940978003611157994008376946790926134265892713435875093859066271199872 */
        tmpvar1.mul_assign(&tmpvar12); /*  376 : 5509916337786726285081482940978003611157994008376946790926134265892713435875093859066271199899 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  377 : 176317322809175241122607454111296115557055808268062297309636296508566829948003003490120678396768 */
        tmpvar1.mul_assign(&tmpvar7); /*  382 : 176317322809175241122607454111296115557055808268062297309636296508566829948003003490120678396791 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  383 : 5642154329893607715923438531561475697825785864577993513908361488274138558336096111683861708697312 */
        tmpvar1.mul_assign(&tmpvar4); /*  388 : 5642154329893607715923438531561475697825785864577993513908361488274138558336096111683861708697333 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  389 : 90274469278297723454775016504983611165212573833247896222533783812386216933377537786941787339157328 */
        tmpvar1.mul_assign(&tmpvar2); /*  393 : 90274469278297723454775016504983611165212573833247896222533783812386216933377537786941787339157331 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  394 : 5777566033811054301105601056318951114573604725327865358242162163992717883736162418364274389706069184 */
        tmpvar1.mul_assign(&tmpvar15); /*  400 : 5777566033811054301105601056318951114573604725327865358242162163992717883736162418364274389706069189 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  401 : 369764226163907475270758467604412871332710702420983382927498378495533944559114394775313560941188428096 */
        tmpvar1.mul_assign(&tmpvar14); /*  407 : 369764226163907475270758467604412871332710702420983382927498378495533944559114394775313560941188428105 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  408 : 5916227618622519604332135481670605941323371238735734126839974055928543112945830316405016975059014849680 */
        tmpvar1.mul_assign(&tmpvar2); /*  412 : 5916227618622519604332135481670605941323371238735734126839974055928543112945830316405016975059014849683 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  413 : 94659641897960313669314167706729695061173939819771746029439584894856689807133285062480271600944237594928 */
        tmpvar1.mul_assign(&tmpvar2); /*  417 : 94659641897960313669314167706729695061173939819771746029439584894856689807133285062480271600944237594931 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  418 : 24232868325877840299344426932922801935660528593861566983536533733083312590626120975994949529841724824302336 */
        tmpvar1.mul_assign(&tmpvar14); /*  426 : 24232868325877840299344426932922801935660528593861566983536533733083312590626120975994949529841724824302345 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  427 : 775451786428090889579021661853529661941136915003570143473169079458666002900035871231838384954935194377675040 */
        tmpvar1.mul_assign(&tmpvar10); /*  432 : 775451786428090889579021661853529661941136915003570143473169079458666002900035871231838384954935194377675055 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  433 : 49628914331397816933057386358625898364232762560228489182282821085354624185602295758837656637115852440171203520 */
        tmpvar1.mul_assign(&tmpvar3); /*  439 : 49628914331397816933057386358625898364232762560228489182282821085354624185602295758837656637115852440171203527 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  440 : 1588125258604730141857836363476028747655448401927311653833050274731347973939273464282805012387707278085478512864 */
        tmpvar1.mul_assign(&tmpvar10); /*  445 : 1588125258604730141857836363476028747655448401927311653833050274731347973939273464282805012387707278085478512879 */
        for _ in 0..12 {
            *tmpvar1 = tmpvar1.square();
        } /*  446 : 6504961059244974661049697744797813750396716654294268534100173925299601301255264109702369330740049011038119988752384 */
        tmpvar1.mul_assign(&tmpvar9); /*  458 : 6504961059244974661049697744797813750396716654294268534100173925299601301255264109702369330740049011038119988752401 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  459 : 104079376947919594576795163916765020006347466468708296545602782804793620820084225755237909291840784176609919820038416 */
        tmpvar1.mul_assign(&tmpvar5); /*  463 : 104079376947919594576795163916765020006347466468708296545602782804793620820084225755237909291840784176609919820038429 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  464 : 3330540062333427026457445245336480640203118926998665489459289049753395866242695224167613097338905093651517434241229728 */
        tmpvar1.mul_assign(&tmpvar5); /*  469 : 3330540062333427026457445245336480640203118926998665489459289049753395866242695224167613097338905093651517434241229741 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  470 : 213154563989339329693276495701534760972999611327914591325394499184217335439532494346727238229689925993697115791438703424 */
        tmpvar1.mul_assign(&tmpvar2); /*  476 : 213154563989339329693276495701534760972999611327914591325394499184217335439532494346727238229689925993697115791438703427 */
        for _ in 0..9 {
            *tmpvar1 = tmpvar1.square();
        } /*  477 : 109135136762541736802957565799185797618175800999892270758601983582319275745040637105524345973601242108772923285216616154624 */
        tmpvar1.mul_assign(&tmpvar6); /*  486 : 109135136762541736802957565799185797618175800999892270758601983582319275745040637105524345973601242108772923285216616154649 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  487 : 3492324376401335577694642105573945523781625631996552664275263474634216823841300387376779071155239747480733545126931716948768 */
        tmpvar1.mul_assign(&tmpvar6); /*  492 : 3492324376401335577694642105573945523781625631996552664275263474634216823841300387376779071155239747480733545126931716948793 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  493 : 223508760089685476972457094756732513522024040447779370513616862376589876725843224792113860553935343838766946888123629884722752 */
        tmpvar1.mul_assign(&tmpvar2); /*  499 : 223508760089685476972457094756732513522024040447779370513616862376589876725843224792113860553935343838766946888123629884722755 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  500 : 14304560645739870526237254064430880865409538588657879712871479192101752110453966386695287075451862005681084600839912312622256320 */
        tmpvar1.mul_assign(&tmpvar2); /*  506 : 14304560645739870526237254064430880865409538588657879712871479192101752110453966386695287075451862005681084600839912312622256323 */
        for _ in 0..9 {
            *tmpvar1 = tmpvar1.square();
        } /*  507 : 7323935050618813709433474080988611003089683757392834412990197346356097080552430789987986982631353346908715315630035104062595237376 */
        tmpvar1.mul_assign(&tmpvar7); /*  516 : 7323935050618813709433474080988611003089683757392834412990197346356097080552430789987986982631353346908715315630035104062595237399 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  517 : 937463686479208154807484682366542208395479520946282804862745260333580426310711141118462333776813228404315560400644493320012190387072 */
        tmpvar1.mul_assign(&tmpvar10); /*  524 : 937463686479208154807484682366542208395479520946282804862745260333580426310711141118462333776813228404315560400644493320012190387087 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  525 : 59997675934669321907679019671458701337310689340562099511215696661349147283885513031581589361716046617876195865641247572480780184773568 */
        tmpvar1.mul_assign(&tmpvar6); /*  531 : 59997675934669321907679019671458701337310689340562099511215696661349147283885513031581589361716046617876195865641247572480780184773593 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  532 : 1919925629909418301045728629486678442793942058897987184358902293163172713084336417010610859574913491772038267700519922319384965912754976 */
        tmpvar1.mul_assign(&tmpvar14); /*  537 : 1919925629909418301045728629486678442793942058897987184358902293163172713084336417010610859574913491772038267700519922319384965912754985 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  538 : 245750480628405542533853264574294840677624583538942359597939493524886107274795061377358190025588926946820898265666550056881275636832638080 */
        tmpvar1.mul_assign(&tmpvar7); /*  545 : 245750480628405542533853264574294840677624583538942359597939493524886107274795061377358190025588926946820898265666550056881275636832638103 */
        for _ in 0..2 {
            *tmpvar1 = tmpvar1.square();
        } /*  546 : 983001922513622170135413058297179362710498334155769438391757974099544429099180245509432760102355707787283593062666200227525102547330552412 */
        tmpvar1.mul_assign(tmpvar0); /*  548 : 983001922513622170135413058297179362710498334155769438391757974099544429099180245509432760102355707787283593062666200227525102547330552413 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  549 : 251648492163487275554665742924077916853887573543876976228290041369483373849390142850414786586203061193544599824042547258246426252116621417728 */
        tmpvar1.mul_assign(&tmpvar13); /*  557 : 251648492163487275554665742924077916853887573543876976228290041369483373849390142850414786586203061193544599824042547258246426252116621417739 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  558 : 4026375874615796408874651886785246669662201176702031619652640661911733981590242285606636585379248979096713597184680756131942820033865942683824 */
        tmpvar1.mul_assign(&tmpvar15); /*  562 : 4026375874615796408874651886785246669662201176702031619652640661911733981590242285606636585379248979096713597184680756131942820033865942683829 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  563 : 515376111950821940335955441508511573716761750617860047315538004724701949643551012557649482928543869324379340439639136784888680964334840663530112 */
        tmpvar1.mul_assign(&tmpvar3); /*  570 : 515376111950821940335955441508511573716761750617860047315538004724701949643551012557649482928543869324379340439639136784888680964334840663530119 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  571 : 131936284659410416726004593026178962871491008158172172112777729209523699108749059214758267629707230547041111152547619016931502326869719209863710464 */
        tmpvar1.mul_assign(&tmpvar14); /*  579 : 131936284659410416726004593026178962871491008158172172112777729209523699108749059214758267629707230547041111152547619016931502326869719209863710473 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  580 : 16887844436404533340928587907350907247550849044246038030435549338819033485919879579489058256602525510021262227526095234167232297839324058862554940544 */
        tmpvar1.mul_assign(&tmpvar5); /*  587 : 16887844436404533340928587907350907247550849044246038030435549338819033485919879579489058256602525510021262227526095234167232297839324058862554940557 */
        for _ in 0..10 {
            *tmpvar1 = tmpvar1.square();
        } /*  588 : 17293152702878242141110874017127329021492069421307942943166002522950690289581956689396795654760986122261772520986721519787245872987467836275256259130368 */
        tmpvar1.mul_assign(&tmpvar14); /*  598 : 17293152702878242141110874017127329021492069421307942943166002522950690289581956689396795654760986122261772520986721519787245872987467836275256259130377 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  599 : 1106761772984207497031095937096149057375492442963708348362624161468844178533245228121394921904703111824753441343150177266383735871197941521616400584344128 */
        tmpvar1.mul_assign(&tmpvar13); /*  605 : 1106761772984207497031095937096149057375492442963708348362624161468844178533245228121394921904703111824753441343150177266383735871197941521616400584344139 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  606 : 70832753470989279809990139974153539672031516349677334295207946334006027426127694599769275001900999156784220245961611345048559095756668257383449637398024896 */
        tmpvar1.mul_assign(&tmpvar5); /*  612 : 70832753470989279809990139974153539672031516349677334295207946334006027426127694599769275001900999156784220245961611345048559095756668257383449637398024909 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  613 : 4533296222143313907839368958345826539010017046379349394893308565376385755272172454385233600121663946034190095741543126083107782128426768472540776793473594176 */
        tmpvar1.mul_assign(&tmpvar11); /*  619 : 4533296222143313907839368958345826539010017046379349394893308565376385755272172454385233600121663946034190095741543126083107782128426768472540776793473594207 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  620 : 145065479108586045050859806667066449248320545484139180636585874092044344168709518540327475203893246273094083063729380034659449028109656591121304857391155014624 */
        tmpvar1.mul_assign(&tmpvar6); /*  625 : 145065479108586045050859806667066449248320545484139180636585874092044344168709518540327475203893246273094083063729380034659449028109656591121304857391155014649 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  626 : 18568381325899013766510055253384505503785029821969815121482991883781676053594818373161916826098335522956042632157360644436409475598036043663527021746067841875072 */
        tmpvar1.mul_assign(&tmpvar10); /*  633 : 18568381325899013766510055253384505503785029821969815121482991883781676053594818373161916826098335522956042632157360644436409475598036043663527021746067841875087 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  634 : 594188202428768440528321768108304176121120954303034083887455740281013633715034187941181338435146736734593364229035540621965103219137153397232864695874170940002784 */
        tmpvar1.mul_assign(&tmpvar5); /*  639 : 594188202428768440528321768108304176121120954303034083887455740281013633715034187941181338435146736734593364229035540621965103219137153397232864695874170940002797 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  640 : 76056089910882360387625186317862934543503482150788362737594334755969745115524376056471211319698782302027950621316549199611533212049555634845806681071893880320358016 */
        tmpvar1.mul_assign(&tmpvar11); /*  647 : 76056089910882360387625186317862934543503482150788362737594334755969745115524376056471211319698782302027950621316549199611533212049555634845806681071893880320358047 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  648 : 2433794877148235532404005962171613905392111428825227607603018712191031843696780033807078762230361033664894419882129574387569062785585780315065813794300604170251457504 */
        tmpvar1.mul_assign(&tmpvar3); /*  653 : 2433794877148235532404005962171613905392111428825227607603018712191031843696780033807078762230361033664894419882129574387569062785585780315065813794300604170251457511 */
        for _ in 0..8 {
            *tmpvar1 = tmpvar1.square();
        } /*  654 : 623051488549948296295425526315933159780380525779258267546372790320904151986375688654612163130972424618212971489825171043217680073109959760656848331340954667584373122816 */
        tmpvar1.mul_assign(&tmpvar12); /*  662 : 623051488549948296295425526315933159780380525779258267546372790320904151986375688654612163130972424618212971489825171043217680073109959760656848331340954667584373122843 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  663 : 39875295267196690962907233684219722225944353649872529122967858580537865727128044073895178440382235175565630175348810946765931524679037424682038293205821098725399879861952 */
        tmpvar1.mul_assign(&tmpvar8); /*  669 : 39875295267196690962907233684219722225944353649872529122967858580537865727128044073895178440382235175565630175348810946765931524679037424682038293205821098725399879861981 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  670 : 2552018897100588221626062955790062222460438633591841863869942949154423406536194820729291420184463051236200331222323900593019617579458395179650450765172550318425592311166784 */
        tmpvar1.mul_assign(&tmpvar2); /*  676 : 2552018897100588221626062955790062222460438633591841863869942949154423406536194820729291420184463051236200331222323900593019617579458395179650450765172550318425592311166787 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  677 : 326658418828875292368136058341127964474936145099755758575352697491766196036632937053349301783611270558233642396457459275906511050170674582995257697942086440758475815829348736 */
        tmpvar1.mul_assign(&tmpvar13); /*  684 : 326658418828875292368136058341127964474936145099755758575352697491766196036632937053349301783611270558233642396457459275906511050170674582995257697942086440758475815829348747 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  685 : 41812277610096037423121415467664379452791826572768737097645145278946073092689015942828710628302242631453906226746554787316033414421846346623392985336587064417084904426156639616 */
        tmpvar1.mul_assign(&tmpvar13); /*  692 : 41812277610096037423121415467664379452791826572768737097645145278946073092689015942828710628302242631453906226746554787316033414421846346623392985336587064417084904426156639627 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  693 : 2675985767046146395079770589930520284978676900657199174249289297852548677932097020341037480211343528413049998511779506388226138522998166183897151061541572122693433883274024936128 */
        tmpvar1.mul_assign(&tmpvar2); /*  699 : 2675985767046146395079770589930520284978676900657199174249289297852548677932097020341037480211343528413049998511779506388226138522998166183897151061541572122693433883274024936131 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  700 : 85631544545476684642552658877776649119317660821030373575977257531281557693827104650913199366762992909217599952376944204423236432735941317884708833969330307926189884264768797956192 */
        tmpvar1.mul_assign(&tmpvar3); /*  705 : 85631544545476684642552658877776649119317660821030373575977257531281557693827104650913199366762992909217599952376944204423236432735941317884708833969330307926189884264768797956199 */
        for _ in 0..10 {
            *tmpvar1 = tmpvar1.square();
        } /*  706 : 87686701614568125073973922690843288698181284680735102541800711712032315078478955162535116151565304739038822351233990865329394107121603909513941845984594235316418441487123249107147776 */
        tmpvar1.mul_assign(&tmpvar12); /*  716 : 87686701614568125073973922690843288698181284680735102541800711712032315078478955162535116151565304739038822351233990865329394107121603909513941845984594235316418441487123249107147803 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  717 : 1402987225833090001183582763053492619170900554891761640668811387392517041255663282600561858425044875824621157619743853845270305713945662552223069535753507765062695063793971985714364848 */
        tmpvar1.mul_assign(tmpvar0); /*  721 : 1402987225833090001183582763053492619170900554891761640668811387392517041255663282600561858425044875824621157619743853845270305713945662552223069535753507765062695063793971985714364849 */
        for _ in 0..9 {
            *tmpvar1 = tmpvar1.square();
        } /*  722 : 718329459626542080605994374683388221015501084104581960022431430344968725122899600691487671513622976422206032701308853168778396525540179226738211602305795975712099872662513656685754802688 */
        tmpvar1.mul_assign(&tmpvar9); /*  731 : 718329459626542080605994374683388221015501084104581960022431430344968725122899600691487671513622976422206032701308853168778396525540179226738211602305795975712099872662513656685754802705 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  732 : 45973085416098693158783639979736846144992069382693245441435611542077998407865574444255210976871870491021186092883766602801817377634571470511245542547570942445574391850400874027888307373120 */
        tmpvar1.mul_assign(&tmpvar10); /*  738 : 45973085416098693158783639979736846144992069382693245441435611542077998407865574444255210976871870491021186092883766602801817377634571470511245542547570942445574391850400874027888307373135 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  739 : 5884554933260632724324305917406316306558984880984735416503758277385983796206793528864667005039599422850711819889122125158632624337225148225439429446089080633033522156851311875569703343761280 */
        tmpvar1.mul_assign(&tmpvar11); /*  746 : 5884554933260632724324305917406316306558984880984735416503758277385983796206793528864667005039599422850711819889122125158632624337225148225439429446089080633033522156851311875569703343761311 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  747 : 188305757864340247178377789357002121809887516191511533328120264876351481478617392923669344161267181531222778236451908005076243978791204743214061742274850580257072709019241980018230507000361952 */
        tmpvar1.mul_assign(&tmpvar4); /*  752 : 188305757864340247178377789357002121809887516191511533328120264876351481478617392923669344161267181531222778236451908005076243978791204743214061742274850580257072709019241980018230507000361973 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  753 : 3012892125829443954854044629712033948958200259064184533249924238021623703657878286778709506580274904499564451783230528081219903660659275891424987876397609284113163344307871680291688112005791568 */
        tmpvar1.mul_assign(&tmpvar10); /*  757 : 3012892125829443954854044629712033948958200259064184533249924238021623703657878286778709506580274904499564451783230528081219903660659275891424987876397609284113163344307871680291688112005791583 */
        for _ in 0..7 {
            *tmpvar1 = tmpvar1.square();
        } /*  758 : 385650192106168826221317712603140345466649633160215620255990302466767834068208420707674816842275187775944249828253507594396147668564387314102398448178893988366484908071407575077336078336741322624 */
        tmpvar1.mul_assign(&tmpvar8); /*  765 : 385650192106168826221317712603140345466649633160215620255990302466767834068208420707674816842275187775944249828253507594396147668564387314102398448178893988366484908071407575077336078336741322653 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  766 : 12340806147397402439082166803300491054932788261126899848191689678936570690182669462645594138952806008830215994504112243020676725394060394051276750341724607627727517058285042402474754506775722324896 */
        tmpvar1.mul_assign(&tmpvar4); /*  771 : 12340806147397402439082166803300491054932788261126899848191689678936570690182669462645594138952806008830215994504112243020676725394060394051276750341724607627727517058285042402474754506775722324917 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  772 : 394905796716716878050629337705615713757849224356060795142134069725970262085845422804659012446489792282566911824131591776661655212609932609640856010935187444087280545865121356879192144216823114397344 */
        tmpvar1.mul_assign(&tmpvar4); /*  777 : 394905796716716878050629337705615713757849224356060795142134069725970262085845422804659012446489792282566911824131591776661655212609932609640856010935187444087280545865121356879192144216823114397365 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  778 : 12636985494934940097620138806579702840251175179393945444548290231231048386747053529749088398287673353042141178372210936853172966803517843508507392349925998210792977467683883420134148614938339660715680 */
        tmpvar1.mul_assign(&tmpvar9); /*  783 : 12636985494934940097620138806579702840251175179393945444548290231231048386747053529749088398287673353042141178372210936853172966803517843508507392349925998210792977467683883420134148614938339660715697 */
        for _ in 0..4 {
            *tmpvar1 = tmpvar1.square();
        } /*  784 : 202191767918959041561922220905275245444018802870303127112772643699696774187952856475985414372602773648674258853955374989650767468856285496136118277598815971372687639482942134722146377839013434571451152 */
        tmpvar1.mul_assign(&tmpvar5); /*  788 : 202191767918959041561922220905275245444018802870303127112772643699696774187952856475985414372602773648674258853955374989650767468856285496136118277598815971372687639482942134722146377839013434571451165 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  789 : 12940273146813378659963022137937615708417203383699400135217449196780593548028982814463066519846577513515152566653143999337649118006802271752711569766324222167852008926908296622217368181696859812572874560 */
        tmpvar1.mul_assign(&tmpvar8); /*  795 : 12940273146813378659963022137937615708417203383699400135217449196780593548028982814463066519846577513515152566653143999337649118006802271752711569766324222167852008926908296622217368181696859812572874589 */
        *tmpvar1 = tmpvar1.square(); /*  796 : 25880546293626757319926044275875231416834406767398800270434898393561187096057965628926133039693155027030305133306287998675298236013604543505423139532648444335704017853816593244434736363393719625145749178 */
        tmpvar1.mul_assign(tmpvar0); /*  797 : 25880546293626757319926044275875231416834406767398800270434898393561187096057965628926133039693155027030305133306287998675298236013604543505423139532648444335704017853816593244434736363393719625145749179 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  798 : 1656354962792112468475266833656014810677402033113523217307833497187915974147709800251272514540361921729939528531602431915219087104870690784347080930089500437485057142644261967643823127257198056009327947456 */
        tmpvar1.mul_assign(&tmpvar3); /*  804 : 1656354962792112468475266833656014810677402033113523217307833497187915974147709800251272514540361921729939528531602431915219087104870690784347080930089500437485057142644261967643823127257198056009327947463 */
        for _ in 0..10 {
            *tmpvar1 = tmpvar1.square();
        } /*  805 : 1696107481899123167718673237663759166133659681908247774523221501120425957527254835457303054889330607851458077216360890281184345195387587363171410872411648447984698514067724254867274882311370809353551818202112 */
        tmpvar1.mul_assign(&tmpvar7); /*  815 : 1696107481899123167718673237663759166133659681908247774523221501120425957527254835457303054889330607851458077216360890281184345195387587363171410872411648447984698514067724254867274882311370809353551818202135 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  816 : 108550878841543882733995087210480586632554219642127857569486176071707261281744309469267395512917158902493316941847096977995798092504805591242970295834345500671020704900334352311505592467927731798627316364936640 */
        tmpvar1.mul_assign(&tmpvar4); /*  822 : 108550878841543882733995087210480586632554219642127857569486176071707261281744309469267395512917158902493316941847096977995798092504805591242970295834345500671020704900334352311505592467927731798627316364936661 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  823 : 6947256245858808494975685581470757544483470057096182884447115268589264722031635806033113312826698169759572284278214206591731077920307557839550098933398112042945325113621398547936357917947374835112148247355946304 */
        tmpvar1.mul_assign(&tmpvar6); /*  829 : 6947256245858808494975685581470757544483470057096182884447115268589264722031635806033113312826698169759572284278214206591731077920307557839550098933398112042945325113621398547936357917947374835112148247355946329 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  830 : 444624399734963743678443877214128482846942083654155704604615377189712942210024691586119252020908682864612626193805709221870788986899683701731206331737479170748500807271769507067926906748631989447177487830780565056 */
        tmpvar1.mul_assign(&tmpvar5); /*  836 : 444624399734963743678443877214128482846942083654155704604615377189712942210024691586119252020908682864612626193805709221870788986899683701731206331737479170748500807271769507067926906748631989447177487830780565069 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  837 : 28455961583037679595420408141704222902204293353865965094695384140141628301441580261511632129338155703335208076403565390199730495161579756910797205231198666927904051665393248452347322031912447324619359221169956164416 */
        tmpvar1.mul_assign(&tmpvar4); /*  843 : 28455961583037679595420408141704222902204293353865965094695384140141628301441580261511632129338155703335208076403565390199730495161579756910797205231198666927904051665393248452347322031912447324619359221169956164437 */
        for _ in 0..23 {
            *tmpvar1 = tmpvar1.square();
        } /*  844 : 238705906983162543355580399100765177871214152862586865721082456961065184302499251714358569373223087638243353151383559860752580829556389241459968722180074986980751351032731127113348364375477010926880553717580063640645533696 */
        tmpvar1.mul_assign(&tmpvar3); /*  867 : 238705906983162543355580399100765177871214152862586865721082456961065184302499251714358569373223087638243353151383559860752580829556389241459968722180074986980751351032731127113348364375477010926880553717580063640645533703 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  868 : 15277178046922402774757145542448971383757705783205559406149277245508171795359952109718948439886277608847574601688547831088165173091608911453437998219524799166768086466094792135254295320030528699320355437925124073001314156992 */
        tmpvar1.mul_assign(&tmpvar3); /*  874 : 15277178046922402774757145542448971383757705783205559406149277245508171795359952109718948439886277608847574601688547831088165173091608911453437998219524799166768086466094792135254295320030528699320355437925124073001314156999 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  875 : 488869697501516888792228657358367084280246585062577900996776871856261497451518467511006350076360883483122387254033530594821285538931485166510015943024793573336578766915033348328137450240976918378251374013603970336042053023968 */
        tmpvar1.mul_assign(&tmpvar2); /*  880 : 488869697501516888792228657358367084280246585062577900996776871856261497451518467511006350076360883483122387254033530594821285538931485166510015943024793573336578766915033348328137450240976918378251374013603970336042053023971 */
        for _ in 0..6 {
            *tmpvar1 = tmpvar1.square();
        } /*  881 : 31287660640097080882702634070935493393935781444004985663793719798800735836897181920704406404887096542919832784258145958068562274491615050656641020353586788693541041082562134293000796815422522776208087936870654101506691393534144 */
        tmpvar1.mul_assign(&tmpvar3); /*  887 : 31287660640097080882702634070935493393935781444004985663793719798800735836897181920704406404887096542919832784258145958068562274491615050656641020353586788693541041082562134293000796815422522776208087936870654101506691393534151 */
        for _ in 0..5 {
            *tmpvar1 = tmpvar1.square();
        } /*  888 : 1001205140483106588246484290269935788605945006208159541241399033561623546780709821462541004956387089373434649096260670658193992783731681621012512651314777238193313314641988297376025498093520728838658813979860931248214124593092832 */
        tmpvar1.mul_assign(&tmpvar2); /*  893 : 1001205140483106588246484290269935788605945006208159541241399033561623546780709821462541004956387089373434649096260670658193992783731681621012512651314777238193313314641988297376025498093520728838658813979860931248214124593092835 */
    }

    impl_pippenger_sum_of_products!();
}

pub struct G2Compressed([u8; 96]);

impl fmt::Debug for G2Compressed {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G2Compressed {
    fn default() -> Self {
        G2Compressed([0; 96])
    }
}

impl AsRef<[u8]> for G2Compressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G2Compressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

pub struct G2Uncompressed([u8; 192]);

impl fmt::Debug for G2Uncompressed {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G2Uncompressed {
    fn default() -> Self {
        G2Uncompressed([0; 192])
    }
}

impl AsRef<[u8]> for G2Uncompressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G2Uncompressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl Group for G2Projective {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        loop {
            let x = Fp2::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + 4)
            let p = ((x.square() * x) + B).sqrt().map(|y| G2Affine {
                x,
                y: if flip_sign { -y } else { y },
                infinity: 0.into(),
            });

            if p.is_some().into() {
                let p = p.unwrap().to_curve().clear_cofactor();

                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    fn identity() -> Self {
        Self::identity()
    }

    fn generator() -> Self {
        Self::generator()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }
}

#[cfg(feature = "alloc")]
impl WnafGroup for G2Projective {
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        const RECOMMENDATIONS: [usize; 11] = [1, 3, 8, 20, 47, 126, 260, 826, 1501, 4555, 84071];

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

impl PrimeGroup for G2Projective {}

impl Curve for G2Projective {
    type AffineRepr = G2Affine;

    fn batch_normalize(p: &[Self], q: &mut [Self::AffineRepr]) {
        Self::batch_normalize(p, q);
    }

    fn to_affine(&self) -> Self::AffineRepr {
        self.into()
    }
}

impl PrimeCurve for G2Projective {
    type Affine = G2Affine;
}

impl PrimeCurveAffine for G2Affine {
    type Scalar = Scalar;
    type Curve = G2Projective;

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
        G2Affine::from_bytes(bytes).map(Self::from)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        G2Affine::from_bytes_unchecked(bytes).map(Self::from)
    }

    fn to_bytes(&self) -> Self::Repr {
        G2Affine::from(self).to_bytes()
    }
}

impl GroupEncoding for G2Affine {
    type Repr = G2Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_compressed(&bytes.0)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_compressed_unchecked(&bytes.0)
    }

    fn to_bytes(&self) -> Self::Repr {
        G2Compressed(self.to_compressed())
    }
}

impl UncompressedEncoding for G2Affine {
    type Uncompressed = G2Uncompressed;

    fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
        Self::from_uncompressed(&bytes.0)
    }

    fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(&bytes.0)
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        G2Uncompressed(self.to_uncompressed())
    }
}

impl serde::Serialize for G2Projective {
    fn serialize<S>(&self, s: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.to_affine().serialize(s)
    }
}

impl<'de> serde::Deserialize<'de> for G2Projective {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(G2Projective::from(G2Affine::deserialize(deserializer)?))
    }
}

#[test]
fn test_is_on_curve() {
    assert!(bool::from(G2Affine::identity().is_on_curve()));
    assert!(bool::from(G2Affine::generator().is_on_curve()));
    assert!(bool::from(G2Projective::identity().is_on_curve()));
    assert!(bool::from(G2Projective::generator().is_on_curve()));

    let z = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xba7a_fa1f_9a6f_e250,
            0xfa0f_5b59_5eaf_e731,
            0x3bdc_4776_94c3_06e7,
            0x2149_be4b_3949_fa24,
            0x64aa_6e06_49b2_078c,
            0x12b1_08ac_3364_3c3e,
        ]),
        c1: Fp::from_raw_unchecked([
            0x1253_25df_3d35_b5a8,
            0xdc46_9ef5_555d_7fe3,
            0x02d7_16d2_4431_06a9,
            0x05a1_db59_a6ff_37d0,
            0x7cf7_784e_5300_bb8f,
            0x16a8_8922_c7a5_e844,
        ]),
    };

    let gen = G2Affine::generator();
    let mut test = G2Projective {
        x: gen.x * z,
        y: gen.y * z,
        z,
    };

    assert!(bool::from(test.is_on_curve()));

    test.x = z;
    assert!(!bool::from(test.is_on_curve()));
}

#[test]
#[allow(clippy::eq_op)]
fn test_affine_point_equality() {
    let a = G2Affine::generator();
    let b = G2Affine::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);
}

#[test]
#[allow(clippy::eq_op)]
fn test_projective_point_equality() {
    let a = G2Projective::generator();
    let b = G2Projective::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);

    let z = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xba7a_fa1f_9a6f_e250,
            0xfa0f_5b59_5eaf_e731,
            0x3bdc_4776_94c3_06e7,
            0x2149_be4b_3949_fa24,
            0x64aa_6e06_49b2_078c,
            0x12b1_08ac_3364_3c3e,
        ]),
        c1: Fp::from_raw_unchecked([
            0x1253_25df_3d35_b5a8,
            0xdc46_9ef5_555d_7fe3,
            0x02d7_16d2_4431_06a9,
            0x05a1_db59_a6ff_37d0,
            0x7cf7_784e_5300_bb8f,
            0x16a8_8922_c7a5_e844,
        ]),
    };

    let mut c = G2Projective {
        x: a.x * z,
        y: a.y * z,
        z,
    };
    assert!(bool::from(c.is_on_curve()));

    assert!(a == c);
    assert!(b != c);
    assert!(c == a);
    assert!(c != b);

    c.y = -c.y;
    assert!(bool::from(c.is_on_curve()));

    assert!(a != c);
    assert!(b != c);
    assert!(c != a);
    assert!(c != b);

    c.y = -c.y;
    c.x = z;
    assert!(!bool::from(c.is_on_curve()));
    assert!(a != b);
    assert!(a != c);
    assert!(b != c);
}

#[test]
fn test_conditionally_select_affine() {
    let a = G2Affine::generator();
    let b = G2Affine::identity();

    assert_eq!(G2Affine::conditional_select(&a, &b, Choice::from(0u8)), a);
    assert_eq!(G2Affine::conditional_select(&a, &b, Choice::from(1u8)), b);
}

#[test]
fn test_conditionally_select_projective() {
    let a = G2Projective::generator();
    let b = G2Projective::identity();

    assert_eq!(
        G2Projective::conditional_select(&a, &b, Choice::from(0u8)),
        a
    );
    assert_eq!(
        G2Projective::conditional_select(&a, &b, Choice::from(1u8)),
        b
    );
}

#[test]
fn test_projective_to_affine() {
    let a = G2Projective::generator();
    let b = G2Projective::identity();

    assert!(bool::from(G2Affine::from(a).is_on_curve()));
    assert!(!bool::from(G2Affine::from(a).is_identity()));
    assert!(bool::from(G2Affine::from(b).is_on_curve()));
    assert!(bool::from(G2Affine::from(b).is_identity()));

    let z = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xba7a_fa1f_9a6f_e250,
            0xfa0f_5b59_5eaf_e731,
            0x3bdc_4776_94c3_06e7,
            0x2149_be4b_3949_fa24,
            0x64aa_6e06_49b2_078c,
            0x12b1_08ac_3364_3c3e,
        ]),
        c1: Fp::from_raw_unchecked([
            0x1253_25df_3d35_b5a8,
            0xdc46_9ef5_555d_7fe3,
            0x02d7_16d2_4431_06a9,
            0x05a1_db59_a6ff_37d0,
            0x7cf7_784e_5300_bb8f,
            0x16a8_8922_c7a5_e844,
        ]),
    };

    let c = G2Projective {
        x: a.x * z,
        y: a.y * z,
        z,
    };

    assert_eq!(G2Affine::from(c), G2Affine::generator());
}

#[test]
fn test_affine_to_projective() {
    let a = G2Affine::generator();
    let b = G2Affine::identity();

    assert!(bool::from(G2Projective::from(a).is_on_curve()));
    assert!(!bool::from(G2Projective::from(a).is_identity()));
    assert!(bool::from(G2Projective::from(b).is_on_curve()));
    assert!(bool::from(G2Projective::from(b).is_identity()));
}

#[test]
fn test_doubling() {
    {
        let tmp = G2Projective::identity().double();
        assert!(bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));
    }
    {
        let tmp = G2Projective::generator().double();
        assert!(!bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));

        assert_eq!(
            G2Affine::from(tmp),
            G2Affine {
                x: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0xe9d9_e2da_9620_f98b,
                        0x54f1_1993_46b9_7f36,
                        0x3db3_b820_376b_ed27,
                        0xcfdb_31c9_b0b6_4f4c,
                        0x41d7_c127_8635_4493,
                        0x0571_0794_c255_c064,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xd6c1_d3ca_6ea0_d06e,
                        0xda0c_bd90_5595_489f,
                        0x4f53_52d4_3479_221d,
                        0x8ade_5d73_6f8c_97e0,
                        0x48cc_8433_925e_f70e,
                        0x08d7_ea71_ea91_ef81,
                    ]),
                },
                y: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x15ba_26eb_4b0d_186f,
                        0x0d08_6d64_b7e9_e01e,
                        0xc8b8_48dd_652f_4c78,
                        0xeecf_46a6_123b_ae4f,
                        0x255e_8dd8_b6dc_812a,
                        0x1641_42af_21dc_f93f,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0xf9b4_a1a8_9598_4db4,
                        0xd417_b114_cccf_f748,
                        0x6856_301f_c89f_086e,
                        0x41c7_7787_8931_e3da,
                        0x3556_b155_066a_2105,
                        0x00ac_f7d3_25cb_89cf,
                    ]),
                },
                infinity: Choice::from(0u8)
            }
        );
    }
}

#[test]
fn test_projective_addition() {
    {
        let a = G2Projective::identity();
        let b = G2Projective::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
    {
        let a = G2Projective::identity();
        let mut b = G2Projective::generator();
        {
            let z = Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xba7a_fa1f_9a6f_e250,
                    0xfa0f_5b59_5eaf_e731,
                    0x3bdc_4776_94c3_06e7,
                    0x2149_be4b_3949_fa24,
                    0x64aa_6e06_49b2_078c,
                    0x12b1_08ac_3364_3c3e,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x1253_25df_3d35_b5a8,
                    0xdc46_9ef5_555d_7fe3,
                    0x02d7_16d2_4431_06a9,
                    0x05a1_db59_a6ff_37d0,
                    0x7cf7_784e_5300_bb8f,
                    0x16a8_8922_c7a5_e844,
                ]),
            };

            b = G2Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = a + b;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G2Projective::generator());
    }
    {
        let a = G2Projective::identity();
        let mut b = G2Projective::generator();
        {
            let z = Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xba7a_fa1f_9a6f_e250,
                    0xfa0f_5b59_5eaf_e731,
                    0x3bdc_4776_94c3_06e7,
                    0x2149_be4b_3949_fa24,
                    0x64aa_6e06_49b2_078c,
                    0x12b1_08ac_3364_3c3e,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x1253_25df_3d35_b5a8,
                    0xdc46_9ef5_555d_7fe3,
                    0x02d7_16d2_4431_06a9,
                    0x05a1_db59_a6ff_37d0,
                    0x7cf7_784e_5300_bb8f,
                    0x16a8_8922_c7a5_e844,
                ]),
            };

            b = G2Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = b + a;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G2Projective::generator());
    }
    {
        let a = G2Projective::generator().double().double(); // 4P
        let b = G2Projective::generator().double(); // 2P
        let c = a + b;

        let mut d = G2Projective::generator();
        for _ in 0..5 {
            d += G2Projective::generator();
        }
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(!bool::from(d.is_identity()));
        assert!(bool::from(d.is_on_curve()));
        assert_eq!(c, d);
    }

    // Degenerate case
    {
        let beta = Fp2 {
            c0: Fp::from_raw_unchecked([
                0xcd03_c9e4_8671_f071,
                0x5dab_2246_1fcd_a5d2,
                0x5870_42af_d385_1b95,
                0x8eb6_0ebe_01ba_cb9e,
                0x03f9_7d6e_83d0_50d2,
                0x18f0_2065_5463_8741,
            ]),
            c1: Fp::zero(),
        };
        let beta = beta.square();
        let a = G2Projective::generator().double().double();
        let b = G2Projective {
            x: a.x * beta,
            y: -a.y,
            z: a.z,
        };
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(b.is_on_curve()));

        let c = a + b;
        assert_eq!(
            G2Affine::from(c),
            G2Affine::from(G2Projective {
                x: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x705a_bc79_9ca7_73d3,
                        0xfe13_2292_c1d4_bf08,
                        0xf37e_ce3e_07b2_b466,
                        0x887e_1c43_f447_e301,
                        0x1e09_70d0_33bc_77e8,
                        0x1985_c81e_20a6_93f2,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x1d79_b25d_b36a_b924,
                        0x2394_8e4d_5296_39d3,
                        0x471b_a7fb_0d00_6297,
                        0x2c36_d4b4_465d_c4c0,
                        0x82bb_c3cf_ec67_f538,
                        0x051d_2728_b67b_f952,
                    ])
                },
                y: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x41b1_bbf6_576c_0abf,
                        0xb6cc_9371_3f7a_0f9a,
                        0x6b65_b43e_48f3_f01f,
                        0xfb7a_4cfc_af81_be4f,
                        0x3e32_dadc_6ec2_2cb6,
                        0x0bb0_fc49_d798_07e3,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x7d13_9778_8f5f_2ddf,
                        0xab29_0714_4ff0_d8e8,
                        0x5b75_73e0_cdb9_1f92,
                        0x4cb8_932d_d31d_af28,
                        0x62bb_fac6_db05_2a54,
                        0x11f9_5c16_d14c_3bbe,
                    ])
                },
                z: Fp2::one()
            })
        );
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
}

#[test]
fn test_mixed_addition() {
    {
        let a = G2Affine::identity();
        let b = G2Projective::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
    {
        let a = G2Affine::identity();
        let mut b = G2Projective::generator();
        {
            let z = Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xba7a_fa1f_9a6f_e250,
                    0xfa0f_5b59_5eaf_e731,
                    0x3bdc_4776_94c3_06e7,
                    0x2149_be4b_3949_fa24,
                    0x64aa_6e06_49b2_078c,
                    0x12b1_08ac_3364_3c3e,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x1253_25df_3d35_b5a8,
                    0xdc46_9ef5_555d_7fe3,
                    0x02d7_16d2_4431_06a9,
                    0x05a1_db59_a6ff_37d0,
                    0x7cf7_784e_5300_bb8f,
                    0x16a8_8922_c7a5_e844,
                ]),
            };

            b = G2Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = a + b;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G2Projective::generator());
    }
    {
        let a = G2Affine::identity();
        let mut b = G2Projective::generator();
        {
            let z = Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xba7a_fa1f_9a6f_e250,
                    0xfa0f_5b59_5eaf_e731,
                    0x3bdc_4776_94c3_06e7,
                    0x2149_be4b_3949_fa24,
                    0x64aa_6e06_49b2_078c,
                    0x12b1_08ac_3364_3c3e,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x1253_25df_3d35_b5a8,
                    0xdc46_9ef5_555d_7fe3,
                    0x02d7_16d2_4431_06a9,
                    0x05a1_db59_a6ff_37d0,
                    0x7cf7_784e_5300_bb8f,
                    0x16a8_8922_c7a5_e844,
                ]),
            };

            b = G2Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = b + a;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G2Projective::generator());
    }
    {
        let a = G2Projective::generator().double().double(); // 4P
        let b = G2Projective::generator().double(); // 2P
        let c = a + b;

        let mut d = G2Projective::generator();
        for _ in 0..5 {
            d += G2Affine::generator();
        }
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(!bool::from(d.is_identity()));
        assert!(bool::from(d.is_on_curve()));
        assert_eq!(c, d);
    }

    // Degenerate case
    {
        let beta = Fp2 {
            c0: Fp::from_raw_unchecked([
                0xcd03_c9e4_8671_f071,
                0x5dab_2246_1fcd_a5d2,
                0x5870_42af_d385_1b95,
                0x8eb6_0ebe_01ba_cb9e,
                0x03f9_7d6e_83d0_50d2,
                0x18f0_2065_5463_8741,
            ]),
            c1: Fp::zero(),
        };
        let beta = beta.square();
        let a = G2Projective::generator().double().double();
        let b = G2Projective {
            x: a.x * beta,
            y: -a.y,
            z: a.z,
        };
        let a = G2Affine::from(a);
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(b.is_on_curve()));

        let c = a + b;
        assert_eq!(
            G2Affine::from(c),
            G2Affine::from(G2Projective {
                x: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x705a_bc79_9ca7_73d3,
                        0xfe13_2292_c1d4_bf08,
                        0xf37e_ce3e_07b2_b466,
                        0x887e_1c43_f447_e301,
                        0x1e09_70d0_33bc_77e8,
                        0x1985_c81e_20a6_93f2,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x1d79_b25d_b36a_b924,
                        0x2394_8e4d_5296_39d3,
                        0x471b_a7fb_0d00_6297,
                        0x2c36_d4b4_465d_c4c0,
                        0x82bb_c3cf_ec67_f538,
                        0x051d_2728_b67b_f952,
                    ])
                },
                y: Fp2 {
                    c0: Fp::from_raw_unchecked([
                        0x41b1_bbf6_576c_0abf,
                        0xb6cc_9371_3f7a_0f9a,
                        0x6b65_b43e_48f3_f01f,
                        0xfb7a_4cfc_af81_be4f,
                        0x3e32_dadc_6ec2_2cb6,
                        0x0bb0_fc49_d798_07e3,
                    ]),
                    c1: Fp::from_raw_unchecked([
                        0x7d13_9778_8f5f_2ddf,
                        0xab29_0714_4ff0_d8e8,
                        0x5b75_73e0_cdb9_1f92,
                        0x4cb8_932d_d31d_af28,
                        0x62bb_fac6_db05_2a54,
                        0x11f9_5c16_d14c_3bbe,
                    ])
                },
                z: Fp2::one()
            })
        );
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
}

#[test]
#[allow(clippy::eq_op)]
fn test_projective_negation_and_subtraction() {
    let a = G2Projective::generator().double();
    assert_eq!(a + (-a), G2Projective::identity());
    assert_eq!(a + (-a), a - a);
}

#[test]
fn test_affine_negation_and_subtraction() {
    let a = G2Affine::generator();
    assert_eq!(G2Projective::from(a) + (-a), G2Projective::identity());
    assert_eq!(G2Projective::from(a) + (-a), G2Projective::from(a) - a);
}

#[test]
fn test_projective_scalar_multiplication() {
    let g = G2Projective::generator();
    let a = Scalar::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Scalar::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!((g * a) * b, g * c);
}

#[test]
fn test_affine_scalar_multiplication() {
    let g = G2Affine::generator();
    let a = Scalar::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Scalar::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!(G2Affine::from(g * a) * b, g * c);
}

#[test]
fn test_is_torsion_free() {
    let a = G2Affine {
        x: Fp2 {
            c0: Fp::from_raw_unchecked([
                0x89f5_50c8_13db_6431,
                0xa50b_e8c4_56cd_8a1a,
                0xa45b_3741_14ca_e851,
                0xbb61_90f5_bf7f_ff63,
                0x970c_a02c_3ba8_0bc7,
                0x02b8_5d24_e840_fbac,
            ]),
            c1: Fp::from_raw_unchecked([
                0x6888_bc53_d707_16dc,
                0x3dea_6b41_1768_2d70,
                0xd8f5_f930_500c_a354,
                0x6b5e_cb65_56f5_c155,
                0xc96b_ef04_3477_8ab0,
                0x0508_1505_5150_06ad,
            ]),
        },
        y: Fp2 {
            c0: Fp::from_raw_unchecked([
                0x3cf1_ea0d_434b_0f40,
                0x1a0d_c610_e603_e333,
                0x7f89_9561_60c7_2fa0,
                0x25ee_03de_cf64_31c5,
                0xeee8_e206_ec0f_e137,
                0x0975_92b2_26df_ef28,
            ]),
            c1: Fp::from_raw_unchecked([
                0x71e8_bb5f_2924_7367,
                0xa5fe_049e_2118_31ce,
                0x0ce6_b354_502a_3896,
                0x93b0_1200_0997_314e,
                0x6759_f3b6_aa5b_42ac,
                0x1569_44c4_dfe9_2bbb,
            ]),
        },
        infinity: Choice::from(0u8),
    };
    assert!(!bool::from(a.is_torsion_free()));

    assert!(bool::from(G2Affine::identity().is_torsion_free()));
    assert!(bool::from(G2Affine::generator().is_torsion_free()));
}

#[cfg(feature = "endo")]
#[test]
fn test_mul_by_x() {
    // multiplying by `x` a point in G2 is the same as multiplying by
    // the equivalent scalar.
    let generator = G2Projective::generator();
    let x = if crate::BLS_X_IS_NEGATIVE {
        -Scalar::from(crate::BLS_X)
    } else {
        Scalar::from(crate::BLS_X)
    };
    assert_eq!(generator.mul_by_x(), generator * x);

    let point = G2Projective::generator() * Scalar::from(42);
    assert_eq!(point.mul_by_x(), point * x);
}

#[cfg(feature = "endo")]
#[test]
fn test_psi() {
    let generator = G2Projective::generator();

    let z = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x0ef2ddffab187c0a,
            0x2424522b7d5ecbfc,
            0xc6f341a3398054f4,
            0x5523ddf409502df0,
            0xd55c0b5a88e0dd97,
            0x066428d704923e52,
        ]),
        c1: Fp::from_raw_unchecked([
            0x538bbe0c95b4878d,
            0xad04a50379522881,
            0x6d5c05bf5c12fb64,
            0x4ce4a069a2d34787,
            0x59ea6c8d0dffaeaf,
            0x0d42a083a75bd6f3,
        ]),
    };

    // `point` is a random point in the curve
    let point = G2Projective {
        x: Fp2 {
            c0: Fp::from_raw_unchecked([
                0xee4c8cb7c047eaf2,
                0x44ca22eee036b604,
                0x33b3affb2aefe101,
                0x15d3e45bbafaeb02,
                0x7bfc2154cd7419a4,
                0x0a2d0c2b756e5edc,
            ]),
            c1: Fp::from_raw_unchecked([
                0xfc224361029a8777,
                0x4cbf2baab8740924,
                0xc5008c6ec6592c89,
                0xecc2c57b472a9c2d,
                0x8613eafd9d81ffb1,
                0x10fe54daa2d3d495,
            ]),
        } * z,
        y: Fp2 {
            c0: Fp::from_raw_unchecked([
                0x7de7edc43953b75c,
                0x58be1d2de35e87dc,
                0x5731d30b0e337b40,
                0xbe93b60cfeaae4c9,
                0x8b22c203764bedca,
                0x01616c8d1033b771,
            ]),
            c1: Fp::from_raw_unchecked([
                0xea126fe476b5733b,
                0x85cee68b5dae1652,
                0x98247779f7272b04,
                0xa649c8b468c6e808,
                0xb5b9a62dff0c4e45,
                0x1555b67fc7bbe73d,
            ]),
        },
        z: z.square() * z,
    };
    assert!(bool::from(point.is_on_curve()));

    // psi2(P) = psi(psi(P))
    assert_eq!(generator.psi2(), generator.psi().psi());
    assert_eq!(point.psi2(), point.psi().psi());
    // psi(P) is a morphism
    assert_eq!(generator.double().psi(), generator.psi().double());
    assert_eq!(point.psi() + generator.psi(), (point + generator).psi());
    // psi(P) behaves in the same way on the same projective point
    let mut normalized_point = [G2Affine::identity()];
    G2Projective::batch_normalize(&[point], &mut normalized_point);
    let normalized_point = G2Projective::from(normalized_point[0]);
    assert_eq!(point.psi(), normalized_point.psi());
    assert_eq!(point.psi2(), normalized_point.psi2());
}

#[test]
fn test_clear_cofactor() {
    let z = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x0ef2ddffab187c0a,
            0x2424522b7d5ecbfc,
            0xc6f341a3398054f4,
            0x5523ddf409502df0,
            0xd55c0b5a88e0dd97,
            0x066428d704923e52,
        ]),
        c1: Fp::from_raw_unchecked([
            0x538bbe0c95b4878d,
            0xad04a50379522881,
            0x6d5c05bf5c12fb64,
            0x4ce4a069a2d34787,
            0x59ea6c8d0dffaeaf,
            0x0d42a083a75bd6f3,
        ]),
    };

    // `point` is a random point in the curve
    let point = G2Projective {
        x: Fp2 {
            c0: Fp::from_raw_unchecked([
                0xee4c8cb7c047eaf2,
                0x44ca22eee036b604,
                0x33b3affb2aefe101,
                0x15d3e45bbafaeb02,
                0x7bfc2154cd7419a4,
                0x0a2d0c2b756e5edc,
            ]),
            c1: Fp::from_raw_unchecked([
                0xfc224361029a8777,
                0x4cbf2baab8740924,
                0xc5008c6ec6592c89,
                0xecc2c57b472a9c2d,
                0x8613eafd9d81ffb1,
                0x10fe54daa2d3d495,
            ]),
        } * z,
        y: Fp2 {
            c0: Fp::from_raw_unchecked([
                0x7de7edc43953b75c,
                0x58be1d2de35e87dc,
                0x5731d30b0e337b40,
                0xbe93b60cfeaae4c9,
                0x8b22c203764bedca,
                0x01616c8d1033b771,
            ]),
            c1: Fp::from_raw_unchecked([
                0xea126fe476b5733b,
                0x85cee68b5dae1652,
                0x98247779f7272b04,
                0xa649c8b468c6e808,
                0xb5b9a62dff0c4e45,
                0x1555b67fc7bbe73d,
            ]),
        },
        z: z.square() * z,
    };

    assert!(bool::from(point.is_on_curve()));
    assert!(!bool::from(G2Affine::from(point).is_torsion_free()));
    let cleared_point = point.clear_cofactor();

    assert!(bool::from(cleared_point.is_on_curve()));
    assert!(bool::from(G2Affine::from(cleared_point).is_torsion_free()));

    // the generator (and the identity) are always on the curve,
    // even after clearing the cofactor
    let generator = G2Projective::generator();
    assert!(bool::from(generator.clear_cofactor().is_on_curve()));
    let id = G2Projective::identity();
    assert!(bool::from(id.clear_cofactor().is_on_curve()));

    // test the effect on q-torsion points multiplying by h_eff modulo |Scalar|
    // h_eff % q = 0x2b116900400069009a40200040001ffff
    let h_eff_modq: [u8; 32] = [
        0xff, 0xff, 0x01, 0x00, 0x04, 0x00, 0x02, 0xa4, 0x09, 0x90, 0x06, 0x00, 0x04, 0x90, 0x16,
        0xb1, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00,
    ];
    assert_eq!(generator.clear_cofactor(), generator.multiply(&h_eff_modq));
    assert_eq!(
        cleared_point.clear_cofactor(),
        cleared_point.multiply(&h_eff_modq)
    );
}

#[test]
fn test_batch_normalize() {
    let a = G2Projective::generator().double();
    let b = a.double();
    let c = b.double();

    for a_identity in (0..1).map(|n| n == 1) {
        for b_identity in (0..1).map(|n| n == 1) {
            for c_identity in (0..1).map(|n| n == 1) {
                let mut v = [a, b, c];
                if a_identity {
                    v[0] = G2Projective::identity()
                }
                if b_identity {
                    v[1] = G2Projective::identity()
                }
                if c_identity {
                    v[2] = G2Projective::identity()
                }

                let mut t = [
                    G2Affine::identity(),
                    G2Affine::identity(),
                    G2Affine::identity(),
                ];
                let expected = [
                    G2Affine::from(v[0]),
                    G2Affine::from(v[1]),
                    G2Affine::from(v[2]),
                ];

                G2Projective::batch_normalize(&v[..], &mut t[..]);

                assert_eq!(&t[..], &expected[..]);
            }
        }
    }
}

#[cfg(feature = "hashing")]
#[test]
fn test_hash() {
    use crate::hash_to_field::ExpandMsgXmd;
    use std::convert::TryFrom;
    const DST: &'static [u8] = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";

    let tests: [(&[u8], &str); 1] = [
        // (b"", "05cb8437535e20ecffaef7752baddf98034139c38452458baeefab379ba13dff5bf5dd71b72418717047f5b0f37da03d0141ebfbdca40eb85b87142e130ab689c673cf60f1a3e98d69335266f30d9b8d4ac44c1038e9dcdd5393faf5c41fb78a12424ac32561493f3fe3c260708a12b7c620e7be00099a974e259ddc7d1f6395c3c811cdd19f1e8dbf3e9ecfdcbab8d60503921d7f6a12805e72940b963c0cf3471c7b2a524950ca195d11062ee75ec076daf2d4bc358c4b190c0c98064fdd92"),
        (b"abc", "139cddbccdc5e91b9623efd38c49f81a6f83f175e80b06fc374de9eb4b41dfe4ca3a230ed250fbe3a2acf73a41177fd802c2d18e033b960562aae3cab37a27ce00d80ccd5ba4b7fe0e7a210245129dbec7780ccc7954725f4168aff2787776e600aa65dae3c8d732d10ecd2c50f8a1baf3001578f71c694e03866e9f3d49ac1e1ce70dd94a733534f106d4cec0eddd161787327b68159716a37440985269cf584bcb1e621d3a7202be6ea05c4cfe244aeb197642555a0645fb87bf7466b2ba48"),
        // (b"abcdef0123456789", "190d119345b94fbd15497bcba94ecf7db2cbfd1e1fe7da034d26cbba169fb3968288b3fafb265f9ebd380512a71c3f2c121982811d2491fde9ba7ed31ef9ca474f0e1501297f68c298e9f4c0028add35aea8bb83d53c08cfc007c1e005723cd00bb5e7572275c567462d91807de765611490205a941a5a6af3b1691bfe596c31225d3aabdf15faff860cb4ef17c7c3be05571a0f8d3c08d094576981f4a3b8eda0a8e771fcdcc8ecceaf1356a6acf17574518acb506e435b639353c2e14827c8"),
        // (b"q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", "0934aba516a52d8ae479939a91998299c76d39cc0c035cd18813bec433f587e2d7a4fef038260eef0cef4d02aae3eb9119a84dd7248a1066f737cc34502ee5555bd3c19f2ecdb3c7d9e24dc65d4e25e50d83f0f77105e955d78f4762d33c17da09bcccfa036b4847c9950780733633f13619994394c23ff0b32fa6b795844f4a0673e20282d07bc69641cee04f5e566214f81cd421617428bc3b9fe25afbb751d934a00493524bc4e065635b0555084dd54679df1536101b2c979c0152d09192"),
        // (b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "11fca2ff525572795a801eed17eb12785887c7b63fb77a42be46ce4a34131d71f7a73e95fee3f812aea3de78b4d0156901a6ba2f9a11fa5598b2d8ace0fbe0a0eacb65deceb476fbbcb64fd24557c2f4b18ecfc5663e54ae16a84f5ab7f6253403a47f8e6d1763ba0cad63d6114c0accbef65707825a511b251a660a9b3994249ae4e63fac38b23da0c398689ee2ab520b6798718c8aed24bc19cb27f866f1c9effcdbf92397ad6448b5c9db90d2b9da6cbabf48adc1adf59a1a28344e79d57e")
    ];

    for (msg, exp) in &tests {
        let a = G2Projective::hash::<ExpandMsgXmd<sha2::Sha256>>(msg, DST);
        let d = <[u8; 192]>::try_from(hex::decode(exp).unwrap().as_slice()).unwrap();
        let e = G2Affine::from_uncompressed(&d).unwrap();
        assert_eq!(a.to_affine(), e);
    }
}
