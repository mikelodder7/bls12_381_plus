//! # `bls12_381`
//!
//! This crate provides an implementation of the BLS12-381 pairing-friendly elliptic
//! curve construction.
//!
//! * **This implementation has not been reviewed or audited. Use at your own risk.**
//! * This implementation targets Rust `1.36` or later.
//! * This implementation does not require the Rust standard library.
//! * All operations are constant time unless explicitly noted.

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
// Catch documentation errors caused by code changes.
#![deny(rustdoc::broken_intra_doc_links)]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
// #![deny(unsafe_code)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::many_single_char_names)]
// This lint is described at
// https://rust-lang.github.io/rust-clippy/master/index.html#suspicious_arithmetic_impl
// In our library, some of the arithmetic involving extension fields will necessarily
// involve various binary operators, and so this lint is triggered unnecessarily.
#![allow(clippy::suspicious_arithmetic_impl)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(any(test, feature = "std"))]
#[macro_use]
extern crate std;

pub use elliptic_curve;
pub use elliptic_curve_013;
pub use ff;
pub use ff_013;
#[cfg(feature = "groups")]
pub use group;
#[cfg(feature = "groups")]
pub use group_013;

#[cfg(test)]
#[cfg(feature = "groups")]
mod tests;

#[macro_use]
mod util;

/// Notes about how the BLS12-381 elliptic curve is designed, specified
/// and implemented by this library.
pub mod notes {
    pub mod design;
    pub mod serialization;
}

mod compat_014;
mod scalar;

#[cfg(target_arch = "wasm32")]
pub use scalar::run_test_wasm;
pub use scalar::{Scalar, ScalarLe};

#[cfg(all(feature = "groups", not(feature = "expose-fields")))]
mod fp;
#[cfg(feature = "expose-fields")]
pub mod fp;
#[cfg(all(feature = "groups", not(feature = "expose-fields")))]
mod fp2;
#[cfg(feature = "expose-fields")]
pub mod fp2;
#[cfg(feature = "groups")]
mod g1;
#[cfg(feature = "groups")]
mod g2;

#[cfg(feature = "groups")]
pub use g1::{G1Affine, G1Projective};
#[cfg(feature = "expose-fields")]
pub use g1::{G1Compressed, G1Uncompressed};
#[cfg(feature = "groups")]
pub use g2::{G2Affine, G2Projective};
#[cfg(feature = "expose-fields")]
pub use g2::{G2Compressed, G2Uncompressed};

mod fp12;
mod fp6;

use elliptic_curve_013::{
    bigint::{ArrayEncoding, U384},
    consts::U48,
    point::PointCompression,
    Curve, FieldBytes, FieldBytesEncoding, PrimeCurve,
};

/// An engine for operations generic G1 operations
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd)]
pub struct Bls12381G1;

unsafe impl Send for Bls12381G1 {}
unsafe impl Sync for Bls12381G1 {}

/// An engine for operations generic G2 operations
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd)]
pub struct Bls12381G2;

impl Curve for Bls12381G1 {
    type FieldBytesSize = U48;
    type Uint = U384;
    const ORDER: U384 = U384::from_be_hex("0000000000000000000000000000000073eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
}

impl PrimeCurve for Bls12381G1 {}

impl PointCompression for Bls12381G1 {
    const COMPRESS_POINTS: bool = true;
}

impl FieldBytesEncoding<Bls12381G1> for U384 {
    fn decode_field_bytes(field_bytes: &FieldBytes<Bls12381G1>) -> Self {
        U384::from_be_byte_array(*field_bytes)
    }

    fn encode_field_bytes(&self) -> FieldBytes<Bls12381G1> {
        self.to_be_byte_array()
    }
}

impl Curve for Bls12381G2 {
    type FieldBytesSize = U48;
    type Uint = U384;
    const ORDER: U384 = U384::from_be_hex("0000000000000000000000000000000073eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
}

impl PrimeCurve for Bls12381G2 {}

impl PointCompression for Bls12381G2 {
    const COMPRESS_POINTS: bool = true;
}

impl FieldBytesEncoding<Bls12381G2> for U384 {
    fn decode_field_bytes(field_bytes: &FieldBytes<Bls12381G2>) -> Self {
        U384::from_be_byte_array(*field_bytes)
    }

    fn encode_field_bytes(&self) -> FieldBytes<Bls12381G2> {
        self.to_be_byte_array()
    }
}

/// The BLS parameter x for BLS12-381 is -0xd201000000010000
#[cfg(feature = "groups")]
const BLS_X: u64 = 0xd201_0000_0001_0000;
#[cfg(feature = "groups")]
const BLS_X_IS_NEGATIVE: bool = true;

#[cfg(feature = "pairings")]
mod pairings;

#[cfg(feature = "pairings")]
pub use pairings::{pairing, Bls12, Gt, MillerLoopResult};

#[cfg(feature = "pairings")]
pub use pairings::{multi_miller_loop, G2Prepared};

mod isogeny;
