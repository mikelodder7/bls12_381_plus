# 0.9.0-rc0

- Release candidate for `0.9.0`.
- Pinned core crypto trait dependencies to exact latest compatible patch releases.
- Bumped `criterion` to `0.8.2`.
- Use `elliptic-curve-tools` for faster small/medium-size sum-of-products with `alloc`.
- Added aliased `elliptic-curve`, `ff`, and `group` `0.14` compatibility while preserving `0.13` traits for `pairing`.

# 0.8.11

- Added invert to Gt
- Added From<Vec<u8>> and TryFrom<Vec<u8>> for types

# 0.8.6

- Added debug_asserts to check lengths
- Fixed issue with 32-bit builds

# 0.8.3

- Add LowerHex and UpperHex impls for G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar
- Add from hex methods for G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar
- Removed heapless dependency

# 0.8.2

- Fix a bug in hash to g1

# 0.8.1

- Fix a serialization issue where is\_torsion\_free was not working for valid points

# 0.8.0

- Bumped dependencies to `ff 0.13`, `group 0.13`, `pairing 0.23`, `elliptic-curve 0.13`
- Change hash2curve to use elliptic-curve traits where possible.
- `Mul` implementations for `{G1Affine, G1Projective, G2Affine, G2Projective}`
    with `Scalar` on the left-hand side. These enable code to be written matching
    equations that use the additive notation `[a] P`.
- Zeroization support for `Gt`.

# 0.7.0

## Changed
- MSRV bumped to `1.56.0`
- Bumped dependencies to `ff 0.12`, `group 0.12`, `pairing 0.22`.

# 0.6.1

## Changed
- G2 arithmetic is now 25-30% faster across the board.
- Pairings are now 10-15% faster.

# 0.6.0

## Fixed
- `bls12_381::Gt::default()` now returns `Gt::identity()` instead of a nonsensical value.

## Added
- Zeroization support for most public types, behind the `zeroize` feature flag.
- `bls12_381::MillerLoopResult` trait implementations:
  - `Default`
  - `AddAssign<MillerLoopResult>`
  - `AddAssign<&MillerLoopResult>`

## Changed
- Bumped dependencies to `ff 0.11`, `group 0.11`, `pairing 0.21`.

## Removed
- The deprecated `endo` feature flag.

# 0.5.0

## Added
- `bits` feature flag (on by default) that exposes an `ff::PrimeFieldBits` implementation
  on `bls12_381::Scalar`.
- `experimental` feature flag, for features have no backwards-compatibility guarantees and
  may change at any time.
- Hashing to curves ([Internet Draft v11](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11)),
  behind the `experimental` feature flag.

## Changed
- Bumped dependencies to `ff 0.10`, `group 0.10`, `pairing 0.20`.
- MSRV is now 1.51.0.

# 0.4.0

## Changed
- Bumped dependencies to `ff 0.9`, `group 0.9`, `pairing 0.19`, `rand_core 0.6`.
- The `endo` feature is on by default. The `endo` feature flag itself is deprecated, and
  will be removed in a future minor release.
- MSRV is now 1.47.0.

# 0.3.1

# Fixed
* The crate now compiles for non-64-bit targets, such as the `wasm32-*` targets.

# 0.3.0

# Changed
* Migrated to `ff 0.8`, `group 0.8`, and `pairing 0.18`.
* MSRV is now 1.44.0.
* Switched to complete addition formulas for G1/G2.

# 0.2.0

This release adds implementations of the `ff`, `group`, and `pairing` traits (with the
latter two being gated by the `groups` and `pairings` feature flags respectively).
Additional trait implementations (for standard traits) have been added where the `ff`,
`group`, and `pairing` trait bounds require them.

## Added
* `bls12_381::Bls12`, a `pairing::Engine` for BLS12-381 pairing operations. It implements
  the following traits:
  * `pairing::{Engine, MultiMillerLoop}`
* New trait implementations for `bls12_381::G1Projective`:
  * `group::{Curve, Group, GroupEncoding, WnafGroup}`
  * `group::prime::{PrimeCurve, PrimeGroup}`
* New trait implementations for `bls12_381::G1Affine`:
  * `group::{GroupEncoding, UncompressedEncoding}`
  * `group::prime::PrimeCurveAffine`
  * `pairing::PairingCurveAffine`
* New trait implementations for `bls12_381::G2Projective`:
  * `group::{Curve, Group, GroupEncoding, WnafGroup}`
  * `group::prime::{PrimeCurve, PrimeGroup}`
* New trait implementations for `bls12_381::G2Affine`:
  * `group::{GroupEncoding, UncompressedEncoding}`
  * `group::prime::PrimeCurveAffine`
  * `pairing::PairingCurveAffine`
* New trait implementations for `bls12_381::Gt`:
  * `group::Group`
* New trait implementations for `bls12_381::MillerLoopResult`:
  * `pairing::MillerLoopResult`
* New trait implementations for `bls12_381::Scalar`:
  * `ff::{Field, PrimeField}`

# 0.1.1

Added `clear_cofactor` methods to `G1Projective` and `G2Projective`. If the crate feature `endo`
is enabled the G2 cofactor clearing will use the curve endomorphism technique described by
[Budroni-Pintore](https://ia.cr/2017/419). If the crate feature `endo` is _not_ enabled then
the code will simulate the effects of the Budroni-Pintore cofactor clearing in order to keep
the API consistent. In September 2020, when patents US7110538B2 and US7995752B2 expire, the
endo feature will be made default. However, for now it must be explicitly enabled.

# 0.1.0

Initial release.
