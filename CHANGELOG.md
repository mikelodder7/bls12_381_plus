# Changelog

All notable changes to `bls12_381_plus` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.10.0] - 2026-11-03

### Added

- `G1Affine::sum_of_products` and `G2Affine::sum_of_products` (with `alloc`), ported from [mikelodder7/blstrs#2] by @andrewwhitehead.

### Security

- Removed `crypto-bigint` `0.5` ([RUSTSEC-2026-0387]) by dropping the `elliptic-curve` `0.13` dependency, and updated `crypto-bigint` `0.7` to `0.7.10`.

### Breaking

- Removed the `elliptic_curve_013` re-export and all `elliptic-curve` `0.13` trait impls (`Curve`, `PrimeCurve`, `FieldBytesEncoding`, `AffineCoordinates`, `MulByGenerator`, `LinearCombination`, `ScalarPrimitive`, `FromUintUnchecked`, `Reduce`, `Invert`, `IsHigh`). The `ff`/`group` `0.13` traits stay because `pairing` needs them.
- `hash`/`encode` now take an `X: hash2curve::ExpandMsg<U16>` from `hash2curve` `0.14`, for example `hash2curve::ExpandMsgXmd<sha2::Sha256>` with `sha2` `0.11`.
- `Scalar` conversions to and from `U256`/`U384`/`U512` now use `crypto-bigint` `0.7` types. `From<Scalar> for U384`/`U512` now produce the correct value, and `From<U384>`/`From<[u8; 48]>` reduce modulo the group order instead of truncating.
- The `elliptic-curve` `0.14` `AffineCoordinates::FieldRepr` for `G1Affine`/`G2Affine` is now `elliptic_curve::array::Array` instead of `generic_array::GenericArray`.

### Changed

- Added `hash2curve` `0.14` for `expand_message`, and bumped `sha2` to `0.11`.
- G1 and G2 now share one isogeny evaluation helper. Output is unchanged and checked against the RFC 9380 test vectors.

[mikelodder7/blstrs#2]: https://github.com/mikelodder7/blstrs/pull/2
[RUSTSEC-2026-0387]: https://rustsec.org/advisories/RUSTSEC-2026-0387

