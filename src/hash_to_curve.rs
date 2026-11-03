//! Shared helpers for hashing to the BLS12-381 fields and groups (RFC 9380).

use core::num::NonZero;
use core::ops::{Add, Mul};
use elliptic_curve::consts::U16;
use hash2curve::{ExpandMsg, Expander};
use subtle::Choice;

/// The target security level of BLS12-381 in bytes (128 bits).
pub(crate) type SecurityLevel = U16;

/// The `sgn0` function from RFC 9380 section 4.1.
pub(crate) trait Sgn0 {
    /// Returns the parity of the field element: 1 is negative, 0 is non-negative.
    fn sgn0(&self) -> Choice;
}

/// Expands `msg` and `dst` into `N` uniformly random bytes with `X`.
///
/// # Panics
///
/// Panics if `dst` is empty or `N` is zero or larger than `X` can produce.
pub(crate) fn expand_message<X, const N: usize>(msg: &[u8], dst: &[u8]) -> [u8; N]
where
    X: ExpandMsg<SecurityLevel>,
{
    let dst = [dst];
    let mut okm = [0u8; N];
    u16::try_from(N)
        .ok()
        .and_then(NonZero::new)
        .and_then(|len| X::expand_message(&[msg], &dst, len).ok())
        .and_then(|mut expander| expander.fill_bytes(&mut okm).ok())
        .expect("expand_message requires a non-empty DST and a valid output length");
    okm
}

/// Evaluates the polynomial with coefficients `k` at the powers of `x` in `xs`.
pub(crate) fn evaluate_iso<F>(xs: &[F], k: &[F], zero: F) -> F
where
    F: Copy + Add<Output = F> + Mul<Output = F>,
{
    xs.iter().zip(k).fold(zero, |acc, (x, k)| acc + *x * *k)
}
