//! Benchmarks using non-wasm32 architectures only
#![cfg(not(target_arch = "wasm32"))]

#[macro_use]
extern crate criterion;

use bls12_381_plus::*;

use criterion::Criterion;
use elliptic_curve::hash2curve::ExpandMsgXmd;
use std::hint::black_box;

fn criterion_benchmark(c: &mut Criterion) {
    // Hash-to-curve
    {
        let msg = b"bench message";
        let g1_dst = b"BLS12381G1_XMD:SHA-256_SSWU_RO_BENCH";
        let g1_encode_dst = b"BLS12381G1_XMD:SHA-256_SSWU_NU_BENCH";
        let g2_dst = b"BLS12381G2_XMD:SHA-256_SSWU_RO_BENCH";
        let g2_encode_dst = b"BLS12381G2_XMD:SHA-256_SSWU_NU_BENCH";

        c.bench_function("G1Projective hash to curve", move |b| {
            b.iter(|| {
                G1Projective::hash::<ExpandMsgXmd<sha2::Sha256>>(black_box(msg), black_box(g1_dst))
            })
        });
        c.bench_function("G1Projective encode to curve", move |b| {
            b.iter(|| {
                G1Projective::encode::<ExpandMsgXmd<sha2::Sha256>>(
                    black_box(msg),
                    black_box(g1_encode_dst),
                )
            })
        });
        c.bench_function("G2Projective hash to curve", move |b| {
            b.iter(|| {
                G2Projective::hash::<ExpandMsgXmd<sha2::Sha256>>(black_box(msg), black_box(g2_dst))
            })
        });
        c.bench_function("G2Projective encode to curve", move |b| {
            b.iter(|| {
                G2Projective::encode::<ExpandMsgXmd<sha2::Sha256>>(
                    black_box(msg),
                    black_box(g2_encode_dst),
                )
            })
        });
    }

    // Pairings
    {
        let g = G1Affine::generator();
        let h = G2Affine::generator();
        c.bench_function("full pairing", move |b| {
            b.iter(|| pairing(black_box(&g), black_box(&h)))
        });
        c.bench_function("G2 preparation for pairing", move |b| {
            b.iter(|| G2Prepared::from(h))
        });
        let prep = G2Prepared::from(h);
        c.bench_function("miller loop for pairing", move |b| {
            b.iter(|| multi_miller_loop(&[(&g, &prep)]))
        });
        let prep = G2Prepared::from(h);
        let r = multi_miller_loop(&[(&g, &prep)]);
        c.bench_function("final exponentiation for pairing", move |b| {
            b.iter(|| r.final_exponentiation())
        });
    }
    // G1Affine
    {
        let name = "G1Affine";
        let a = G1Affine::generator();
        let s = Scalar::from(1u64);
        let compressed = [0u8; 48];
        let uncompressed = [0u8; 96];
        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} subgroup check", name), move |b| {
            b.iter(|| black_box(a).is_torsion_free())
        });
        c.bench_function(
            &format!("{} deserialize compressed point", name),
            move |b| b.iter(|| G1Affine::from_compressed(black_box(&compressed))),
        );
        c.bench_function(
            &format!("{} deserialize uncompressed point", name),
            move |b| b.iter(|| G1Affine::from_uncompressed(black_box(&uncompressed))),
        );
    }

    // G1Projective
    {
        let name = "G1Projective";
        let a = G1Projective::GENERATOR;
        let a_affine = G1Affine::generator();
        let s = Scalar::from(1u64);

        const N: usize = 10000;
        let v = vec![G1Projective::GENERATOR; N];
        let mut q = vec![G1Affine::identity(); N];

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} to affine", name), move |b| {
            b.iter(|| G1Affine::from(black_box(a)))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("{} mixed addition", name), move |b| {
            b.iter(|| black_box(a).add_mixed(&a_affine))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G1Projective::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }

    // G2Affine
    {
        let name = "G2Affine";
        let a = G2Affine::generator();
        let s = Scalar::from(1u64);
        let compressed = [0u8; 96];
        let uncompressed = [0u8; 192];
        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} subgroup check", name), move |b| {
            b.iter(|| black_box(a).is_torsion_free())
        });
        c.bench_function(
            &format!("{} deserialize compressed point", name),
            move |b| b.iter(|| G2Affine::from_compressed(black_box(&compressed))),
        );
        c.bench_function(
            &format!("{} deserialize uncompressed point", name),
            move |b| b.iter(|| G2Affine::from_uncompressed(black_box(&uncompressed))),
        );
    }

    // G2Projective
    {
        let name = "G2Projective";
        let a = G2Projective::GENERATOR;
        let a_affine = G2Affine::generator();
        let s = Scalar::from(1u64);

        const N: usize = 10000;
        let v = vec![G2Projective::GENERATOR; N];
        let mut q = vec![G2Affine::identity(); N];

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} to affine", name), move |b| {
            b.iter(|| G2Affine::from(black_box(a)))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("{} mixed addition", name), move |b| {
            b.iter(|| black_box(a).add_mixed(&a_affine))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G2Projective::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
