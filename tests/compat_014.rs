#[cfg(feature = "pairings")]
use bls12_381_plus::Gt;
use bls12_381_plus::Scalar;
use bls12_381_plus::{elliptic_curve_014, ff_014, group_014};
use bls12_381_plus::{Bls12381G1, Bls12381G2, G1Affine, G1Projective, G2Affine, G2Projective};

fn assert_prime_field<F: ff_014::PrimeField>() {}

fn assert_prime_curve<C: elliptic_curve_014::PrimeCurve>() {}

fn assert_group<G: group_014::Group>() {}

fn assert_curve<C: group_014::Curve>() {}

fn assert_affine_coordinates<A: elliptic_curve_014::point::AffineCoordinates>() {}

#[cfg(feature = "alloc")]
fn assert_wnaf_group<G: group_014::WnafGroup>() {}

#[test]
fn scalar_implements_ff_014() {
    assert_prime_field::<Scalar>();

    let value = Scalar::from(42u64);
    let repr = <Scalar as ff_014::PrimeField>::to_repr(&value);
    let roundtrip = <Scalar as ff_014::PrimeField>::from_repr(repr).unwrap();
    assert_eq!(value, roundtrip);
}

#[cfg(feature = "bits")]
#[test]
fn scalar_implements_ff_014_bits() {
    fn assert_prime_field_bits<F: ff_014::PrimeFieldBits>() {}

    assert_prime_field_bits::<Scalar>();

    let value = Scalar::from(42u64);
    let bits = <Scalar as ff_014::PrimeFieldBits>::to_le_bits(&value);
    assert!(!bits[0]);
    assert!(bits[1]);
    assert!(!bits[2]);
    assert!(bits[3]);
    assert!(!bits[4]);
    assert!(bits[5]);
}

#[test]
fn curve_markers_implement_elliptic_curve_014() {
    assert_prime_curve::<Bls12381G1>();
    assert_prime_curve::<Bls12381G2>();
}

#[test]
fn groups_implement_group_014() {
    assert_group::<G1Projective>();
    assert_group::<G2Projective>();
    assert_curve::<G1Projective>();
    assert_curve::<G2Projective>();
    assert_affine_coordinates::<G1Affine>();
    assert_affine_coordinates::<G2Affine>();

    #[cfg(feature = "alloc")]
    {
        assert_wnaf_group::<G1Projective>();
        assert_wnaf_group::<G2Projective>();
        assert_eq!(
            <G1Projective as group_014::WnafGroup>::recommended_wnaf_for_num_scalars(1),
            4
        );
        assert_eq!(
            <G2Projective as group_014::WnafGroup>::recommended_wnaf_for_num_scalars(1),
            4
        );
    }

    let g1 = <G1Projective as group_014::Group>::generator();
    let g2 = <G2Projective as group_014::Group>::generator();
    assert!(!bool::from(
        <G1Projective as group_014::Group>::is_identity(&g1)
    ));
    assert!(!bool::from(
        <G2Projective as group_014::Group>::is_identity(&g2)
    ));
}

#[cfg(feature = "pairings")]
#[test]
fn gt_implements_group_014() {
    assert_group::<Gt>();

    let gt = <Gt as group_014::Group>::generator();
    assert!(!bool::from(<Gt as group_014::Group>::is_identity(&gt)));
}
