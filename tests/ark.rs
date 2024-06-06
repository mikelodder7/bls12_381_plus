#[cfg(feature = "ark")]
mod ark_tests {
    use ark_poly::{
        univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain,
    };
    use bls12_381_plus::{G1Projective, Scalar};
    use ff::Field;

    #[test]
    fn fft() {
        let tau = Scalar::random(&mut rand::rngs::OsRng);
        let mut t_poly_coeffs = vec![tau; 5];
        t_poly_coeffs[1] = t_poly_coeffs[0] * tau;
        t_poly_coeffs[2] = t_poly_coeffs[1] * tau;
        t_poly_coeffs[3] = t_poly_coeffs[2] * tau;
        t_poly_coeffs[4] = t_poly_coeffs[3] * tau;
        let t_poly = DensePolynomial::<Scalar>::from_coefficients_slice(&t_poly_coeffs);
        let poly = DensePolynomial::<Scalar>::rand(5, &mut rand::thread_rng());

        let domain = GeneralEvaluationDomain::<Scalar>::new(8).unwrap();
        let aux_domain = GeneralEvaluationDomain::<Scalar>::new(16).unwrap();

        let t_evals = aux_domain.fft(t_poly.coeffs());
        let d_evals = aux_domain.fft(poly.coeffs());

        let dt_evals = t_evals
            .iter()
            .zip(d_evals.iter())
            .map(|(t, d)| G1Projective::GENERATOR * t * d)
            .collect::<Vec<_>>();

        let dt_poly = aux_domain.ifft(&dt_evals);
        domain.fft(&dt_poly[domain.size()..]);
    }
}
