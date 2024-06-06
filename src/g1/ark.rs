use crate::{G1Projective, Scalar};
use ark_ec::{AffineRepr, Group as ArkGroup};
use ark_ff::UniformRand;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, Read, SerializationError, Valid, Validate,
    Write,
};
use group::cofactor::CofactorGroup;
use group::Group;
use num_traits::Zero;
use rand::Rng;

impl CanonicalSerialize for G1Projective {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        match compress {
            Compress::No => {
                let _ = writer
                    .write(&self.to_uncompressed())
                    .map_err(|e| SerializationError::IoError(e))?;
            }
            Compress::Yes => {
                let _ = writer
                    .write(&self.to_compressed())
                    .map_err(|e| SerializationError::IoError(e))?;
            }
        }
        Ok(())
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        match compress {
            Compress::No => G1Projective::UNCOMPRESSED_BYTES,
            Compress::Yes => G1Projective::COMPRESSED_BYTES,
        }
    }
}

impl CanonicalDeserialize for G1Projective {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        // We only allow valid points to be deserialized
        // ignore validate
        let g1 = match compress {
            Compress::No => {
                let mut bytes = [0u8; G1Projective::UNCOMPRESSED_BYTES];
                reader
                    .read(&mut bytes)
                    .map_err(|e| SerializationError::IoError(e))?;
                G1Projective::from_uncompressed(&bytes)
            }
            Compress::Yes => {
                let mut bytes = [0u8; G1Projective::COMPRESSED_BYTES];
                reader
                    .read(&mut bytes)
                    .map_err(|e| SerializationError::IoError(e))?;
                G1Projective::from_compressed(&bytes)
            }
        };
        match validate {
            Validate::Yes => {
                Option::<G1Projective>::from(g1).ok_or(SerializationError::InvalidData)
            }
            Validate::No => Ok(g1.unwrap()),
        }
    }
}

impl Valid for G1Projective {
    fn check(&self) -> Result<(), SerializationError> {
        if bool::from(self.is_torsion_free() & self.is_on_curve()) {
            Ok(())
        } else {
            Err(SerializationError::InvalidData)
        }
    }
}

impl UniformRand for G1Projective {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        G1Projective::random(rng)
    }
}

impl Zero for G1Projective {
    fn zero() -> Self {
        G1Projective::IDENTITY
    }

    fn is_zero(&self) -> bool {
        self.is_identity().into()
    }
}

impl ArkGroup for G1Projective {
    type ScalarField = Scalar;

    fn generator() -> Self {
        G1Projective::GENERATOR
    }

    fn double_in_place(&mut self) -> &mut Self {
        *self = self.double();
        self
    }

    fn mul_bigint(&self, other: impl AsRef<[u64]>) -> Self {
        let mut bytes = [0u8; 32];
        for (i, b) in other.as_ref().iter().enumerate() {
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&b.to_le_bytes())
        }
        G1Projective::multiply(self, &bytes)
    }
}
