use crate::{G2Projective, Scalar};
use ark_ec::Group as ArkGroup;
use ark_ff::UniformRand;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, Read, SerializationError, Valid, Validate,
    Write,
};
use group::cofactor::CofactorGroup;
use group::Group;
use num_traits::Zero;
use rand::Rng;

impl CanonicalSerialize for G2Projective {
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
            Compress::No => G2Projective::UNCOMPRESSED_BYTES,
            Compress::Yes => G2Projective::COMPRESSED_BYTES,
        }
    }
}

impl CanonicalDeserialize for G2Projective {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        // We only allow valid points to be deserialized
        // ignore validate
        let g1 = match compress {
            Compress::No => {
                let mut bytes = [0u8; G2Projective::UNCOMPRESSED_BYTES];
                reader
                    .read(&mut bytes)
                    .map_err(|e| SerializationError::IoError(e))?;
                G2Projective::from_uncompressed(&bytes)
            }
            Compress::Yes => {
                let mut bytes = [0u8; G2Projective::COMPRESSED_BYTES];
                reader
                    .read(&mut bytes)
                    .map_err(|e| SerializationError::IoError(e))?;
                G2Projective::from_compressed(&bytes)
            }
        };
        match validate {
            Validate::Yes => {
                Option::<G2Projective>::from(g1).ok_or(SerializationError::InvalidData)
            }
            Validate::No => Ok(g1.unwrap()),
        }
    }
}

impl Valid for G2Projective {
    fn check(&self) -> Result<(), SerializationError> {
        if bool::from(self.is_torsion_free() & self.is_on_curve()) {
            Ok(())
        } else {
            Err(SerializationError::InvalidData)
        }
    }
}

impl UniformRand for G2Projective {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        G2Projective::random(rng)
    }
}

impl Zero for G2Projective {
    fn zero() -> Self {
        G2Projective::IDENTITY
    }

    fn is_zero(&self) -> bool {
        self.is_identity().into()
    }
}

impl ArkGroup for G2Projective {
    type ScalarField = Scalar;

    fn generator() -> Self {
        G2Projective::GENERATOR
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
        G2Projective::multiply(self, &bytes)
    }
}
