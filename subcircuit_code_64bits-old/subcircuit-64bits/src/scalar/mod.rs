mod ristretto255;

pub type Scalar = ristretto255::Scalar;
pub type ScalarBytes = curve25519_dalek::scalar::Scalar;

pub const S_MODULUS: [u8; 32] = [
  0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58, 0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10,
];

pub trait ScalarFromPrimitives {
  fn to_scalar(self) -> Scalar;
}

impl ScalarFromPrimitives for usize {
  #[inline]
  fn to_scalar(self) -> Scalar {
    (0..self).map(|_i| Scalar::one()).sum()
  }
}

impl ScalarFromPrimitives for bool {
  #[inline]
  fn to_scalar(self) -> Scalar {
    if self {
      Scalar::one()
    } else {
      Scalar::zero()
    }
  }
}

pub trait ScalarBytesFromScalar {
  fn decompress_scalar(s: &Scalar) -> ScalarBytes;
  fn decompress_vector(s: &[Scalar]) -> Vec<ScalarBytes>;
}

impl ScalarBytesFromScalar for Scalar {
  fn decompress_scalar(s: &Scalar) -> ScalarBytes {
    ScalarBytes::from_bytes_mod_order(s.to_bytes())
  }

  fn decompress_vector(s: &[Scalar]) -> Vec<ScalarBytes> {
    (0..s.len())
      .map(|i| Scalar::decompress_scalar(&s[i]))
      .collect::<Vec<ScalarBytes>>()
  }
}
