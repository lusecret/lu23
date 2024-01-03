use super::booleanity::{BooleanityProof, BooleanitySystem};
use super::groups::group::CompressedGroup;
use super::groups::transcript::{AppendToTranscript, ProofTranscript};
use super::groups::input_mapping::{InputMappingInternalProof};
use super::montgomery;
use super::polycommit::univariate_poly::UnivariatePolyCommit;
use super::scalar::Scalar;
use super::subcircuit::{SubCircuitProof, SubCircuitProofShort, SubCircuitSystem};
use super::subcircuit_errors::SubCircuitError;
use num_bigint::BigUint;
use merlin::Transcript;
use rand::rngs::OsRng;
use rand_core::{CryptoRng, RngCore};
use std::mem;

#[derive(Debug)]
pub struct BoolSubCircuitSystem {
  pub cir: SubCircuitSystem,
  pub bs: BooleanitySystem,
}

impl BoolSubCircuitSystem {
  pub fn new(label: &[u8], n: usize, bool_degrees: usize) -> Self {
    let cir = SubCircuitSystem::new(label, n);

    let label_static: &'static [u8] = unsafe { mem::transmute(label) };
    let bs = BooleanitySystem::new(bool_degrees, label_static);

    BoolSubCircuitSystem { cir, bs }
  }

  fn protocol_name() -> &'static [u8] {
    b"sub-protocol for circuit validatioin"
  }

  fn split_a(a: u32) -> Vec<u64> {
    let mut b: Vec<u64> = Vec::with_capacity(32);

    for i in 0..32 {
      let bit = (a >> i) & 1;
      b.push(bit as u64);
    }
    b
  }

  fn split_alpha<Rng: RngCore + CryptoRng>(alpha: u64, rng: &mut Rng) -> Vec<u64> {
    let betam: Vec<u64> = (0..31).map(|_| montgomery::random64_mod(rng)).collect();

    let betam_2pow: Vec<u128> = betam
      .iter()
      .enumerate()
      .map(|(i, &bm)| montgomery::mont_mul_mod(montgomery::POW_TWO[i], bm) as u128)
      .collect();

    let beta1_31m: u128 = betam_2pow.iter().sum();
    let beta1_31m: u64 = (beta1_31m % (montgomery::get_Q() as u128)) as u64;

    let mut beta_lastm: u64;
    let alpham: u64 = montgomery::enter_mont(alpha);
    if alpham > beta1_31m {
      beta_lastm = alpham - beta1_31m;
    } else {
      beta_lastm = alpham + montgomery::get_Q() - beta1_31m;
    }
    beta_lastm = montgomery::mont_mul_mod(beta_lastm, montgomery::INV_POW_OF_TWO);
    let beta_last: u64 = montgomery::back_from_mont(beta_lastm);

    let mut beta: Vec<u64> = betam
      .into_iter()
      .map(|bm| montgomery::back_from_mont(bm))
      .collect();
    beta.push(beta_last);

    beta
  }
}

#[derive(Debug, PartialEq)]
pub struct BoolSubCircuitProof {
  cir_proof: SubCircuitProof,
  bool_commit: UnivariatePolyCommit,
  bool_proof: BooleanityProof,
  b_prime: Vec<u64>,
}

#[derive(Clone, Debug, serde_derive::Deserialize, serde_derive::Serialize)]
pub struct BoolSubCircuitProofShort {
  cir_proof: SubCircuitProofShort,
  bool_commit: UnivariatePolyCommit,
  bool_proof: BooleanityProof,
}

impl BoolSubCircuitSystem {
  pub fn prove(
    &self,
    ao_vec: &Vec<u64>,
    alphao_vec: &Vec<u64>,
    v_vec: &Vec<Scalar>,
    phi: Scalar,
    n: usize,
    transcript: &mut Transcript,
  ) -> Result<BoolSubCircuitProof, SubCircuitError> {

    if !ao_vec.iter().all(|&a| a <= u32::MAX as u64) {
      return Err(SubCircuitError::IllegalParameters);
    }

    let (tau_vecs, mu_vecs, r, epsilon, M, R, Cs) = self
      .cir
      .prove_stage1(ao_vec, alphao_vec, v_vec, phi, n, transcript)?;

    let b_vec: Vec<u64> = ao_vec
      .iter()
      .flat_map(|&a| Self::split_a(a as u32))
      .collect();

    let mut csprng: OsRng = OsRng;
    let beta_vec: Vec<u64> = alphao_vec
      .iter()
      .flat_map(|&al| Self::split_alpha(al, &mut csprng))
      .collect();

    let (bool_setup, bool_blinds) = self.bs.setup(&b_vec, &beta_vec).unwrap();

    for C in bool_setup.comms.C.iter() {
      C.append_to_transcript(b"poly commitment in bool-subcircuit", transcript);
    }
    bool_setup
      .comms
      .U
      .append_to_transcript(b"poly commitment in bool-subcircuit", transcript);

    let x = transcript.challenge_mont(Self::protocol_name());

    let cir_proof = self.cir.prove_stage2(
      ao_vec, alphao_vec, v_vec, phi, tau_vecs, mu_vecs, x, r, epsilon, M, R, Cs, transcript,
    )?;

    let bool_proof = self
      .bs
      .prove(&bool_setup, bool_blinds, &b_vec, &beta_vec, x, transcript);
    let b_prime = montgomery::compute_aprime(&b_vec, &beta_vec, x);

    Ok(BoolSubCircuitProof {
      cir_proof,
      bool_commit: bool_setup.comms,
      bool_proof,
      b_prime,
    })
  }

  pub fn verify(
    &self,
    pedersen: &Vec<CompressedGroup>,
    n: usize,
    proof: &BoolSubCircuitProof,
    transcript: &mut Transcript,
  ) -> Result<(), SubCircuitError> {
    transcript.append_protocol_name(SubCircuitSystem::protocol_name());
    proof
      .cir_proof
      .R
      .append_to_transcript(b"point R", transcript);

    for C in proof.cir_proof.Cs.iter() {
      C.append_to_transcript(b"point C", transcript);
    }

    for C in proof.bool_commit.C.iter() {
      C.append_to_transcript(b"poly commitment in bool-subcircuit", transcript);
    }
    proof
      .bool_commit
      .U
      .append_to_transcript(b"poly commitment in bool-subcircuit", transcript);

    let x = transcript.challenge_mont(Self::protocol_name());

    self
      .cir
      .internal_verify(pedersen, n, x, &proof.cir_proof, transcript)?;

    self
      .bs
      .verify(&proof.bool_commit, &proof.bool_proof, x, transcript)?;

    let alen: usize = proof.cir_proof.ao_prime.len();
    if alen * 32 != proof.b_prime.len() {
      return Err(SubCircuitError::BooleanAPrimeBPrimeNotMatch);
    };

    let mut index = 0;
    while index + 32 <= proof.b_prime.len() {
      let chunk = &proof.b_prime[index..index + 32];
      let mut sum: u128 = 0;

      for (i, &num) in chunk.iter().enumerate() {
        sum += (num as u128) << (i as u32);
      }

      let sum: u64 = (sum % (montgomery::get_Q() as u128)) as u64;
      let a_exp = proof.cir_proof.ao_prime[index >> 5];
      if sum != a_exp {
        return Err(SubCircuitError::BooleanAPrimeBPrimeNotMatch);
      }

      index += 32;
    }

    Ok(())
  }

  pub fn simulate_write_read_file(&self, proof: BoolSubCircuitProof, chunksize: usize) -> BoolSubCircuitProof {
    let y_prime_length = (61 + 61 + chunksize.ilog2()) / 8 + 2; //y_prime in bytes; tau_i and x are 61-bits, y_prime is at most 61+61+chunksize bits
  //  println!("chunksize: {}", chunksize);
  //  println!("y_prime_length: {}", y_prime_length);

    let short_y_primes = proof.cir_proof
      .ys_prime
      .iter()
      .map(|x| x.to_bytes()[0..y_prime_length as usize].try_into().unwrap())
      .collect();
    let e_vec_bytes: Vec<Vec<u8>> = proof.cir_proof
      .mapping_proof
      .e_vec
      .iter()
      .map(|u| u.to_bytes_le())
      .collect();

    let SubCircuitProof_short: SubCircuitProofShort = SubCircuitProofShort {
      ys_prime: short_y_primes, // 17 * 1024 =
      ao_prime: proof.cir_proof.ao_prime.clone(), //30
      M: proof.cir_proof.M.clone(),
      R: proof.cir_proof.R.clone(),   //1
      Cs: proof.cir_proof.Cs.clone(), // 1023
      S: proof.cir_proof.S.clone(),   //31
      T: proof.cir_proof.T.clone(),   //31
      e_vec: e_vec_bytes,   //31
      tr: proof.cir_proof.mapping_proof.tr.clone(),
    };

    let short_proof_out: BoolSubCircuitProofShort = BoolSubCircuitProofShort {
      cir_proof: SubCircuitProof_short,
      bool_commit: proof.bool_commit.clone(),
      bool_proof: proof.bool_proof.clone(),
    };

    let encoded: Vec<u8> = bincode::serialize(&short_proof_out).unwrap();

    let expected_proof_size = short_proof_out.cir_proof.ys_prime.len() * 17 // ys_prime
    + short_proof_out.cir_proof.ao_prime.len() * 8  // ao_prime
    + 1 * 32                            // B,  32bytes
    + short_proof_out.cir_proof.Cs.len() * 32       // Cs, 32bytes each
    + short_proof_out.cir_proof.S.len() * 32        // S,  32bytes each
    + short_proof_out.cir_proof.T.len() * 32        // T,  32bytes each
    + short_proof_out.cir_proof.e_vec.len() * 40    // mapping_proof.e_vec, 40bytes each for n= 2^20
    + 2 * 32 // mapping_proof.tr, 64bytes (Schnorr signature)
    + short_proof_out.bool_commit.C.len() * 32
    + 1 * 32 //short_proof_out.bool_commit.U, 32bytes
    + 2 * 8 // short_proof_out.bool_proof.y1 + y2, 8bytes each
    + short_proof_out.bool_proof.b_prime.len() * 8 
    + short_proof_out.bool_proof.poly_proof.t_col.len() * 32
    + 1 * 32; //short_proof_out.bool_proof.poly_proof.tau 32bytes //


    println!("Actual proof size (bytes): {}", expected_proof_size);
    println!("Proof size after encoding (bytes): {}", encoded.len());

    let decoded: BoolSubCircuitProofShort = bincode::deserialize(&encoded[..]).unwrap();

    let y_prime_in: Vec<Scalar> = decoded.cir_proof
      .ys_prime
      .iter()
      .map(|x| Scalar::from_bytes(&SubCircuitSystem::concat_index2(x.to_vec())).unwrap())
      .collect();
    let e_vec_biguint: Vec<BigUint> = decoded.cir_proof
      .e_vec
      .iter()
      .map(|b| BigUint::from_bytes_le(b))
      .collect();

    let proofSubCircuitIn: SubCircuitProof = SubCircuitProof {
      ys_prime: y_prime_in, // 17 * 1024 =
      ao_prime: decoded.cir_proof.ao_prime.clone(), //30
      M: decoded.cir_proof.M.clone(),             //1
      R: decoded.cir_proof.R.clone(),             //1
      Cs: decoded.cir_proof.Cs.clone(),           // 1023
      S: decoded.cir_proof.S.clone(),             //31
      T: decoded.cir_proof.T.clone(),             //31
      mapping_proof: InputMappingInternalProof {
        e_vec: e_vec_biguint,
        tr: decoded.cir_proof.tr.clone(),
      },
    };

    let proofIn = BoolSubCircuitProof {
      cir_proof: proofSubCircuitIn,
      bool_commit: decoded.bool_commit,
      bool_proof: decoded.bool_proof.clone(), //lazy short-cut
      b_prime: decoded.bool_proof.b_prime.clone(), //lazy short-cut
    };

    assert_eq!(proofIn.cir_proof.ys_prime, proof.cir_proof.ys_prime);
    assert_eq!(proofIn.cir_proof.ao_prime, proof.cir_proof.ao_prime);
    assert_eq!(proofIn.cir_proof.M, proof.cir_proof.M);
    assert_eq!(proofIn.cir_proof.R, proof.cir_proof.R);
    assert_eq!(proofIn.cir_proof.Cs, proof.cir_proof.Cs);
    assert_eq!(proofIn.cir_proof.S, proof.cir_proof.S);
    assert_eq!(proofIn.cir_proof.T, proof.cir_proof.T);
    assert_eq!(proofIn.cir_proof.mapping_proof, proof.cir_proof.mapping_proof);
    proofIn
  }
}



#[cfg(test)]
mod tests {

  use super::super::montgomery::random64_mod;
  use super::*;
  use rand::rngs::OsRng;
  use super::super::subcircuit::HYPOTHESIS_SIZE;

  #[test]
  fn test_binary_circuit() {
    let n = (2_usize).pow(20 as u32); //validate values are even numbers 2, 4, 8,....20
    
    println!("Circuit size: {}", n);

    let sc = BoolSubCircuitSystem::new(b"unit test", n, HYPOTHESIS_SIZE * 32);

    let mut prover_transcript = Transcript::new(b"unit test");

    let mut csprng: OsRng = OsRng;

    let a: Vec<u64> = vec![csprng.next_u32() as u64; HYPOTHESIS_SIZE];
    let alpha: Vec<u64> = (0..HYPOTHESIS_SIZE)
      .map(|_| random64_mod(&mut csprng))
      .collect();
    let v: Vec<Scalar> = (0..HYPOTHESIS_SIZE)
      .map(|_| Scalar::random(&mut csprng))
      .collect();

    let pedersen: Vec<CompressedGroup> = sc.cir.get_mapping().compute_pederson(&a, &v);

    let phi = Scalar::random(&mut csprng);

    let proof: BoolSubCircuitProof = sc
      .prove(&a, &alpha, &v, phi, n, &mut prover_transcript)
      .unwrap();

    let proof_from_file = sc.simulate_write_read_file(proof, (n as f64).sqrt() as usize);

    
    //assert_eq!(proof_from_file, proof);
    
    let mut verifier_transcript = Transcript::new(b"unit test");

    assert_eq!(
      sc.verify(&pedersen, n, &proof_from_file, &mut verifier_transcript),
      Ok(())
    );
    
  }
}