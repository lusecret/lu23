#![feature(int_log)]
use criterion::*;
use libsubcircuit_64bits::groups::group::CompressedGroup;
use libsubcircuit_64bits::montgomery::random64_mod;
use libsubcircuit_64bits::scalar::Scalar;
use libsubcircuit_64bits::subcircuit;
use libsubcircuit_64bits::bool_subcircuit;
use merlin::Transcript;
use rand::rngs::OsRng;
use rand_core::RngCore;

fn subcircuit_prove_benchmark(c: &mut Criterion) {
  for &s in [10, 20].iter() {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let sc = subcircuit::SubCircuitSystem::new(b"proof", num_degrees);

    let mut csprng: OsRng = OsRng;
    let a: Vec<u64> = vec![random64_mod(&mut csprng); subcircuit::HYPOTHESIS_SIZE];
    let alpha: Vec<u64> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| random64_mod(&mut csprng))
      .collect();
    let v: Vec<Scalar> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| Scalar::random(&mut csprng))
      .collect();
    let phi = Scalar::random(&mut csprng);

    let mut prover_transcript = Transcript::new(b"proof");

    let name = format!("subcircuit_prove_{}", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        let _ = sc.prove(
          black_box(&a),
          black_box(&alpha),
          black_box(&v),
          black_box(phi),
          black_box(num_degrees),
          black_box(&mut prover_transcript),
        );
      });
    });
    group.finish();
  }
}

fn subcircuit_verify_benchmark(c: &mut Criterion) {
  for &s in [10, 20].iter() {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let sc = subcircuit::SubCircuitSystem::new(b"proof", num_degrees);

    let mut csprng: OsRng = OsRng;
    let a: Vec<u64> = vec![random64_mod(&mut csprng); subcircuit::HYPOTHESIS_SIZE];
    let alpha: Vec<u64> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| random64_mod(&mut csprng))
      .collect();
    let v: Vec<Scalar> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| Scalar::random(&mut csprng))
      .collect();
    let pedersen: Vec<CompressedGroup> = sc.get_mapping().compute_pederson(&a, &v);
    let phi = Scalar::random(&mut csprng);

    let mut prover_transcript = Transcript::new(b"proof");
    let proof = sc
      .prove(&a, &alpha, &v, phi, num_degrees, &mut prover_transcript)
      .unwrap();

    let name = format!("subcircuit_verify_{}", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        let mut verifier_transcript = Transcript::new(b"proof");
        assert_eq!(
          sc.verify(
            black_box(&pedersen),
            black_box(num_degrees),
            black_box(&proof),
            black_box(&mut verifier_transcript)
          ),
          Ok(())
        );
      });
    });
    group.finish();
  }
}

fn binary_circuit_prove_benchmark(c: &mut Criterion) {
  for &s in [20].iter() {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let sc = bool_subcircuit::BoolSubCircuitSystem::new(b"proof", num_degrees, subcircuit::HYPOTHESIS_SIZE * 32);

    let mut csprng: OsRng = OsRng;
    let a: Vec<u64> = vec![csprng.next_u32() as u64; subcircuit::HYPOTHESIS_SIZE];
    let alpha: Vec<u64> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| random64_mod(&mut csprng))
      .collect();
    let v: Vec<Scalar> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| Scalar::random(&mut csprng))
      .collect();
    let phi = Scalar::random(&mut csprng);

    let mut prover_transcript = Transcript::new(b"proof");

    let name = format!("binary_circuit_prove for {} sequential multiplications", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        let _ = sc.prove(
          black_box(&a),
          black_box(&alpha),
          black_box(&v),
          black_box(phi),
          black_box(num_degrees),
          black_box(&mut prover_transcript),
        );
      });
    });
    group.finish();
  }
}

fn binary_circuit_verify_benchmark(c: &mut Criterion) {
  for &s in [20].iter() {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let sc = bool_subcircuit::BoolSubCircuitSystem::new(b"proof", num_degrees, subcircuit::HYPOTHESIS_SIZE * 32);

    let mut csprng: OsRng = OsRng;
    let a: Vec<u64> = vec![csprng.next_u32() as u64; subcircuit::HYPOTHESIS_SIZE];
    let alpha: Vec<u64> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| random64_mod(&mut csprng))
      .collect();
    let v: Vec<Scalar> = (0..subcircuit::HYPOTHESIS_SIZE)
      .map(|_| Scalar::random(&mut csprng))
      .collect();
    let pedersen: Vec<CompressedGroup> = sc.cir.get_mapping().compute_pederson(&a, &v);
    let phi = Scalar::random(&mut csprng);

    let mut prover_transcript = Transcript::new(b"proof");
    let proof = sc
      .prove(&a, &alpha, &v, phi, num_degrees, &mut prover_transcript)
      .unwrap();

    let name = format!("binary_circuit_verify for {} sequential multiplications", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        let mut verifier_transcript = Transcript::new(b"proof");
        assert_eq!(
          sc.verify(
            black_box(&pedersen),
            black_box(num_degrees),
            black_box(&proof),
            black_box(&mut verifier_transcript)
          ),
          Ok(())
        );
      });
    });
    group.finish();
  }
}

fn set_duration() -> Criterion {
  Criterion::default().sample_size(10)
}

criterion_group! {
name = benches_snark;
config = set_duration();
//targets = subcircuit_prove_benchmark, subcircuit_verify_benchmark, boolsubcircuit_prove_benchmark, boolsubcircuit_verify_benchmark
targets = binary_circuit_prove_benchmark, binary_circuit_verify_benchmark
}

criterion_main!(benches_snark);
