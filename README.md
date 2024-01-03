# Submission 17

## Get code compiled

1. Open your terminal and CD to poly-mul-over-64bits-prime folder and type `make lib` command to create the C library libpolyntt64bits.dylib (Mac) or libpolyntt64bits.so (Linux). While our protocol is coded in Rust, the comput_(sub)circuit_keys operations are mostly coded in C to leverage open-source code for fast NTT implementation ( https://github.com/IBM/optimized-number-theoretic- transform-implementations).

3. Make sure the  .so or .dylib file you created from the earlier step is in the path of rustc.
   
4. Goto subcircuit-64bits directory and type `cargo bench` to benchmark the prover/verifier runtime and type `cargo test test_binary_circuit -- --show-output ` command to get the communication cost printed on the screen

 
## Benchmark - Prover and Verifier Runtime

As mentioned before, the benchmark is significantly better than that reported in the submission mainly because `"-C", "target-cpu=native",` flag wasn't turned on in the earlier benchmark testing (and some minor optimizations).

![16871704294520_ pic](https://github.com/lusecret/lu23/assets/8139291/166a618c-9707-4a7a-9535-597225b5753a)

commits_verify_benchmark/binary_circuit_prove for 1048576 sequential multiplications 

time:   [2.5878 s 2.6443 s 2.7350 s]

Found 2 outliers among 10 measurements (20.00%)
  1 (10.00%) high mild
  1 (10.00%) high severe

Benchmarking commits_verify_benchmark/binary_circuit_verify for 1048576 sequential multiplications: Collecting 10 samples in                                                         

commits_verify_benchmark/binary_circuit_verify for 1048576 sequential multiplications     

time:   [16.967 ms 17.030 ms 17.096 ms]

## Benchmark - Communication Cost

The "actual" benchmark is computed by adding up all elements (group and field) used to create the proof transcript, this is also the theoretical (achievable) communication cost. The "encoded" proof size includes additional overheads to allow encoded streams to quickly convert back to a rust struct instance (We used bincode crate to serialize/deserialize the proof transcript). 

![Screenshot 2024-01-03 at 11 45 35 PM](https://github.com/lusecret/lu23/assets/8139291/2f0a2925-52b0-4743-a20c-501239072f43)

---- bool_subcircuit::tests::test_binary_circuit stdout ----
Circuit size: 1048576

Actual proof size (bytes): 64408

Proof size after encoding (bytes): 73957

## Random circuit

The random circuit we use in benchmarking adds-up all input bits (b') and then performs `n` sequential multiplications on itself.
