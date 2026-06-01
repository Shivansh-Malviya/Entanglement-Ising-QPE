# Entanglement in Ising Spin Chains

This repository presents numerical C++ studies of entanglement and ground-state behaviour in one-dimensional transverse-field Ising spin chains. The work focuses on concurrence, derivative-of-concurrence signatures near the quantum phase transition, dense and sparse Hamiltonian construction, and Lanczos-based ground-state estimation.

## Repository Contents

| Path | Contents |
|---|---|
| `src/concurrence_derivative_dense.cpp` | Dense-matrix concurrence sweep and numerical derivative with respect to the field ratio. |
| `src/lanczos_dense.cpp` | Dense Lanczos eigensolver study for transverse-field Ising Hamiltonians. |
| `src/lanczos_sparse.cpp` | Sparse-matrix Lanczos eigensolver study for larger Hilbert spaces. |
| `src/concurrence_sparse_trace.cpp` | Sparse concurrence and derivative sweep using an explicit reduced-density trace routine. |
| `src/partial_trace_demo.cpp` | Small Armadillo partial-trace prototype used while developing the concurrence calculation. |
| `results/concurrence_derivative/` | Dense concurrence and derivative output tables. |
| `results/lanczos_dense/` | Dense Lanczos ground-state energy sweep outputs. |
| `results/lanczos_sparse/` | Sparse Lanczos ground-state energy sweep outputs. |
| `results/sparse_trace/` | Sparse trace-based concurrence and derivative output tables. |
| `plotting/sparse_trace_concurrence.gnuplot` | Gnuplot script for the sparse trace concurrence/derivative data family. |

## Technical Scope

The code builds transverse-field Ising Hamiltonians from Pauli operator tensor products, applies periodic boundary terms, extracts ground-state vectors through direct diagonalisation or Lanczos iterations, and evaluates two-spin concurrence from reduced density matrices. The retained result tables record representative parameter sweeps across chain sizes and solver configurations.

The C++ programs use Armadillo for linear algebra and QIClib for quantum-information routines such as subsystem traces and spin matrices. QIClib is referenced as a normal system include in the source files. Vendored dependency copies, compiled binaries, and transient run outputs are not part of this repository.

## Result Families

- `con*.dat` and `dcon*.dat`: concurrence and derivative curves from dense diagonalisation.
- `Lan_gev_*.dat`: dense Lanczos ground-state energy estimates for selected chain sizes and iteration counts.
- `LG_sp*.dat`: sparse Lanczos ground-state energy estimates.
- `sdt_con*.dat` and `sdt_dcon*.dat`: sparse trace-based concurrence and derivative curves.

## License

This project is licensed under the MIT License.
