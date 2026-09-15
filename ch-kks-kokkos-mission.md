Mission: make MOOSE's phase_field module run on GPUs through the Kokkos backend.

First checkpoint, done on `gg-kokkos` (2026-09-09): the grain growth model, multiple
non-conserved order parameters evolving by Allen-Cahn dynamics. Three Kokkos objects
(KokkosACGrGrPoly, KokkosACInterface, KokkosGBEvolution) with hand-coded Jacobians,
gated by Jacobian verification, analytical benchmarks and CPU-versus-GPU comparison
on identical decks, then measured on one B200 against one CPU node (Campaign GPU-01).

Second checkpoint, this branch (`ch-kks-kokkos`, cut from `gg-kokkos` on 2026-09-15):
Cahn-Hilliard in its split forms and the KKS model. The vehicle is the spinodal
decomposition deck of the phase-field tutorial. What is different from the first
checkpoint: the free energy and its derivatives must be written in closed form for the
device (Kokkos-MOOSE has no derivative materials), the linear solve is the bottleneck
rather than assembly, and the solver has to be block-structured (FieldSplit or a Schur
form) because Jacobi-class preconditioners and whole-system AMG fail on the implicit
step. Correctness gates come first, as before: Jacobian verification, CPU-versus-Kokkos
comparison on identical no-flux decks, free-energy decay and mass conservation, then the
KKS examples. Design choices favor patterns that generalize across the phase_field
module over model-specific shortcuts. The design document lives in the research
repository (`port/DESIGN-kokkos-ch-kks.md`).

This file is a branch note and will be removed before any upstream contribution.
