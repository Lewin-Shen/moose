Mission: make MOOSE's phase_field module run on GPUs through the Kokkos backend.

First checkpoint: the grain growth model — multiple non-conserved order parameters
evolving by Allen-Cahn dynamics. The approach follows decisions from a prior
capability survey: assembled CSR matrices rather than matrix-free, GAMG
preconditioning, and automatic differentiation computed inline in kernels.
Correctness gates come first: Jacobian verification, analytical benchmarks
(circular-grain and half-loop area decay), and CPU-versus-GPU comparison on
identical input decks. Design choices favor patterns that generalize across the
phase_field module over model-specific shortcuts.

This file is a branch-mechanics test artifact and will be removed before any
upstream contribution.
