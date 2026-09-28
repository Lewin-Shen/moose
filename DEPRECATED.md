# This branch is deprecated (2026-09-27)

`ch-kks-kokkos` was cut from `gg-kokkos` on 2026-09-15, so it carries the preliminary Kokkos grain growth port
(`b660ac01b9`) beneath the Cahn-Hilliard port (`13d336e580`). The two ports are unrelated code, and stacking them
made the Cahn-Hilliard work impossible to review or submit on its own.

The Cahn-Hilliard work continues on `ch-kokkos`: the same port commit cherry-picked onto the pinned upstream commit
`6956830688`, with no grain growth commit and no branch note. The grain growth port stays on `gg-kokkos`.

This branch is kept as a record and receives no new work. Its commit `13d336e580` is the build every row of the
second GPU campaign ran on, and the parent of the campaign branch `gpu03-gg-predeploy`. Do not delete it while those
records cite it.
