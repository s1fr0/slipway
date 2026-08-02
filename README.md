# Slipway

Companion code for **Slipway: Accessing Finite Subspace Trails in Poseidon**.
This repository contains two distinct tools:

| Command | Purpose |
|---|---|
| `search.py` | Find a fresh matrix, control, and CICO-2 solution on a configurable reduced-round instance. |
| `verify_paper_rf8_rp20.py` | Verify the supplied paper witness at width 16, RF=8, RP=20. |

The verifier receives the matrix, constants, control, root, input, and output
from `paper_rf8_rp20.json`. **It checks an existing answer; it does not search
for or independently rediscover the paper's solution.** Its runtime is a
verification time, not the cost of finding the full-round solution.

CICO-2 means that the first two input coordinates and the first two output
coordinates are zero. Search experiments choose a matrix after the round
constants have been fixed. They demonstrate the paper's construction under
that matrix-selection model, not an attack on a fixed standardized matrix.

## Setup

Requirements: **SageMath** (tested with 10.7), Python 3, and an internet connection
for the initial dependency download. No GPU or extra Python packages are needed.

Run all commands below from this directory, which is intended to be the root
of the standalone repository:

```sh
python3 reference.py
```

This installs the unmodified
[poseidon-tools](https://github.com/khovratovich/poseidon-tools/tree/60075da7c0521d9493749a035b1f30d4eda37138)
dependency at commit `60075da7c0521d9493749a035b1f30d4eda37138`.
The installer verifies the archive's SHA-256 and keeps the source and its
license in `.deps/`. Subsequent runs work offline. Upstream has no Python
packaging metadata, so the code is imported from this pinned source directory.

All required Slipway source and fixture data are included here. The surrounding
paper repository is not needed.

## Find a fresh solution

Start with the small KoalaBear example:

```sh
sage -python search.py --rf 4 --rp 5 --seed 0 --output results/kb_rf4_rp5.json
```

This uses the KoalaBear prime `2130706433`, width 6, a cubic S-box, four
full rounds and five partial rounds. It generates constants using upstream
`Poseidon(p, 3, width, rf, rp).round_constants` **before** sampling any matrix.
The seed controls the matrix search, not the constants.

The search constructs a matrix and constrained-input family, checks every
nonempty square minor for MDS, and proves the prefix and finite-trail
identities by exact algebra. It then solves the two output equations jointly
using a resultant, recovers the control and input, and verifies the complete
permutation with upstream `Poseidon.permutation()`. No successful matrix,
control, root, or solution is supplied to this search.

An observed seed-0 run produced:

```text
Matrix attempts: 9; certified matrices: 2; resultant degree: 243
Control: 89474717; inactive partial rounds: 4; degrees: [27, 27]
CICO-2 solution X=2068155477
Input:  [0, 0, 1987899505, 1438241491, 9062786, 842432415]
Output: [0, 0, 1172977448, 1974871881, 1645683826, 1405796349]
```

The JSON report records the parameters, upstream revision, matrix, constants,
family, discovered control, complete input/output witness, search counts,
and timings. The output directory is created automatically; `results/`
is ignored by Git.

More examples:

```sh
# One additional partial round over KoalaBear.
sage -python search.py --rf 4 --rp 6 --output results/kb_rf4_rp6.json

# Six full rounds; the default width increases to nine.
sage -python search.py --rf 6 --rp 7

# A smaller field for experiments.
sage -python search.py --prime 257 --rf 4 --rp 5

sage -python search.py --help
```

### Parameters and computational limits

| Option | Default | Meaning |
|---|---|---|
| `--prime` | 2130706433 | Prime field modulus; must exceed twice the width and be 2 modulo 3. |
| `--rf` | 4 | Total full rounds; even and at least 4. |
| `--rp` | 5 | Partial rounds; nonnegative. |
| `--width` | `3*rf/2` | State width; must be at least this value for the construction. |
| `--seed` | 0 | Reproducible matrix-search seed. |
| `--max-attempts` | 10000 | Maximum matrix candidates to try. |
| `--max-degree` | 81 | Maximum predicted output degree accepted by the joint solver. |
| `--output` | none | Save the report to a JSON file. |

Changing `rf` also changes the default width. Pass `--width` explicitly
when comparing round counts at a fixed width, subject to the construction's
minimum-width requirement.

The construction absorbs the first `rf/2` full rounds and makes at least
`min(rp, width-2)` partial S-boxes inactive. The output-degree bound is

```text
3 ** (rf/2 + max(0, rp - (width-2)))
```

Exceeding `--max-degree` fails before searching. Raising the limit can make
elimination much more expensive. It is a degree limit, not a time or memory
limit. Likewise, `--max-attempts` limits candidates, not elapsed time.
The search retries when a matrix yields no solution in the parameterized
family or elimination is degenerate. Exhausting the attempt limit exits with
status 1; success is not guaranteed for every seed.

MDS verification checks `binomial(2*width, width)-1` minors: 923 at width 6,
48,619 at width 9, and 2,704,155 at width 12. This can dominate runtime at
larger widths.

### Interpreting timings

Illustrative local seed-0 runs with SageMath 10.7, KoalaBear, and width 6:

| RF | RP | Absorbed full / inactive partial | Output degree | Complete search time |
|---|---|---|---|---|
| 4 | 5 | 2 / 4 | 27 | about 0.3 s |
| 4 | 6 | 2 / 4 | 81 | about 34 s |

These include constant generation, all matrix attempts and checks, solving,
and complete verification. They exclude Sage startup and writing the report.
They are individual observations, not general runtime guarantees or a
controlled benchmark; machine load and seed affect the result.

For a reproducible timing report, retain the JSON and record the hardware,
Sage version, all parameters, and seed. Report the matrix attempts as well as
elapsed time. Verification of a supplied paper witness is a different task
and should not be presented as the time needed to find that solution.

## Verify the supplied RF=8, RP=20 paper witness

```sh
sage -python verify_paper_rf8_rp20.py
```

This checks the fixed KoalaBear instance with width 16, eight full rounds, and
twenty partial rounds. The included fixture explicitly supplies:

- the matrix, round constants, and input-family data;
- control `1203195526` and root `1404166674`;
- the complete 16-coordinate input and output reported in the paper.

The verifier checks input recovery against the supplied input, evaluates every
permutation round with the reference dependency, and compares the entire output
with the supplied output. It also checks the symbolic four-round identity,
`dim K14 = 2`, characteristic-polynomial irreducibility through `M^64`,
the fourteen inactive partial S-boxes, both degree-59,049 output polynomials,
and their evaluation at the **supplied** root.

It performs **no matrix search, control search, or root extraction**.
It does not recheck the 601,080,389 MDS minors of this 16-wide matrix.
Success prints `VERIFIED` and exits with status 0; a failed check exits with
status 1. Verification time is printed separately.

The fixture corresponds to the paper's “Concrete CICO-2 Solution” section and
“CICO-2 Matrix, Parameters, and Solution” appendix. Its reported constants
come from the original Poseidon Grain procedure and differ from the pinned
dependency's default generator. Those fixed values are passed explicitly
to the reference permutation. **New searches always generate their constants
through the dependency.**

Fresh searches at the paper's full parameters are outside the intended desktop
examples. The repository does not include a CICO-3 solver, a zero-test challenge
solver, or a GPU/NTT implementation.

## Tests and repository layout

```sh
sage -python -m unittest discover -v
```

Tests cover fresh searches with different fields, seeds and round counts,
upstream constants and permutation verification, polynomial identities,
control recovery, invalid parameters, and rejection of modified paper witnesses.
They also check that the paper verifier does not invoke either search solver.
The routine RF=8 construction test skips exhaustive MDS enumeration; smaller
search tests exercise the complete MDS check.

| File | Role |
|---|---|
| `search.py` | Fresh matrix and solution search; command-line parameters. |
| `verify_paper_rf8_rp20.py` | Verification-only command for the supplied paper witness. |
| `paper_rf8_rp20.json` | Self-contained paper instance and explicit witness. |
| `core.py` | Shared algebra, certificate checks, and search solvers. |
| `reference.py` | Pinned dependency installer and loader. |
| `test_slipway.py` | Regression and integration tests. |

Downloaded dependencies, Python caches, and generated results are excluded
by `.gitignore`.

## Licensing

Licensing for this repository's own code is **undecided**.
The third-party poseidon-tools dependency retains its Apache-2.0 license.
