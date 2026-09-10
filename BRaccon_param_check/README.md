# BRaccoon Parameters check

This folder contains a standalone Sage certificate for the BRaccoon concrete
parameters. It reproduces the security checks.
## Dependencies

- SageMath
- A copy of `lattice-estimator-main` is bundled in this folder and is used by
  default.

## Run

From the project root:

```bash
sage BRaccon_param_check/verify_braccoon_parameters.sage
```

If one wants to use their own lattice estimator instead of the bundled
copy, they can pass it explicitly:

```bash
sage verify_braccoon_parameters.sage --estimator-path /path/to/lattice-estimator-main
```

The Sage file is not a parameter search. It contains the three concrete rows
from the current paper and calls the lattice estimator directly on the fixed
instances used for:

- verification-key MLWE,
- `w'`-MLWE,
- signature MSIS,
- binary MLWE-BGV IND-CPA security,
- binding, hiding, and simulated-setup security of the extractable commitment.

It also recomputes the current rounded-signature reduction bound, the shifted
commitment entropy conditions, BGV and extraction no-wrap inequalities, the
mixed-R1CS field check, signature/public-key/ciphertext/commitment sizes, the
SHA3-256/SHAKE256 nonlinear constraint count, and the communication table.
The extractable commitment uses binary opening dimension `rho_ext=15`, the
smallest tested value whose hiding-MLWE estimate exceeds 128 bits for all three
rows.

The approximately 110 KB size of each ZK-LaBRADOR proof is an analytical input
from the cited implementation paper; it is not derived by this script. The BGV
correctness calculation includes the flooding bound `B_flood=2^40`. This
certifies decryption no-wrap under that bound, but not by itself the statistical
translation property of the complete flooding-coin distribution required for
chosen-key circuit privacy.

## Useful options

```bash
sage verify_braccoon_parameters.sage --full-estimator
sage verify_braccoon_parameters.sage --summary-only
sage verify_braccoon_parameters.sage --json-out results.json
```

All sizes are reported in decimal KB, i.e. bits divided by 8000.
