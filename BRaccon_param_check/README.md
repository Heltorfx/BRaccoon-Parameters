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

The Sage file is not a parameter search.  It contains the concrete parameter
rows from the paper and calls the lattice estimator directly on the resulting
fixed MLWE/MSIS instances:

- verification-key MLWE,
- `w'`-MLWE,
- signature MSIS.

It also recomputes the BRaccoon bounds, signature/public-key sizes, and the R1CS
constraint estimate for `L_1,2`, using the current convention with one
in-circuit Poseidon hash `H(mu,w)`, where `mu=H(msg,vk)` is a `2*lambda`-bit
digest.

## Useful options

```bash
sage verify_braccoon_parameters.sage --no-r1cs
sage verify_braccoon_parameters.sage --full-estimator
sage verify_braccoon_parameters.sage --summary-only
sage verify_braccoon_parameters.sage --json-out results.json
```

All sizes are reported in decimal KB, i.e. bits divided by 8000.
