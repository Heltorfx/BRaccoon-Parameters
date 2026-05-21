#!/usr/bin/env sage
#
# Standalone verification script for the BRaccoon parameter tables.
#
# Dependencies:
#   - SageMath
#   - the lattice-estimator repository, passed with --estimator-path
#
# Example from the project root:
#   sage reviewer_sage/verify_braccoon_parameters.sage --estimator-path lattice-estimator-main
#
# Example if this file and lattice-estimator-main are shipped side-by-side:
#   sage verify_braccoon_parameters.sage --estimator-path lattice-estimator-main

import argparse
import json
import math
import pprint
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from sage.all import next_prime, oo


SCRIPT_DIR = Path(__file__).resolve().parent
TABLE_Q = 302231454903657293688833


@dataclass(frozen=True)
class ParameterSet:
    name: str
    max_signatures: int
    ring_degree: int = 256
    q: int = TABLE_Q
    q_bits: int = 78
    module_dim: int = 12
    secret_dim: int = 12
    challenge_weight: int = 23
    nu_t: int = 52
    nu_w: int = 55
    sigma_t: int = 2**10
    sigma_wprime: int = 2**17
    sigma_w: int = 2**57


PARAMETER_SETS = (
    ParameterSet("Q_s=2^20", max_signatures=2**20, secret_dim=14),
    ParameterSet("Q_s=2^32", max_signatures=2**32, secret_dim=16),
    ParameterSet("Q_s=2^64", max_signatures=2**64, secret_dim=22),
)


def resolve_path(path_text):
    path = Path(path_text)
    candidates = []
    if path.is_absolute():
        candidates.append(path)
    else:
        # Prefer the estimator bundled next to this certificate.  Reviewers can
        # still override it with --estimator-path /path/to/their/estimator.
        candidates.append(SCRIPT_DIR / path)
        candidates.append(Path.cwd() / path)
        candidates.append(SCRIPT_DIR.parent / path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def load_estimator(estimator_path):
    root = resolve_path(estimator_path)
    if not root.exists():
        raise FileNotFoundError(
            "Estimator path not found: {}. Pass --estimator-path /path/to/lattice-estimator-main".format(
                root
            )
        )
    sys.path.insert(0, str(root))
    from estimator import LWE, SIS, ND

    return LWE, SIS, ND


def ntt_prime_at_least(value, ring_degree):
    modulus = 2 * ring_degree
    candidate = max(2, int(value))
    while True:
        candidate = int(next_prime(candidate - 1))
        if candidate % modulus == 1 and math.gcd(candidate - 1, modulus) != 1:
            return candidate
        candidate += 1


def tail_bound_inf(rejection_rate):
    tail = 0.01
    log_rej = math.log(rejection_rate)
    while math.log(2.0) - (tail**2) / 2.0 > log_rej:
        tail += 0.01
    return tail


def tail_bound_l2(ring_degree, rejection_rate):
    tail = 1.01
    log_rej = math.log(rejection_rate)
    while ring_degree * math.log(tail) + (ring_degree / 2.0) * (1.0 - tail**2) > log_rej:
        tail += 0.01
    return tail


def compute_profile(params, q, target_bits):
    n = params.ring_degree
    k = params.module_dim
    ell = params.secret_dim
    omega = params.challenge_weight
    sessions = params.max_signatures
    nu_t = params.nu_t
    nu_w = params.nu_w
    sigma_t = params.sigma_t
    sigma_wprime = params.sigma_wprime
    sigma_w = params.sigma_w

    # Same notation as the paper/scripts: n=t=1 for the Raccoon threshold sizes.
    threshold_n = 1.0
    threshold_t = 1.0
    rounding_lhs = omega * (2**nu_t) + 2 ** (nu_w + 1)

    b_hmlwe = sessions * omega * (
        1.0 + n * (target_bits + 1.0 + 2.0 * math.log2(n)) / math.sqrt(sessions)
    )
    sigma_mlwe = math.sqrt(
        1.0
        / (
            2.0 / (threshold_n * sigma_t**2)
            + 2.0 * b_hmlwe / (threshold_t * sigma_wprime**2)
        )
    )

    module_sqrt = math.sqrt(n * (k + ell))
    moddim_sqrt = math.sqrt(n * k)
    b2_star = (
        math.exp(0.25)
        * (threshold_n * omega * sigma_t + threshold_t * sigma_wprime + sigma_w)
        * module_sqrt
        + rounding_lhs * moddim_sqrt
    )
    b_stmsis = b2_star + math.sqrt(omega) + rounding_lhs * moddim_sqrt
    b_msis = b_stmsis - omega
    q_required = math.ceil(2.0 * b_msis + 1.0)
    q_ntt_required = ntt_prime_at_least(q_required, n)

    rejection_rate = 0.01
    tail_inf = tail_bound_inf(rejection_rate)
    b_inf = tail_inf * sigma_t * math.sqrt(n) + tail_inf * sigma_wprime + tail_inf * sigma_w
    challenge_bits = n
    response_bits = n * ell * math.ceil(math.log2(b_inf))
    hint_bits = n * k * math.ceil(math.log2(b2_star / ((2**nu_w) * moddim_sqrt)))
    signature_bits = challenge_bits + response_bits + hint_bits
    public_key_bits = n * k * (math.log2(q) - nu_t) + target_bits

    return {
        "B_HMLWE": b_hmlwe,
        "sigma_MLWE": sigma_mlwe,
        "B_2_star": b2_star,
        "B_STMSIS": b_stmsis,
        "B_MSIS": b_msis,
        "q_required_bits": math.log2(q_required),
        "q_ntt_required": q_ntt_required,
        "q_ntt_required_bits": math.log2(q_ntt_required),
        "q_satisfies_msis": b_msis < (q - 1) / 2.0,
        "signature_kb": signature_bits / 8000.0,
        "public_key_kb": public_key_bits / 8000.0,
    }


def min_log2_rop(result):
    if not result:
        return None
    values = result.values() if hasattr(result, "values") else result
    bits = []
    for item in values:
        try:
            rop = item.get("rop", None)
        except AttributeError:
            rop = None
        if rop is None:
            continue
        try:
            rop_float = float(rop)
        except Exception:
            continue
        if rop_float > 0:
            bits.append(math.log2(rop_float))
    return min(bits) if bits else None


def estimate_lwe(LWE, ND, n, q, xs_sigma, xe_sigma, m, tag, rough):
    lwe_params = LWE.Parameters(
        n=n,
        q=q,
        Xs=ND.DiscreteGaussian(xs_sigma),
        Xe=ND.DiscreteGaussian(xe_sigma),
        m=m,
        tag=tag,
    )
    fn = LWE.estimate.rough if rough else LWE.estimate
    result = fn(lwe_params, quiet=True)
    return min_log2_rop(result), result


def estimate_sis(SIS, n, q, m, length_bound, tag, rough):
    sis_params = SIS.Parameters(
        n=n,
        q=q,
        m=m,
        length_bound=length_bound,
        norm=2,
        tag=tag,
    )
    fn = SIS.estimate.rough if rough else SIS.estimate
    result = fn(sis_params, quiet=True)
    return min_log2_rop(result), result


def ceil_log2_int(value):
    if value <= 1:
        return 0
    return int(math.ceil(math.log2(float(value))))


def log2_binomial(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2)


def poseidon_rounds(field_bits, state_size, alpha, security_bits):
    def satisfied(round_f, round_p):
        sec = security_bits
        t = state_size
        log_alpha = math.log(alpha, 2)
        rf1 = 6 if sec <= (math.floor(field_bits - ((alpha - 1) / 2.0)) * (t + 1)) else 10
        rf2 = (
            1
            + math.ceil(math.log(2, alpha) * min(sec, field_bits))
            + math.ceil(math.log(t, alpha))
            - round_p
        )
        rf3 = math.log(2, alpha) * min(sec, field_bits) - round_p
        rf4 = t - 1 + math.log(2, alpha) * min(sec / float(t + 1), field_bits / 2.0) - round_p
        rf5 = (t - 2 + (sec / float(2 * log_alpha)) - round_p) / float(t - 1)
        rf_max = max(math.ceil(rf1), math.ceil(rf2), math.ceil(rf3), math.ceil(rf4), math.ceil(rf5))

        r_temp = math.floor(t / 3.0)
        over = int(round((round_f - 1) * t + round_p + r_temp + r_temp * (round_f / 2.0) + round_p + alpha))
        under = int(round(r_temp * (round_f / 2.0) + round_p + alpha))
        binom_log = log2_binomial(over, under)
        if not math.isfinite(binom_log):
            binom_log = sec + 1
        cost_gb4 = math.ceil(2 * binom_log)
        return round_f >= rf_max and cost_gb4 >= sec

    best = None
    min_cost = float("inf")
    max_cost_rf = 0
    for round_p in range(1, 500):
        for round_f in range(4, 100):
            if round_f % 2 != 0 or not satisfied(round_f, round_p):
                continue
            rf = round_f + 2
            rp = int(math.ceil(float(round_p) * 1.075))
            cost = state_size * rf + rp
            if cost < min_cost or (cost == min_cost and rf < max_cost_rf):
                best = (rf, rp)
                min_cost = cost
                max_cost_rf = rf
    if best is None:
        raise ValueError("Could not find Poseidon rounds.")
    return best


def poseidon_constraints(field_bits, security_bits, state_size, alpha, capacity):
    cap = capacity if capacity is not None else int(math.ceil((2.0 * security_bits) / field_bits))
    if cap <= 0 or cap >= state_size:
        raise ValueError("Invalid Poseidon capacity for state size.")
    full_rounds, partial_rounds = poseidon_rounds(field_bits, state_size, alpha, security_bits)
    sbox_cost = ceil_log2_int(alpha)
    per_perm = (
        2 * state_size * (full_rounds + partial_rounds)
        + state_size * full_rounds * sbox_cost
        + partial_rounds * sbox_cost
    )
    return {
        "field_bits": field_bits,
        "state_size": state_size,
        "capacity": cap,
        "rate": state_size - cap,
        "alpha": alpha,
        "full_rounds": full_rounds,
        "partial_rounds": partial_rounds,
        "permutation_constraints": int(per_perm),
    }


def digest_element_count(security_bits, field_bits):
    # mu=H(msg,vk) is a 2*lambda-bit digest, i.e. 256 bits for lambda=128.
    return max(1, int(math.ceil((2.0 * security_bits) / field_bits)))


def rounded_value_inf_bound(q, nu):
    return math.ceil(((q - 1) / 2.0) / (2**nu)) + 1


def packed_element_count(num_values, value_bits, field_modulus):
    bits = max(1, int(value_bits))
    pack = 1
    while (1 << ((pack + 1) * bits)) < field_modulus:
        pack += 1
    return int(math.ceil(num_values / pack)), pack


def auto_hash_state_candidates(field_bits, security_bits, input_len, capacity):
    cap = capacity if capacity is not None else int(math.ceil((2.0 * security_bits) / field_bits))
    candidates = {max(cap + 2, 8), 16, 32, 64, 128, 256, 512}
    if input_len > 0:
        for target_perms in range(1, 9):
            candidates.add(cap + int(math.ceil(input_len / target_perms)))
    max_auto_state = max(256, min(512, input_len + cap))
    return tuple(sorted(state for state in candidates if cap < state <= max_auto_state))


def select_poseidon_hash(field_bits, target_bits, input_len, args):
    hash_states = args.r1cs_hash_states
    if not hash_states:
        hash_states = auto_hash_state_candidates(
            field_bits,
            target_bits,
            input_len,
            args.r1cs_hash_capacity,
        )

    best_hash = None
    for state in hash_states:
        hp = poseidon_constraints(
            field_bits,
            target_bits,
            state,
            args.r1cs_hash_alpha,
            args.r1cs_hash_capacity,
        )
        perms = int(math.ceil(input_len / hp["rate"]))
        constraints = perms * hp["permutation_constraints"]
        candidate = (constraints, perms, hp)
        if best_hash is None or candidate[0] < best_hash[0]:
            best_hash = candidate
    return best_hash


def estimate_r1cs(params, q, target_bits, args):
    n = params.ring_degree
    k = params.module_dim
    ell = params.secret_dim
    field_bits = max(1, int(math.ceil(math.log2(q))))

    bdlop_rows_a = args.r1cs_bdlop_rows_a
    bdlop_message_len = args.r1cs_bdlop_message_len
    if bdlop_message_len is None:
        bdlop_message_len = ell + 2
    bdlop_randomness_len = args.r1cs_bdlop_randomness
    gamma_bit_bound = 0  # binary opening, so coefficients are already boolean.
    bdlop_constraints = n * (
        bdlop_rows_a
        + 2 * bdlop_message_len
        + 2 * bdlop_randomness_len
        + gamma_bit_bound * bdlop_randomness_len
    )

    mu_elements = args.r1cs_mu_elements
    if mu_elements is None:
        mu_elements = digest_element_count(target_bits, field_bits)

    msg_elements = args.r1cs_msg_elements
    vk_elements = args.r1cs_vk_elements
    if vk_elements is None:
        vk_elements = digest_element_count(target_bits, field_bits)

    digest_hash_input_len = msg_elements + vk_elements
    digest_hash_constraints, digest_hash_perms, digest_hash_poseidon = select_poseidon_hash(
        field_bits,
        target_bits,
        digest_hash_input_len,
        args,
    )

    w_bound = rounded_value_inf_bound(q, params.nu_w)
    w_coeff_bits = ceil_log2_int(2 * w_bound + 1)
    auto_w_elements, w_pack_coeffs = packed_element_count(k * n, w_coeff_bits, q)
    w_elements = args.r1cs_w_elements if args.r1cs_w_elements is not None else auto_w_elements

    challenge_hash_input_len = mu_elements + w_elements
    challenge_hash_constraints, challenge_hash_perms, challenge_hash_poseidon = select_poseidon_hash(
        field_bits,
        target_bits,
        challenge_hash_input_len,
        args,
    )

    xof_poseidon = poseidon_constraints(
        field_bits,
        target_bits,
        args.r1cs_xof_state,
        args.r1cs_xof_alpha,
        args.r1cs_xof_capacity,
    )
    xof_output_elements = params.challenge_weight + int(math.ceil(params.challenge_weight / 8.0))
    xof_perms = int(math.ceil(xof_output_elements / xof_poseidon["rate"]))
    xof_constraints = xof_perms * xof_poseidon["permutation_constraints"]
    sample_constraints = (
        xof_output_elements * (field_bits + 1)
        + params.challenge_weight * (3 * n + 3)
        + params.challenge_weight
    )
    total = (
        bdlop_constraints
        + digest_hash_constraints
        + challenge_hash_constraints
        + xof_constraints
        + sample_constraints
    )

    return {
        "total_constraints": int(total),
        "total_log2": math.log2(total),
        "bdlop_constraints": int(bdlop_constraints),
        "hash_vk_msg_constraints": int(digest_hash_constraints),
        "hash_mu_w_constraints": int(challenge_hash_constraints),
        "sample_in_ball_constraints": int(sample_constraints),
        "poseidon_xof_constraints": int(xof_constraints),
        "mu_elements": int(mu_elements),
        "msg_elements": int(msg_elements),
        "vk_elements": int(vk_elements),
        "vk_seed_bits": int(2 * target_bits),
        "w_elements": int(w_elements),
        "w_coeff_bits": int(w_coeff_bits),
        "w_pack_coeffs": int(w_pack_coeffs),
        "hash_vk_msg_input_len": int(digest_hash_input_len),
        "hash_mu_w_input_len": int(challenge_hash_input_len),
        "hash_vk_msg_permutations": int(digest_hash_perms),
        "hash_mu_w_permutations": int(challenge_hash_perms),
        "xof_output_elements": int(xof_output_elements),
        "xof_permutations": int(xof_perms),
        "bdlop_rows_a": int(bdlop_rows_a),
        "bdlop_message_len": int(bdlop_message_len),
        "bdlop_randomness_len": int(bdlop_randomness_len),
        "hash_vk_msg_poseidon": digest_hash_poseidon,
        "hash_mu_w_poseidon": challenge_hash_poseidon,
        "xof_poseidon": xof_poseidon,
    }


def fmt_bits(value):
    if value is None:
        return "n/a"
    if math.isinf(value):
        return "inf bits"
    return "{:.2f} bits".format(value)


def fmt_power(value):
    if value <= 0:
        return "0"
    logv = math.log2(float(value))
    if abs(logv - round(logv)) < 1e-12 and logv >= 10:
        return "2^{}".format(int(round(logv)))
    return "2^{:.2f} (~{:,.4g})".format(logv, float(value))


def ascii_estimator_output(value):
    text = pprint.pformat(value, width=100)
    replacements = {
        "≈": "~",
        "β": "beta",
        "β'": "beta'",
        "δ": "delta",
        "ζ": "zeta",
        "η": "eta",
        "σ": "sigma",
        "α": "alpha",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def parse_int_or_auto(text):
    if text is None:
        return None
    value = str(text).strip().lower()
    if value in ("", "auto"):
        return None
    return int(eval(value, {"__builtins__": {}}, {}))


def parse_range_or_auto(text):
    value = str(text).strip().lower()
    if value in ("", "auto"):
        return ()
    if "," in value:
        return tuple(int(part.strip()) for part in value.split(",") if part.strip())
    if ":" in value:
        parts = [int(part.strip()) for part in value.split(":")]
        if len(parts) == 2:
            return tuple(range(parts[0], parts[1] + 1))
        if len(parts) == 3:
            return tuple(range(parts[0], parts[1] + 1, parts[2]))
    return (int(value),)


def verify_parameter_set(params, LWE, SIS, ND, args):
    q = params.q
    profile = compute_profile(params, q, args.target)
    rough = not args.full_estimator

    vk_mlwe, vk_mlwe_raw = estimate_lwe(
        LWE,
        ND,
        n=params.ring_degree * params.secret_dim,
        q=q,
        xs_sigma=profile["sigma_MLWE"],
        xe_sigma=profile["sigma_MLWE"],
        m=params.ring_degree * params.module_dim,
        tag=params.name + " verification-key MLWE",
        rough=rough,
    )
    wprime_mlwe, wprime_mlwe_raw = estimate_lwe(
        LWE,
        ND,
        n=params.ring_degree * params.secret_dim,
        q=q,
        xs_sigma=params.sigma_wprime,
        xe_sigma=params.sigma_w,
        m=params.ring_degree * params.module_dim,
        tag=params.name + " wprime MLWE",
        rough=rough,
    )
    signature_msis, signature_msis_raw = estimate_sis(
        SIS,
        n=params.ring_degree * params.module_dim,
        q=q,
        m=params.ring_degree * (params.module_dim + params.secret_dim + 1),
        length_bound=profile["B_MSIS"],
        tag=params.name + " signature MSIS",
        rough=rough,
    )

    security_values = [vk_mlwe, wprime_mlwe, signature_msis]
    minimum = min(value for value in security_values if value is not None)
    r1cs = None if args.no_r1cs else estimate_r1cs(params, q, args.target, args)
    passes = minimum >= args.target and profile["q_satisfies_msis"]

    return {
        "parameters": asdict(params),
        "q": q,
        "q_log2": math.log2(q),
        "profile": profile,
        "security": {
            "verification_key_mlwe": vk_mlwe,
            "wprime_mlwe": wprime_mlwe,
            "signature_msis": signature_msis,
            "minimum": minimum,
            "passes": passes,
        },
        "estimator_outputs": {
            "verification_key_mlwe": ascii_estimator_output(vk_mlwe_raw),
            "wprime_mlwe": ascii_estimator_output(wprime_mlwe_raw),
            "signature_msis": ascii_estimator_output(signature_msis_raw),
        },
        "r1cs_l1_2": r1cs,
    }


def print_result(result, target):
    p = result["parameters"]
    profile = result["profile"]
    security = result["security"]
    status = "PASS" if security["passes"] else "FAIL"
    print("{} [{}]".format(p["name"], status))
    print(
        "  algebraic: n={}, log2(q)={:.2f}, k={}, ell={}, omega={}".format(
            p["ring_degree"], result["q_log2"], p["module_dim"], p["secret_dim"], p["challenge_weight"]
        )
    )
    print(
        "  gaussian/rounding: sigma_t={}, sigma_wprime={}, sigma_w={}, nu_t={}, nu_w={}".format(
            fmt_power(p["sigma_t"]),
            fmt_power(p["sigma_wprime"]),
            fmt_power(p["sigma_w"]),
            p["nu_t"],
            p["nu_w"],
        )
    )
    print(
        "  derived: B_HMLWE={}, sigma_MLWE={}, B_MSIS={}, q valid={}".format(
            fmt_power(profile["B_HMLWE"]),
            fmt_power(profile["sigma_MLWE"]),
            fmt_power(profile["B_MSIS"]),
            "yes" if profile["q_satisfies_msis"] else "no",
        )
    )
    print(
        "  sizes: vk={:.2f} KB, sig={:.2f} KB".format(
            profile["public_key_kb"], profile["signature_kb"]
        )
    )
    print(
        "  security: VK-MLWE={}, wprime-MLWE={}, SIG-MSIS={}, min={}, target={:.0f} bits".format(
            fmt_bits(security["verification_key_mlwe"]),
            fmt_bits(security["wprime_mlwe"]),
            fmt_bits(security["signature_msis"]),
            fmt_bits(security["minimum"]),
            target,
        )
    )
    r1cs = result["r1cs_l1_2"]
    if r1cs is not None:
        print(
            "  R1CS L1,2: total={} (2^{:.2f}), BDLOP={}, H(vk,msg)={}, H(mu,w)={}, SampleInBall={}".format(
                r1cs["total_constraints"],
                r1cs["total_log2"],
                r1cs["bdlop_constraints"],
                r1cs["hash_vk_msg_constraints"],
                r1cs["hash_mu_w_constraints"],
                r1cs["sample_in_ball_constraints"],
            )
        )
        print(
            "  R1CS details: msg={} field elements, vk={} field elements, mu={} field elements, w={} field elements, BDLOP shape=({}, {}, {})".format(
                r1cs["msg_elements"],
                r1cs["vk_elements"],
                r1cs["mu_elements"],
                r1cs["w_elements"],
                r1cs["bdlop_rows_a"],
                r1cs["bdlop_message_len"],
                r1cs["bdlop_randomness_len"],
            )
        )
        print(
            "  R1CS hash details: H(vk,msg) input={}, perms={}, state/rate={}/{}; H(mu,w) input={}, perms={}, state/rate={}/{}".format(
                r1cs["hash_vk_msg_input_len"],
                r1cs["hash_vk_msg_permutations"],
                r1cs["hash_vk_msg_poseidon"]["state_size"],
                r1cs["hash_vk_msg_poseidon"]["rate"],
                r1cs["hash_mu_w_input_len"],
                r1cs["hash_mu_w_permutations"],
                r1cs["hash_mu_w_poseidon"]["state_size"],
                r1cs["hash_mu_w_poseidon"]["rate"],
            )
        )
    if not result.get("summary_only", False):
        outputs = result["estimator_outputs"]
        print("  raw lattice-estimator outputs:")
        for label, raw in (
            ("Verification-key MLWE", outputs["verification_key_mlwe"]),
            ("wprime MLWE", outputs["wprime_mlwe"]),
            ("Signature MSIS", outputs["signature_msis"]),
        ):
            print("    {}:".format(label))
            for line in raw.splitlines():
                print("      {}".format(line))
    print("")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Standalone Sage verification of the BRaccoon table parameters."
    )
    parser.add_argument("--estimator-path", default="lattice-estimator-main")
    parser.add_argument("--target", type=float, default=128.0)
    parser.add_argument("--full-estimator", action="store_true")
    parser.add_argument("--json-out", default=None)
    parser.add_argument("--no-r1cs", action="store_true")
    parser.add_argument("--r1cs-bdlop-randomness", type=int, default=13)
    parser.add_argument("--r1cs-bdlop-rows-a", type=int, default=2)
    parser.add_argument("--r1cs-bdlop-message-len", default="auto")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--r1cs-hash-states", default="32")
    parser.add_argument("--r1cs-hash-alpha", type=int, default=3)
    parser.add_argument("--r1cs-hash-capacity", default="auto")
    parser.add_argument("--r1cs-xof-state", type=int, default=11)
    parser.add_argument("--r1cs-xof-alpha", type=int, default=5)
    parser.add_argument("--r1cs-xof-capacity", default="auto")
    parser.add_argument("--r1cs-mu-elements", default="auto")
    parser.add_argument("--r1cs-w-elements", default="auto")
    parser.add_argument("--r1cs-msg-elements", type=int, default=256)
    parser.add_argument("--r1cs-vk-elements", default="auto")
    args = parser.parse_args(argv)
    args.r1cs_bdlop_message_len = parse_int_or_auto(args.r1cs_bdlop_message_len)
    args.r1cs_hash_states = parse_range_or_auto(args.r1cs_hash_states)
    args.r1cs_hash_capacity = parse_int_or_auto(args.r1cs_hash_capacity)
    args.r1cs_xof_capacity = parse_int_or_auto(args.r1cs_xof_capacity)
    args.r1cs_mu_elements = parse_int_or_auto(args.r1cs_mu_elements)
    args.r1cs_w_elements = parse_int_or_auto(args.r1cs_w_elements)
    args.r1cs_vk_elements = parse_int_or_auto(args.r1cs_vk_elements)
    return args


def main(argv=None):
    args = parse_args(argv)
    LWE, SIS, ND = load_estimator(args.estimator_path)

    print("Standalone BRaccoon parameter verification")
    print("target: {:.0f} bits".format(args.target))
    print("estimator: {}".format("full" if args.full_estimator else "rough"))
    print("estimator path: {}".format(resolve_path(args.estimator_path)))
    print("fixed q: {} (log2(q)={:.2f})".format(TABLE_Q, math.log2(TABLE_Q)))
    if args.no_r1cs:
        print("R1CS: skipped")
    else:
        print("R1CS: L_1,2 with in-circuit hashes H(vk,msg) and H(mu,w)")
        print("R1CS convention: msg has {} field elements, vk is a 2*lambda-bit seed".format(
            args.r1cs_msg_elements,
        ))
    print("")

    results = [verify_parameter_set(params, LWE, SIS, ND, args) for params in PARAMETER_SETS]
    for result in results:
        result["summary_only"] = args.summary_only
        print_result(result, args.target)

    all_pass = all(result["security"]["passes"] for result in results)
    print("overall:", "PASS" if all_pass else "FAIL")

    if args.json_out:
        out_path = Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
        print("wrote {}".format(out_path))

    return int(not all_pass)


sys.exit(int(main()))
