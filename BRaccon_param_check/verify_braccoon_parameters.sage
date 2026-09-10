#!/usr/bin/env sage
"""Fixed-parameter certificate for the current BRaccoon tables.

This verifier recomputes the published numerical quantities and runs the
bundled lattice-estimator on the associated fixed LWE/SIS instances.  It is
not a parameter search.
"""

import argparse
import json
import math
import pprint
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from sage.all import next_prime


SCRIPT_DIR = Path(__file__).resolve().parent
SIGNATURE_Q = 2**66 - 5
MIXED_R1CS_MODULUS = 2**256 + 1
KECCAK_CONSTRAINTS = 24 * 1600
LABRADOR_PROOF_KB = 110.0


@dataclass(frozen=True)
class ParameterSet:
    name: str
    max_signatures: int
    secret_dim: int
    expected_bgv_bits: int
    expected_ciphertext_kb: float
    expected_ext_commitment_kb: float
    expected_total_kb: float
    ring_degree: int = 256
    q_sig: int = SIGNATURE_Q
    module_dim: int = 12
    challenge_weight: int = 23
    nu_t: int = 52
    nu_w: int = 55
    sigma_t: int = 2**10
    sigma_wprime: int = 2**17
    sigma_w: int = 2**57
    bgv_rows: int = 10
    bgv_secret_dim: int = 18
    bgv_margin: int = 4
    bgv_flood_bound: int = 2**40
    ext_rows: int = 18
    ext_randomness: int = 15

    @property
    def q(self):
        """Compatibility alias used by compare_reduction_bounds.py."""
        return self.q_sig


PARAMETER_SETS = (
    ParameterSet("Q_s=2^20", 2**20, 14, 68, 23.89, 116.74, 720.43),
    ParameterSet("Q_s=2^32", 2**32, 16, 68, 23.91, 126.76, 778.57),
    ParameterSet("Q_s=2^64", 2**64, 22, 69, 23.95, 156.87, 953.06),
)


def resolve_path(path_text):
    path = Path(path_text)
    candidates = [path] if path.is_absolute() else [
        SCRIPT_DIR / path, Path.cwd() / path, SCRIPT_DIR.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def load_estimator(estimator_path):
    root = resolve_path(estimator_path)
    if not root.exists():
        raise FileNotFoundError("Estimator path not found: {}".format(root))
    sys.path.insert(0, str(root))
    from estimator import LWE, SIS, ND
    return LWE, SIS, ND


def first_prime_at_least(value):
    return int(next_prime(max(2, int(math.ceil(value))) - 1))


def min_log2_rop(result):
    if not result:
        return None
    bits = []
    for item in result.values():
        rop = item.get("rop", None)
        if rop is None:
            continue
        try:
            value = float(rop)
        except Exception:
            continue
        if value > 0:
            bits.append(math.log2(value))
    return min(bits) if bits else None


def estimate_lwe(LWE, ND, *, n, q, m, distribution, secret_width, error_width, tag, rough):
    Xs = ND.Binary if distribution == "binary" else ND.DiscreteGaussian(secret_width)
    Xe = ND.Binary if distribution == "binary" else ND.DiscreteGaussian(error_width)
    params = LWE.Parameters(n=n, q=q, Xs=Xs, Xe=Xe, m=m, tag=tag)
    fn = LWE.estimate.rough if rough else LWE.estimate
    raw = fn(params, quiet=True)
    return min_log2_rop(raw), raw


def estimate_sis(SIS, *, n, q, m, bound, tag, rough):
    params = SIS.Parameters(
        n=n, q=q, m=m, length_bound=math.ceil(bound), norm=2, tag=tag,
    )
    fn = SIS.estimate.rough if rough else SIS.estimate
    raw = fn(params, quiet=True)
    return min_log2_rop(raw), raw


def compute_signature_profile(p, target_bits):
    n, k, ell = p.ring_degree, p.module_dim, p.secret_dim
    omega, q = p.challenge_weight, p.q_sig
    root_all, root_k = math.sqrt(n * (k + ell)), math.sqrt(n * k)
    gaussian = math.exp(0.25)
    b_hmlwe = p.max_signatures * omega * (
        1.0 + n * (target_bits + 1.0 + 2.0 * math.log2(n))
        / math.sqrt(p.max_signatures)
    )
    sigma_mlwe = math.sqrt(
        1.0 / (2.0 / p.sigma_t**2 + 2.0 * b_hmlwe / p.sigma_wprime**2)
    )
    rounding = (omega * 2**p.nu_t + 2 ** (p.nu_w + 1)) * root_k
    b2_star = (
        gaussian * (omega * p.sigma_t + p.sigma_wprime + p.sigma_w) * root_all
        + rounding
    )
    # Current rounded-to-non-rounded reduction, with the exact remainder.
    b_corr = gaussian * (omega * p.sigma_t + p.sigma_wprime) * root_all
    b_w = gaussian * p.sigma_w * root_all
    delta_rounding = (
        omega * (2**p.nu_t - 1) + (2**p.nu_w - 1) + q % (2**p.nu_w)
    ) * root_k
    auxiliary_b = max(b2_star + delta_rounding, b_corr + b_w)
    b_msis = auxiliary_b + math.sqrt(omega)
    entropy_rhs = 2.0 * n * q ** (1.0 / (k + ell) + 2.0 / (n * ell))

    tail = 0.01
    while math.log(2.0) - tail**2 / 2.0 > math.log(0.01):
        tail += 0.01
    b_inf = tail * p.sigma_t * math.sqrt(n) + tail * p.sigma_wprime + tail * p.sigma_w
    signature_bits = (
        n + n * ell * math.ceil(math.log2(b_inf))
        + n * k * math.ceil(math.log2(b2_star / ((2**p.nu_w) * root_k)))
    )
    public_key_bits = n * k * (math.log2(q) - p.nu_t) + target_bits
    return {
        "B_HMLWE": b_hmlwe, "sigma_MLWE": sigma_mlwe,
        "B_2_star": b2_star, "B_corr": b_corr, "B_w": b_w,
        "delta_rounding": delta_rounding, "B_auxiliary": auxiliary_b,
        "B_MSIS": b_msis,
        "q_satisfies_msis": b_msis < (q - 1) / 2.0,
        "q_half_headroom_bits": math.log2((q - 1) / (2.0 * b_msis)),
        "entropy_rhs": entropy_rhs,
        "entropy_sigma_ok": p.sigma_wprime > entropy_rhs,
        "entropy_rounding_ok": p.nu_w < math.log2(q) - 2.0,
        "min_entropy_bits": n - 1,
        "signature_kb": signature_bits / 8000.0,
        "public_key_kb": public_key_bits / 8000.0,
    }


def compute_profile(p, q, target_bits):
    """Legacy-bound compatibility hook for compare_reduction_bounds.py.

    The main certificate uses compute_signature_profile and the current bound.
    The separate comparison script still needs the former certificate input as
    its baseline, so keep that calculation isolated here.
    """
    n, k, ell, omega = p.ring_degree, p.module_dim, p.secret_dim, p.challenge_weight
    root_all, root_k = math.sqrt(n * (k + ell)), math.sqrt(n * k)
    rounding = (omega * 2**p.nu_t + 2 ** (p.nu_w + 1)) * root_k
    b2_star = (
        math.exp(0.25) * (omega * p.sigma_t + p.sigma_wprime + p.sigma_w) * root_all
        + rounding
    )
    legacy_b_msis = b2_star + math.sqrt(omega) + rounding - omega
    return {"B_MSIS": legacy_b_msis}


def compute_bgv_profile(p):
    n = p.ring_degree
    dec_secret = math.sqrt(n * p.bgv_secret_dim)
    pk_error = enc_randomness = math.sqrt(n * p.bgv_rows)
    enc_e1, enc_e2 = math.sqrt(n * p.bgv_secret_dim), math.sqrt(n)
    enc_noise = pk_error * enc_randomness + enc_e2 + dec_secret * enc_e1
    b_t = math.exp(0.25) * p.sigma_t * math.sqrt(n * p.module_dim)
    b_wprime = math.exp(0.25) * p.sigma_wprime * math.sqrt(
        n * (p.module_dim + p.secret_dim)
    )
    b_plain = p.challenge_weight * b_t + b_wprime
    h = 2 * math.ceil(b_plain) + 1
    eval_noise = b_t * enc_noise + enc_noise
    final_noise = eval_noise + p.bgv_flood_bound
    correctness_lhs = b_plain + h * final_noise
    q_min = math.ceil(2.0 * correctness_lhs * p.bgv_margin)
    q_bgv = first_prime_at_least(q_min)
    ciphertext_kb = (p.bgv_rows + 1) * n * math.log2(q_bgv) / 8000.0
    return {
        "rows": p.bgv_rows, "secret_dim": p.bgv_secret_dim,
        "distribution": "binary", "B_plain": b_plain, "h": h,
        "h_ok": h > 2.0 * b_plain, "B_enc": enc_noise,
        "B_eval": eval_noise, "B_flood": p.bgv_flood_bound,
        "B_final": final_noise, "correctness_lhs": correctness_lhs,
        "q_min": q_min, "q": q_bgv, "q_bits": math.log2(q_bgv),
        "q_table_bits": math.ceil(math.log2(q_bgv)),
        "correctness_ok": correctness_lhs < q_bgv / 2.0,
        "ciphertext_kb": ciphertext_kb,
    }


def compute_extractable_commitment_profile(p, q_bgv):
    n, message_len = p.ring_degree, 2 * p.secret_dim + 1
    opening = math.sqrt(n * p.ext_randomness)
    extraction_noise = opening**2
    message_bound = math.exp(0.25) * p.sigma_w * math.sqrt(
        n * (p.module_dim + p.secret_dim)
    )
    h_ext = 2 * math.ceil(message_bound) + 1
    correctness_lhs = message_bound + h_ext * extraction_noise
    q_ext = first_prime_at_least(math.ceil(2.0 * correctness_lhs) + 1)
    commitment_kb = (
        (p.ext_rows + message_len) * n * math.log2(q_ext) / 8000.0
    )
    # All native-modulus equations are lifted independently into M.
    conservative_expression = max(q_ext * 2, q_bgv * 2, p.q_sig * 2)
    return {
        "rows": p.ext_rows, "message_len": message_len,
        "randomness": p.ext_randomness, "distribution": "binary",
        "opening_bound": opening, "binding_bound": opening,
        "message_bound": message_bound,
        "extraction_noise_bound": extraction_noise,
        "h_ext": h_ext, "h_ext_ok": h_ext > 2.0 * message_bound,
        "correctness_lhs": correctness_lhs, "q": q_ext,
        "q_bits": math.log2(q_ext),
        "q_table_bits": math.ceil(math.log2(q_ext)),
        # Compare integers here: converting the 78-bit prime back to binary64
        # can erase the few low bits by which it exceeds the lower bound.
        "extraction_ok": q_ext > math.ceil(2.0 * correctness_lhs) and math.gcd(h_ext, q_ext) == 1,
        "mixed_r1cs_modulus": MIXED_R1CS_MODULUS,
        "mixed_field_no_wrap_ok": 2.0 * conservative_expression < MIXED_R1CS_MODULUS,
        "commitment_kb": commitment_kb,
    }


def padded_blocks(bit_length, rate, suffix_bits):
    return int(math.ceil((bit_length + suffix_bits + 1) / float(rate)))


def compute_sha3_labrador_profile(p):
    q_w = math.ceil(p.q_sig / float(2**p.nu_w))
    w_bits = math.ceil(math.log2(q_w))
    serialized_w_bits = p.module_dim * p.ring_degree * w_bits
    digest_input_bits, challenge_input_bits = 512, 256 + serialized_w_bits
    digest_perms = padded_blocks(digest_input_bits, 1088, 2)
    challenge_perms = padded_blocks(challenge_input_bits, 1088, 2)
    hash_perms = digest_perms + challenge_perms
    hash_constraints = hash_perms * KECCAK_CONSTRAINTS
    shake_output_bits = p.challenge_weight * 160 + p.challenge_weight
    shake_perms = int(math.ceil(shake_output_bits / 1088.0))
    shake_constraints = shake_perms * KECCAK_CONSTRAINTS
    fisher_yates = p.challenge_weight * (3 * p.ring_degree + 4)
    nonlinear = hash_constraints + shake_constraints + fisher_yates
    return {
        "q_w": q_w, "w_coefficient_bits": w_bits,
        "serialized_w_bits": serialized_w_bits,
        "digest_input_bits": digest_input_bits,
        "challenge_input_bits": challenge_input_bits,
        "sha3_permutations": hash_perms, "sha3_constraints": hash_constraints,
        "shake_output_bits": shake_output_bits,
        "shake_permutations": shake_perms, "shake_constraints": shake_constraints,
        "fisher_yates_constraints": fisher_yates,
        "nonlinear_constraints": nonlinear,
        "pi1_kb": LABRADOR_PROOF_KB, "pi2_kb": LABRADOR_PROOF_KB,
        "proof_size_is_analytical": True,
    }


def close_to(value, expected, tolerance=0.015):
    return abs(value - expected) <= tolerance


def verify_parameter_set(p, LWE, SIS, ND, args):
    rough = not args.full_estimator
    sig, bgv = compute_signature_profile(p, args.target), compute_bgv_profile(p)
    ext = compute_extractable_commitment_profile(p, bgv["q"])
    zk = compute_sha3_labrador_profile(p)

    vk_bits, vk_raw = estimate_lwe(
        LWE, ND, n=p.ring_degree * p.secret_dim, q=p.q_sig,
        m=p.ring_degree * p.module_dim, distribution="gaussian",
        secret_width=sig["sigma_MLWE"], error_width=sig["sigma_MLWE"],
        tag=p.name + " verification-key MLWE", rough=rough,
    )
    wprime_bits, wprime_raw = estimate_lwe(
        LWE, ND, n=p.ring_degree * p.secret_dim, q=p.q_sig,
        m=p.ring_degree * p.module_dim, distribution="gaussian",
        secret_width=p.sigma_wprime, error_width=p.sigma_w,
        tag=p.name + " wprime MLWE", rough=rough,
    )
    sig_bits, sig_raw = estimate_sis(
        SIS, n=p.ring_degree * p.module_dim, q=p.q_sig,
        m=p.ring_degree * (p.module_dim + p.secret_dim + 1),
        bound=sig["B_MSIS"], tag=p.name + " signature MSIS", rough=rough,
    )
    bgv_bits, bgv_raw = estimate_lwe(
        LWE, ND, n=p.ring_degree * p.bgv_secret_dim, q=bgv["q"],
        m=p.ring_degree * p.bgv_rows, distribution="binary",
        secret_width=0.5, error_width=0.5,
        tag=p.name + " BGV IND-CPA MLWE", rough=rough,
    )
    ext_binding_bits, ext_binding_raw = estimate_sis(
        SIS, n=p.ring_degree * p.ext_rows, q=ext["q"],
        m=p.ring_degree * p.ext_randomness, bound=ext["binding_bound"],
        tag=p.name + " extractable commitment binding MSIS", rough=rough,
    )
    ext_hiding_bits, ext_hiding_raw = estimate_lwe(
        LWE, ND, n=p.ring_degree * p.ext_randomness, q=ext["q"],
        m=p.ring_degree * (p.ext_rows + ext["message_len"]), distribution="binary",
        secret_width=0.5, error_width=0.5,
        tag=p.name + " extractable commitment hiding MLWE", rough=rough,
    )
    ext_setup_raw_bits, ext_setup_raw = estimate_lwe(
        LWE, ND, n=p.ring_degree * p.ext_rows, q=ext["q"],
        m=p.ring_degree * p.ext_randomness, distribution="binary",
        secret_width=0.5, error_width=0.5,
        tag=p.name + " extractable commitment setup MLWE", rough=rough,
    )
    ext_setup_bits = (
        None if ext_setup_raw_bits is None
        else ext_setup_raw_bits - math.log2(ext["message_len"])
    )

    ciphertext_count = p.secret_dim + 1
    # The paper multiplies the component sizes after rounding them to two
    # decimals (e.g. 17*20.02), so reproduce that convention exactly.
    communication_total = (
        zk["pi1_kb"] + zk["pi2_kb"]
        + ciphertext_count * round(bgv["ciphertext_kb"], 2)
        + round(ext["commitment_kb"], 2) + 25.34
    )
    expected_sig = {14: 28.00, 16: 31.78, 22: 43.10}[p.secret_dim]
    checks = {
        "signature_modulus_66_bits": math.ceil(math.log2(p.q_sig)) == 66,
        "signature_msis_modulus": sig["q_satisfies_msis"],
        "entropy_sigma": sig["entropy_sigma_ok"],
        "entropy_rounding": sig["entropy_rounding_ok"],
        "public_key_size": close_to(sig["public_key_kb"], 5.39),
        "signature_size": close_to(sig["signature_kb"], expected_sig),
        "bgv_correctness": bgv["correctness_ok"],
        "bgv_table_modulus": bgv["q_table_bits"] == p.expected_bgv_bits,
        "ciphertext_size": close_to(bgv["ciphertext_kb"], p.expected_ciphertext_kb),
        "extractable_commitment_correctness": ext["extraction_ok"],
        "extractable_commitment_modulus": ext["q_table_bits"] == 78,
        "mixed_r1cs_no_wrap": ext["mixed_field_no_wrap_ok"],
        "extractable_commitment_size": close_to(
            ext["commitment_kb"], p.expected_ext_commitment_kb
        ),
        "sha3_constraint_count": zk["nonlinear_constraints"] == 1438556,
        "communication_total": close_to(
            communication_total, p.expected_total_kb, tolerance=0.06
        ),
    }
    security = {
        "verification_key_mlwe": vk_bits, "wprime_mlwe": wprime_bits,
        "signature_msis": sig_bits, "bgv_ind_cpa_mlwe": bgv_bits,
        "ext_binding_msis": ext_binding_bits,
        "ext_hiding_mlwe": ext_hiding_bits,
        "ext_setup_mlwe_after_union_bound": ext_setup_bits,
    }
    present = [value for value in security.values() if value is not None]
    security_minimum = min(present)
    security_passes = all(value >= args.target for value in present)
    return {
        "parameters": asdict(p), "signature": sig, "bgv": bgv,
        "extractable_commitment": ext, "sha3_labrador": zk,
        "communication": {
            "ciphertext_count": ciphertext_count,
            "ciphertexts_kb": ciphertext_count * bgv["ciphertext_kb"],
            "wprime_kb": 25.34, "total_kb": communication_total,
        },
        "checks": checks,
        "security": dict(security, minimum=security_minimum, passes=security_passes),
        "passes": all(checks.values()) and security_passes,
        "estimator_outputs": {
            "verification_key_mlwe": vk_raw, "wprime_mlwe": wprime_raw,
            "signature_msis": sig_raw, "bgv_ind_cpa_mlwe": bgv_raw,
            "ext_binding_msis": ext_binding_raw,
            "ext_hiding_mlwe": ext_hiding_raw, "ext_setup_mlwe": ext_setup_raw,
        },
    }


def fmt_bits(value):
    if value is None:
        return "n/a"
    if math.isinf(value):
        return "inf"
    return "{:.2f}".format(value)


def print_result(result, target, summary_only):
    p, sig = result["parameters"], result["signature"]
    bgv, ext = result["bgv"], result["extractable_commitment"]
    zk, comm, sec = result["sha3_labrador"], result["communication"], result["security"]
    print("{} [{}]".format(p["name"], "PASS" if result["passes"] else "FAIL"))
    print("  signature: log2(q_sig)={:.2f}, k={}, ell={}, B_MSIS=2^{:.2f}, headroom={:.4f} bits".format(
        math.log2(p["q_sig"]), p["module_dim"], p["secret_dim"],
        math.log2(sig["B_MSIS"]), sig["q_half_headroom_bits"],
    ))
    print("  entropy: sigma_w'/required={:.2f}, nu_w={} < {:.2f}, H_inf >= {}".format(
        p["sigma_wprime"] / sig["entropy_rhs"], p["nu_w"],
        math.log2(p["q_sig"]) - 2, sig["min_entropy_bits"],
    ))
    print("  BGV: shape=({},{}), log2(Q)={:.2f} (table {}), ct={:.2f} KB, correctness={}".format(
        p["bgv_rows"], p["bgv_secret_dim"], bgv["q_bits"],
        bgv["q_table_bits"], bgv["ciphertext_kb"],
        "PASS" if bgv["correctness_ok"] else "FAIL",
    ))
    print("  extractable commitment: shape=({}, {}, {}), log2(Q_ext)={:.2f}, size={:.2f} KB, extraction={}".format(
        p["ext_rows"], ext["message_len"], p["ext_randomness"],
        ext["q_bits"], ext["commitment_kb"],
        "PASS" if ext["extraction_ok"] else "FAIL",
    ))
    print("  SHA3/SHAKE: permutations={}/{}, nonlinear constraints={} (2^{:.2f})".format(
        zk["sha3_permutations"], zk["shake_permutations"],
        zk["nonlinear_constraints"], math.log2(zk["nonlinear_constraints"]),
    ))
    print("  sizes: vk={:.2f} KB, sig={:.2f} KB, communication={:.2f} KB".format(
        sig["public_key_kb"], sig["signature_kb"], comm["total_kb"],
    ))
    print("  security: VK={}, w'={}, SIG={}, BGV={}, ext-bind={}, ext-hide={}, ext-setup={}, min={} (target={:.0f})".format(
        fmt_bits(sec["verification_key_mlwe"]), fmt_bits(sec["wprime_mlwe"]),
        fmt_bits(sec["signature_msis"]), fmt_bits(sec["bgv_ind_cpa_mlwe"]),
        fmt_bits(sec["ext_binding_msis"]), fmt_bits(sec["ext_hiding_mlwe"]),
        fmt_bits(sec["ext_setup_mlwe_after_union_bound"]),
        fmt_bits(sec["minimum"]), target,
    ))
    failed = [name for name, ok in result["checks"].items() if not ok]
    print("  numerical checks:", "PASS" if not failed else "FAIL: " + ", ".join(failed))
    if not summary_only:
        print("  raw lattice-estimator outputs:")
        for label, raw in result["estimator_outputs"].items():
            print("    {}: {}".format(label, pprint.pformat(raw, width=100)))
    print("")


def jsonable(value):
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return pprint.pformat(value, width=100)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimator-path", default="lattice-estimator-main")
    parser.add_argument("--target", type=float, default=128.0)
    parser.add_argument("--full-estimator", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--json-out", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    LWE, SIS, ND = load_estimator(args.estimator_path)
    print("BRaccoon fixed-parameter certificate (current paper tables)")
    print("target: {:.0f} bits".format(args.target))
    print("estimator: {}".format("full" if args.full_estimator else "rough"))
    print("estimator path: {}".format(resolve_path(args.estimator_path)))
    print("q_sig: {} (log2={:.2f})".format(SIGNATURE_Q, math.log2(SIGNATURE_Q)))
    print("mixed-R1CS field: 2^256+1")
    print("proof-size convention: 110 KB for each of pi_1 and pi_2 (analytical input)")
    print("BGV numerical profile: binary noise, B_flood=2^40")
    print("")
    results = [verify_parameter_set(p, LWE, SIS, ND, args) for p in PARAMETER_SETS]
    for result in results:
        print_result(result, args.target, args.summary_only)
    all_pass = all(result["passes"] for result in results)
    print("overall numerical certificate:", "PASS" if all_pass else "FAIL")
    print("analytical caveats: LaBRADOR proof size and the full-coin flooding distribution are external assumptions")
    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(jsonable(results), indent=2, sort_keys=True), encoding="utf-8")
        print("wrote {}".format(out))
    return int(not all_pass)


if __name__ == "__main__":
    sys.exit(main())
