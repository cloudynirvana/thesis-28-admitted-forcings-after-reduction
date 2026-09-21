#!/usr/bin/env python3
"""Admission ledger after a rank-preserving stiff-sloppy reduction.

Thesis #19 supplies the gate: a protocol-constant forcing may be admitted
only as a declared input, and a write into kinetic Theta is refused.
Thesis #24 supplies the reduction: slave the fast modifier, keep the
product kappa, and require the shared multi-channel practical rank to
agree on the two parameter vectors.

This script does not read either deposit's results file. The catalogue,
the rates, and the protocols are declared here. The steady map is closed
form. No random draw is used.

Research only. Not a medical device. Not a dose.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
CANDIDATE_PATH = ROOT / "candidates.json"
# SHA-256 of sim/candidates.json as stored. A byte edit raises.
CANDIDATE_PIN = "e4bcbe1791a65e3e884be37290b1eba4da33a4fb851aa374ef5c6275d5457f3e"
SEED = 20260921  # deterministic ledger; no RNG draws

FORCING_SYMBOLS = ("u_h", "u_p", "u_q", "u_w")
THETA_NAMES = ("a", "b", "d_p", "d_q", "h", "k_branch", "k_s", "k_side")
RED_NAMES = ("d_p", "d_q", "k_branch", "k_s", "k_side", "kappa")
FULL_ORDER = ("k_s", "k_branch", "k_side", "d_p", "d_q", "a", "b", "h")

# Decimal strings are the hashed object.
THETA = {
    "a": "1.25",
    "b": "5.00",
    "d_p": "0.45",
    "d_q": "0.35",
    "h": "0.40",
    "k_branch": "0.55",
    "k_s": "0.90",
    "k_side": "0.30",
}

# Synthetic design inputs. Not a cell-line panel.
PROTOCOLS = (
    ("P1", 0.80, 0.25),
    ("P2", 1.10, 0.40),
    ("P3", 1.40, 0.30),
    ("P4", 0.95, 0.55),
)

PROTOCOL_AMPLITUDE = "1"
HORIZON = "40"
SIGMA = 0.08
BOX = (0.25, 4.0)
PRACTICAL_CUT = 1e-3
NUMERICAL_CUT = 1e-8
FD_EPS = 1e-6
GRID_N = 401
# Equal-variance normal-normal witness on the log of h. Not the decision rule.
PRIOR_LOG_SLOPE = 0.25
Y0 = np.array([1.0, 0.50, 0.40, 0.60, 0.08, 0.08], dtype=float)
T_CHECK = 80.0

N_SHARED = 4 * 4
N_READOUT = 4 * 5
# One leaked coordinate: the profile budget is the 95th percentile on 1 degree of freedom.
BUDGET_PROFILE = float(chi2.ppf(0.95, 1))
# Schedule-wide balls. Reported, not used to name the outcome.
BUDGET_SHARED = float(chi2.ppf(0.95, N_SHARED))
BUDGET_READOUT = float(chi2.ppf(0.95, N_READOUT))


def canon_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def digest(obj) -> str:
    return hashlib.sha256(canon_bytes(obj)).hexdigest()


def load_candidates() -> dict:
    raw = CANDIDATE_PATH.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != CANDIDATE_PIN:
        raise SystemExit(f"candidate pin mismatch: {got}")
    return json.loads(raw.decode("utf-8"))


def theta_float(table: dict | None = None) -> dict[str, float]:
    src = THETA if table is None else table
    return {name: float(src[name]) for name in THETA_NAMES}


def kappa_of(th: dict[str, float], u_h: float = 0.0) -> float:
    return (th["h"] + u_h) * th["a"] / th["b"]


def rho_of(th: dict[str, float]) -> float:
    return th["a"] / th["b"]


def steady_full(
    th: dict[str, float],
    v_s: float,
    v_q: float,
    u_p: float = 0.0,
    u_q: float = 0.0,
    u_h: float = 0.0,
    u_w: float = 0.0,
) -> np.ndarray | None:
    """Steady (S, I, P, Q, Z, W). None if a coordinate leaves the positive orthant."""
    kappa = kappa_of(th, u_h)
    rho = rho_of(th)
    quad = _quadratic_intermediate(th["k_branch"], th["k_side"], kappa, v_s, th["k_s"])
    if quad is None:
        return None
    s_star, i_star, branch = quad
    p_star = (branch + u_p) / th["d_p"]
    q_star = (v_q + th["k_side"] * i_star + u_q) / th["d_q"]
    z_star = rho * i_star
    w_star = rho * i_star + u_w / th["b"]
    out = np.array([s_star, i_star, p_star, q_star, z_star, w_star], dtype=float)
    if not np.all(np.isfinite(out)) or np.any(out <= 0.0):
        return None
    return out


def steady_reduced(
    phi: dict[str, float],
    v_s: float,
    v_q: float,
    u_p: float = 0.0,
    u_q: float = 0.0,
) -> np.ndarray | None:
    """Steady (S, I, P, Q) of the reduced field. No modifier forcing, no readout."""
    quad = _quadratic_intermediate(phi["k_branch"], phi["k_side"], phi["kappa"], v_s, phi["k_s"])
    if quad is None:
        return None
    s_star, i_star, branch = quad
    p_star = (branch + u_p) / phi["d_p"]
    q_star = (v_q + phi["k_side"] * i_star + u_q) / phi["d_q"]
    out = np.array([s_star, i_star, p_star, q_star], dtype=float)
    if not np.all(np.isfinite(out)) or np.any(out <= 0.0):
        return None
    return out


def _quadratic_intermediate(
    k_branch: float, k_side: float, kappa: float, v_s: float, k_s: float
) -> tuple[float, float, float] | None:
    if min(k_branch, k_side, kappa, v_s, k_s) <= 0.0:
        return None
    a_coef = k_branch * kappa
    b_coef = k_branch + k_side
    disc = b_coef * b_coef + 4.0 * a_coef * v_s
    if disc <= 0.0:
        return None
    i_star = (-b_coef + math.sqrt(disc)) / (2.0 * a_coef)
    if i_star <= 0.0:
        return None
    s_star = v_s / k_s
    branch = k_branch * (1.0 + kappa * i_star) * i_star
    if s_star <= 0.0 or branch <= 0.0:
        return None
    return s_star, i_star, branch


def phi_from_theta(th: dict[str, float]) -> dict[str, float]:
    return {
        "d_p": th["d_p"],
        "d_q": th["d_q"],
        "k_branch": th["k_branch"],
        "k_s": th["k_s"],
        "k_side": th["k_side"],
        "kappa": kappa_of(th, 0.0),
    }


def scale_table(table: dict[str, float], name: str, factor: float) -> dict[str, float]:
    out = dict(table)
    out[name] = table[name] * factor
    return out


def stack_shared(rows: list[np.ndarray]) -> np.ndarray:
    parts = [np.log(row[:4]) for row in rows]
    return np.concatenate(parts)


def stack_readout(rows: list[np.ndarray]) -> np.ndarray:
    """Shared logs plus log W. Z is not in this schedule."""
    parts = []
    for row in rows:
        parts.append(np.log(np.array([row[0], row[1], row[2], row[3], row[5]], dtype=float)))
    return np.concatenate(parts)


def chi2_distance(left: np.ndarray, right: np.ndarray) -> float:
    delta = (left - right) / SIGMA
    return float(np.dot(delta, delta))


def full_rows(th: dict[str, float], forcing: dict[str, float]) -> list[np.ndarray] | None:
    rows = []
    for _name, v_s, v_q in PROTOCOLS:
        state = steady_full(th, v_s, v_q, **forcing)
        if state is None:
            return None
        rows.append(state)
    return rows


def reduced_rows(phi: dict[str, float], u_p: float = 0.0, u_q: float = 0.0) -> list[np.ndarray] | None:
    rows = []
    for _name, v_s, v_q in PROTOCOLS:
        state = steady_reduced(phi, v_s, v_q, u_p=u_p, u_q=u_q)
        if state is None:
            return None
        rows.append(state)
    return rows


def null_forcing() -> dict[str, float]:
    return {name: 0.0 for name in FORCING_SYMBOLS}


def one_forcing(symbol: str, amplitude: float) -> dict[str, float]:
    forcing = null_forcing()
    forcing[symbol] = amplitude
    return forcing


def search_scale(evaluate) -> dict:
    """Minimum of evaluate(factor) on the closed factor interval."""
    lo, hi = BOX
    best_factor = 1.0
    best_chi = math.inf

    def consider(factor: float) -> None:
        nonlocal best_factor, best_chi
        value = evaluate(factor)
        if value is None or not math.isfinite(value):
            return
        if value < best_chi:
            best_factor = float(factor)
            best_chi = float(value)

    consider(1.0)
    for factor in np.exp(np.linspace(math.log(lo), math.log(hi), GRID_N)):
        consider(float(factor))

    def objective(log_factor: float) -> float:
        value = evaluate(math.exp(log_factor))
        if value is None or not math.isfinite(value):
            return 1e12
        return float(value)

    solved = minimize_scalar(
        objective,
        bounds=(math.log(lo), math.log(hi)),
        method="bounded",
        options={"xatol": 1e-14},
    )
    if solved.success:
        consider(math.exp(float(solved.x)))
    if not math.isfinite(best_chi):
        raise SystemExit("leakage search did not return a finite chi-square")
    return {"factor": best_factor, "chi2": best_chi}


def min_leakage_full(y_star: np.ndarray, baseline: dict[str, float], schedule: str) -> dict:
    found = []
    for name in FULL_ORDER:
        def evaluate(factor: float, name: str = name) -> float | None:
            moved = scale_table(baseline, name, factor)
            rows = full_rows(moved, null_forcing())
            if rows is None:
                return None
            obs = stack_readout(rows) if schedule == "readout" else stack_shared(rows)
            return chi2_distance(y_star, obs)

        hit = search_scale(evaluate)
        hit["coordinate"] = name
        found.append(hit)
    found.sort(key=lambda item: item["chi2"])
    return found[0] | {"runners": [{k: item[k] for k in ("coordinate", "factor", "chi2")} for item in found]}


def min_leakage_reduced(y_star: np.ndarray, phi: dict[str, float]) -> dict:
    found = []
    for name in RED_NAMES:
        def evaluate(factor: float, name: str = name) -> float | None:
            moved = scale_table(phi, name, factor)
            rows = reduced_rows(moved)
            if rows is None:
                return None
            return chi2_distance(y_star, stack_shared(rows))

        hit = search_scale(evaluate)
        hit["coordinate"] = name
        found.append(hit)
    found.sort(key=lambda item: item["chi2"])
    return found[0] | {"runners": [{k: item[k] for k in ("coordinate", "factor", "chi2")} for item in found]}


def fisher_rank(predict, names: list[str], table: dict[str, float]) -> dict:
    centre = predict(table)
    if centre is None:
        raise SystemExit("fisher centre is undefined")
    columns = []
    for name in names:
        up = scale_table(table, name, math.exp(FD_EPS))
        down = scale_table(table, name, math.exp(-FD_EPS))
        yup = predict(up)
        ydown = predict(down)
        if yup is None or ydown is None:
            raise SystemExit(f"fisher step left the domain at {name}")
        columns.append((yup - ydown) / (2.0 * FD_EPS))
    jac = np.column_stack(columns) / SIGMA
    fim = jac.T @ jac
    fim = 0.5 * (fim + fim.T)
    eigenvalues = np.linalg.eigvalsh(fim)
    eigenvalues = np.array(sorted((float(v) for v in eigenvalues), reverse=True))
    leading = eigenvalues[0]
    if leading <= 0.0:
        raise SystemExit("fisher leading eigenvalue is not positive")
    ratios = eigenvalues / leading
    return {
        "names": list(names),
        "eigenvalues": eigenvalues.tolist(),
        "ratios": ratios.tolist(),
        "practical_rank": int(np.sum(ratios > PRACTICAL_CUT)),
        "numerical_rank": int(np.sum(ratios > NUMERICAL_CUT)),
        "dimension": len(names),
    }


def unreduced_plane_mass(eigvecs_smallest: np.ndarray) -> list[float]:
    """Squared mass of vectors on span{e_rho, e_speed} inside the (a, b, h) block."""
    e_kappa = np.array([1.0, -1.0, 1.0]) / math.sqrt(3.0)
    e_speed = np.array([1.0, 1.0, 0.0]) / math.sqrt(2.0)
    e_rho = np.array([1.0, -1.0, -2.0]) / math.sqrt(6.0)
    # Orthonormal check is part of the run, not a result copied from elsewhere.
    frame = np.column_stack([e_kappa, e_speed, e_rho])
    gram = frame.T @ frame
    if np.max(np.abs(gram - np.eye(3))) > 1e-12:
        raise SystemExit("modifier frame is not orthonormal")
    plane = np.column_stack([e_speed, e_rho])
    masses = []
    for vec in eigvecs_smallest:
        block = vec[-3:]
        masses.append(float(np.dot(block, plane @ (plane.T @ block))))
    return masses


def full_fisher_with_plane(predict, table: dict[str, float]) -> dict:
    centre = predict(table)
    columns = []
    for name in FULL_ORDER:
        up = scale_table(table, name, math.exp(FD_EPS))
        down = scale_table(table, name, math.exp(-FD_EPS))
        columns.append((predict(up) - predict(down)) / (2.0 * FD_EPS))
    jac = np.column_stack(columns) / SIGMA
    fim = 0.5 * ((jac.T @ jac) + (jac.T @ jac).T)
    eigenvalues, eigenvectors = np.linalg.eigh(fim)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    leading = float(eigenvalues[0])
    ratios = eigenvalues / leading
    smallest = [eigenvectors[:, -1], eigenvectors[:, -2]]
    packed = fisher_rank(predict, list(FULL_ORDER), table)
    packed["smallest_unreduced_mass"] = unreduced_plane_mass(smallest)
    packed["ratios"] = [float(v) for v in ratios]
    packed["eigenvalues"] = [float(v) for v in eigenvalues]
    return packed


def rhs(t: float, y: np.ndarray, th: dict[str, float], v_s: float, v_q: float, forcing: dict[str, float]) -> np.ndarray:
    s_state, i_state, _p, _q, z_state, _w = y
    h_eff = th["h"] + forcing["u_h"]
    branch = th["k_branch"] * (1.0 + h_eff * z_state)
    return np.array(
        [
            v_s - th["k_s"] * s_state,
            th["k_s"] * s_state - (branch + th["k_side"]) * i_state,
            branch * i_state - th["d_p"] * y[2] + forcing["u_p"],
            v_q + th["k_side"] * i_state - th["d_q"] * y[3] + forcing["u_q"],
            th["a"] * i_state - th["b"] * z_state,
            th["a"] * i_state - th["b"] * y[5] + forcing["u_w"],
        ],
        dtype=float,
    )


def integration_gap(th: dict[str, float], forcing: dict[str, float]) -> float:
    _name, v_s, v_q = PROTOCOLS[0]
    target = steady_full(th, v_s, v_q, **forcing)
    if target is None:
        raise SystemExit("integration target is undefined")
    sol = solve_ivp(
        rhs,
        (0.0, T_CHECK),
        Y0,
        t_eval=np.array([T_CHECK]),
        args=(th, v_s, v_q, forcing),
        method="LSODA",
        rtol=1e-8,
        atol=1e-11,
    )
    if not sol.success:
        raise SystemExit(f"integration failed: {sol.message}")
    return float(np.max(np.abs(sol.y[:, 0] - target)))


def algebraic_gap(state: np.ndarray, th: dict[str, float], v_s: float, v_q: float, forcing: dict[str, float]) -> float:
    return float(np.max(np.abs(rhs(0.0, state, th, v_s, v_q, forcing))))


def null_schedule(symbol: str) -> dict:
    return {
        "amplitude_rule": "null",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": "0"}],
        "provenance": None,
        "symbol": symbol,
    }


def legal_schedule(symbol: str, row: dict, library_sha: str) -> dict:
    return {
        "amplitude_rule": "protocol_constant",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": PROTOCOL_AMPLITUDE}],
        "provenance": {
            "claim": int(row["claim"]),
            "library_sha256": library_sha,
            "pains": int(row["pains"]),
            "row_id": row["id"],
            "score_copied": False,
            "symbol": symbol,
        },
        "symbol": symbol,
    }


def row_by_id(doc: dict) -> dict:
    return {row["id"]: row for row in doc["rows"]}


def failed_checks(doc: dict, call: dict, occupied: set[str]) -> list[str]:
    reasons: list[str] = []
    kind = call["kind"]
    destination = call["destination"]
    rule = call["amplitude_rule"]
    if kind == "anecdote":
        reasons.append("not_a_score_record")
    if kind == "utility":
        reasons.append("utility_not_an_input")
    if rule == "soft_prior" or kind == "soft_prior":
        reasons.append("soft_prior")
    if rule == "soft_weight" or kind == "soft_weight":
        reasons.append("soft_weight")
    if destination in THETA_NAMES:
        reasons.append("theta_destination")
    elif destination not in FORCING_SYMBOLS:
        reasons.append("symbol")
    if kind == "anecdote":
        return reasons

    row = row_by_id(doc)[call["row_id"]]
    library = call.get("library_sha256", doc["library_sha256"])
    engine = call.get("engine", doc["engine"])
    unit = call.get("unit", doc["unit"])
    may = call.get("may_enter_theta", doc["may_enter_theta"])
    claim = call.get("claim", row["claim"])
    pains = call.get("pains", row["pains"])
    if library != doc["library_sha256"]:
        reasons.append("library_sha")
    if engine is not None:
        reasons.append("engine")
    if unit != "arbitrary_surrogate":
        reasons.append("unit")
    if may is not False:
        reasons.append("may_enter_theta")
    if int(claim) != 3:
        reasons.append("claim")
    if int(pains) != 0:
        reasons.append("pains")
    if rule != "protocol_constant":
        reasons.append("amplitude_rule")
    if not reasons and destination in occupied:
        reasons.append("slot_occupied")
    return reasons


def apply_call(doc: dict, ledger: dict, call: dict) -> dict:
    occupied = {
        symbol
        for symbol, schedule in ledger["forcing"].items()
        if schedule["amplitude_rule"] != "null"
    }
    reasons = failed_checks(doc, call, occupied)
    before_theta = digest(ledger["theta"])
    before_forcing = digest(ledger["forcing"])
    record = {
        "id": call["id"],
        "kind": call["kind"],
        "row_id": call.get("row_id"),
        "destination": call["destination"],
        "amplitude_rule": call["amplitude_rule"],
        "reasons": reasons,
        "status": "refused" if reasons else "admitted",
        "theta_sha256": before_theta,
    }
    if record["status"] == "admitted":
        row = row_by_id(doc)[call["row_id"]]
        symbol = call["destination"]
        ledger["forcing"][symbol] = legal_schedule(symbol, row, doc["library_sha256"])
        record["admitted_row"] = row["id"]
        record["admitted_symbol"] = symbol
        if "proxy_score" in ledger["forcing"][symbol]["provenance"]:
            raise SystemExit("provenance recorded a proxy score")
    after_theta = digest(ledger["theta"])
    record["theta_sha256_after"] = after_theta
    record["forcing_sha256_after"] = digest(ledger["forcing"])
    record["theta_unchanged"] = before_theta == after_theta
    record["forcing_unchanged"] = before_forcing == record["forcing_sha256_after"]
    if not record["theta_unchanged"]:
        raise SystemExit(f"{call['id']} changed Theta")
    if record["status"] == "refused" and not record["forcing_unchanged"]:
        raise SystemExit(f"{call['id']} refusal changed the forcing")
    if record["status"] == "admitted" and record["forcing_unchanged"]:
        raise SystemExit(f"{call['id']} admission did not write a forcing")
    ledger["calls"].append(record)
    return record


def build_calls(doc: dict) -> list[dict]:
    library = doc["library_sha256"]
    flipped = library[:-1] + ("0" if library[-1] != "0" else "1")
    return [
        {"id": "R01", "kind": "score", "row_id": "C06", "destination": "u_p", "amplitude_rule": "protocol_constant"},
        {"id": "R02", "kind": "score", "row_id": "C07", "destination": "u_q", "amplitude_rule": "protocol_constant"},
        {"id": "R03", "kind": "score", "row_id": "C08", "destination": "u_h", "amplitude_rule": "protocol_constant"},
        {"id": "R04", "kind": "soft_prior", "row_id": "C01", "destination": "d_p", "amplitude_rule": "soft_prior"},
        {"id": "R05", "kind": "soft_weight", "row_id": "C02", "destination": "k_side", "amplitude_rule": "soft_weight"},
        {"id": "R06", "kind": "utility", "row_id": "C01", "destination": "k_s", "amplitude_rule": "protocol_constant"},
        {"id": "R07", "kind": "score", "row_id": "C01", "destination": "u_p", "amplitude_rule": "score_copy"},
        {"id": "R08", "kind": "anecdote", "destination": "u_p", "amplitude_rule": "protocol_constant"},
        {"id": "R09", "kind": "score", "row_id": "C01", "destination": "u_p", "amplitude_rule": "protocol_constant", "library_sha256": flipped},
        {"id": "R10", "kind": "score", "row_id": "C01", "destination": "u_p", "amplitude_rule": "protocol_constant", "engine": "vina"},
        {"id": "R11", "kind": "score", "row_id": "C01", "destination": "u_p", "amplitude_rule": "protocol_constant", "may_enter_theta": True},
        {"id": "R12", "kind": "score", "row_id": "C01", "destination": "N_ROS", "amplitude_rule": "protocol_constant"},
        {"id": "R13", "kind": "score", "row_id": "C01", "destination": "h", "amplitude_rule": "protocol_constant"},
        {"id": "A01", "kind": "score", "row_id": "C01", "destination": "u_p", "amplitude_rule": "protocol_constant"},
        {"id": "A02", "kind": "score", "row_id": "C02", "destination": "u_q", "amplitude_rule": "protocol_constant"},
        {"id": "A03", "kind": "score", "row_id": "C03", "destination": "u_h", "amplitude_rule": "protocol_constant"},
        {"id": "A04", "kind": "score", "row_id": "C04", "destination": "u_w", "amplitude_rule": "protocol_constant"},
        {"id": "R14", "kind": "score", "row_id": "C05", "destination": "u_p", "amplitude_rule": "protocol_constant"},
        {"id": "R15", "kind": "score", "row_id": "C03", "destination": "a", "amplitude_rule": "protocol_constant"},
    ]


def eligibility(doc: dict) -> list[dict]:
    rows = []
    for row in doc["rows"]:
        call = {
            "id": "elig_" + row["id"],
            "kind": "score",
            "row_id": row["id"],
            "destination": row["symbol"],
            "amplitude_rule": "protocol_constant",
        }
        reasons = failed_checks(doc, call, occupied=set())
        rows.append(
            {
                "id": row["id"],
                "label": row["label"],
                "symbol": row["symbol"],
                "claim": row["claim"],
                "pains": row["pains"],
                "eligible": not reasons,
                "reasons": reasons,
            }
        )
    return rows


def classify(d_reduced: float, d_full: float) -> str:
    """Name the admission from the one-coordinate profile budget.

    Survives: the reduced chart still rejects every single-coordinate leakage.
    Artefact: only the chart that still contains W rejects them.
    Confounded: both charts have a leakage inside the budget.
    """
    survives = d_reduced > BUDGET_PROFILE
    full_separates = d_full > BUDGET_PROFILE
    if survives:
        return "survives"
    if full_separates:
        return "artefact_unreduced"
    return "confounded"


def soft_prior_witness(doc: dict, baseline_sha: str) -> dict:
    row = row_by_id(doc)["C03"]
    score = float(row["proxy_score"])
    h0 = float(THETA["h"])
    log_shift = PRIOR_LOG_SLOPE * (-score)
    posterior = h0 * math.exp(0.5 * log_shift)
    fictional = dict(THETA)
    fictional["h"] = format(posterior, ".12g")
    fictional_sha = digest(fictional)
    if fictional_sha == baseline_sha:
        raise SystemExit("soft-prior witness collided with the kinetic digest")
    return {
        "row_id": "C03",
        "coordinate": "h",
        "proxy_score": row["proxy_score"],
        "log_slope": PRIOR_LOG_SLOPE,
        "posterior_h": fictional["h"],
        "fictional_theta_sha256": fictional_sha,
        "written": False,
    }


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def plot_spectra(spec_full: dict, spec_red: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    full = np.array(spec_full["ratios"])
    red = np.array(spec_red["ratios"])
    ax.plot(np.arange(1, full.size + 1), np.log10(np.maximum(full, 1e-18)), "o-", color="#0072B2", label="Full, 8 rates")
    ax.plot(np.arange(1, red.size + 1), np.log10(np.maximum(red, 1e-18)), "s-", color="#E69F00", label="Reduced, 6 rates")
    ax.axhline(math.log10(PRACTICAL_CUT), color="#333333", lw=0.8, ls="--", label="Practical cut")
    ax.axhline(math.log10(NUMERICAL_CUT), color="#999999", lw=0.8, ls=":", label="Numerical cut")
    ax.set_xlabel("Eigenvalue index")
    ax.set_ylabel("log10 (eigenvalue / leading)")
    ax.set_title("Shared four-channel Fisher spectrum")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_distances(scores: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=False)
    labels = [item["symbol"] for item in scores]
    x = np.arange(len(labels))
    reduced = [item["d_reduced"] / BUDGET_PROFILE for item in scores]
    full = [item["d_full"] / BUDGET_PROFILE for item in scores]
    colors = []
    for item in scores:
        colors.append({"survives": "#0072B2", "artefact_unreduced": "#D55E00", "confounded": "#666666"}[item["outcome"]])
    axes[0].bar(x, reduced, color=colors)
    axes[0].axhline(1.0, color="#333333", lw=0.8, ls="--")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Minimum chi-square / budget")
    axes[0].set_title("Reduced leakage, shared channels")
    axes[1].bar(x, full, color=colors)
    axes[1].axhline(1.0, color="#333333", lw=0.8, ls="--")
    axes[1].set_xticks(x, labels)
    axes[1].set_title("Full leakage, shared channels plus W")
    for ax in axes:
        ax.set_xlabel("Admitted forcing")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_signatures(bundle: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), sharex=True)
    panels = [
        (axes[0, 0], "u_p", "P", 2, "reduced chart"),
        (axes[0, 1], "u_q", "Q", 3, "reduced chart"),
        (axes[1, 0], "u_h", "P", 2, "reduced chart"),
        (axes[1, 1], "u_w", "W", 5, "full chart"),
    ]
    names = [item[0] for item in PROTOCOLS]
    x = np.arange(len(names))
    for ax, symbol, label, index, chart in panels:
        null_y = [bundle["null"][i][index] for i in range(4)]
        force_y = [bundle[symbol][i][index] for i in range(4)]
        leak_y = [bundle["leak_" + symbol][i][index] for i in range(4)]
        ax.plot(x, null_y, "o-", color="#666666", label="Null")
        ax.plot(x, force_y, "s-", color="#0072B2", label="Forcing")
        ax.plot(x, leak_y, "^-", color="#D55E00", label="Best leakage")
        ax.set_xticks(x, names)
        ax.set_ylabel(label)
        ax.set_title(f"{symbol}, {chart}")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Steady value against the best one-coordinate leakage", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_ledger(calls: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    x = np.arange(len(calls))
    admitted = np.array([1.0 if call["status"] == "admitted" else 0.0 for call in calls])
    held = np.array([1.0 if call["theta_unchanged"] else 0.0 for call in calls])
    ax.plot(x, held, "o-", color="#0072B2", label="Theta digest unchanged")
    ax.plot(x, admitted, "s", color="#E69F00", label="Forcing admitted")
    ax.set_ylim(-0.05, 1.15)
    ax.set_yticks([0, 1], ["no", "yes"])
    ax.set_xticks(x, [call["id"] for call in calls], rotation=90, fontsize=7)
    ax.set_xlabel("Call")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def steady_table(th: dict[str, float], forcing: dict[str, float]) -> list[dict]:
    rows = []
    for name, v_s, v_q in PROTOCOLS:
        state = steady_full(th, v_s, v_q, **forcing)
        if state is None:
            raise SystemExit("steady table left the domain")
        rows.append(
            {
                "protocol": name,
                "v_s": v_s,
                "v_q": v_q,
                "S": float(state[0]),
                "I": float(state[1]),
                "P": float(state[2]),
                "Q": float(state[3]),
                "Z": float(state[4]),
                "W": float(state[5]),
            }
        )
    return rows


def leakage_profile_states(th: dict[str, float], phi: dict[str, float], symbol: str, reduced_hit: dict) -> list[np.ndarray]:
    """States drawn for the figure: reduced leakage pushed to a full-looking tuple.

    W and Z are taken from the null full field, because the reduced chart has neither.
    P and Q use the reduced leakage. This is a display of the reduced counterfeit
    on the coordinates that chart owns.
    """
    moved = scale_table(phi, reduced_hit["coordinate"], reduced_hit["factor"])
    red = reduced_rows(moved)
    base = full_rows(th, null_forcing())
    if red is None or base is None:
        raise SystemExit("leakage profile undefined")
    display = []
    for red_row, base_row in zip(red, base):
        display.append(
            np.array([red_row[0], red_row[1], red_row[2], red_row[3], base_row[4], base_row[5]], dtype=float)
        )
    if symbol == "u_w":
        # The best reduced leakage of a readout forcing is the null shared map.
        # Show the full-field best leakage on W as well, so the panel is the full chart.
        return display
    return display


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    doc = load_candidates()
    if doc["library_sha256"] != hashlib.sha256(b"thesis-28 synthetic forcing catalogue v1").hexdigest():
        raise SystemExit("library pin is not the declared sentence digest")
    th = theta_float()
    phi = phi_from_theta(th)
    amplitude = float(PROTOCOL_AMPLITUDE)

    # Identities fixed by the algebra, checked before any ledger outcome is read.
    gaps = []
    for name, v_s, v_q in PROTOCOLS:
        full_null = steady_full(th, v_s, v_q)
        red_null = steady_reduced(phi, v_s, v_q)
        gaps.append(float(np.max(np.abs(full_null[:4] - red_null))))
        only_w = steady_full(th, v_s, v_q, **one_forcing("u_w", amplitude))
        gaps.append(float(np.max(np.abs(only_w[:4] - full_null[:4]))))
        if abs(only_w[5] - (full_null[5] + amplitude / th["b"])) > 1e-9:
            raise SystemExit("readout forcing did not shift W by u_w/b")
        forced_h = steady_full(th, v_s, v_q, **one_forcing("u_h", amplitude))
        phi_match = dict(phi)
        phi_match["kappa"] = kappa_of(th, amplitude)
        red_match = steady_reduced(phi_match, v_s, v_q)
        gaps.append(float(np.max(np.abs(forced_h[:4] - red_match))))
        rhs_gap = algebraic_gap(full_null, th, v_s, v_q, null_forcing())
        if rhs_gap > 1e-8:
            raise SystemExit(f"steady residual {rhs_gap}")
    if max(gaps) > 1e-9:
        raise SystemExit(f"algebraic identity failed: {max(gaps)}")

    factor_h = (th["h"] + amplitude) / th["h"]
    if not (BOX[0] <= factor_h <= BOX[1]):
        raise SystemExit("gain-shift counterfeit lies outside the declared box")

    def predict_shared_full(table: dict[str, float]) -> np.ndarray:
        rows = full_rows(table, null_forcing())
        if rows is None:
            raise SystemExit("shared prediction left the domain")
        return stack_shared(rows)

    def predict_shared_reduced(table: dict[str, float]) -> np.ndarray:
        rows = reduced_rows(table)
        if rows is None:
            raise SystemExit("reduced prediction left the domain")
        return stack_shared(rows)

    def predict_ponly_full(table: dict[str, float]) -> np.ndarray:
        return predict_shared_full(table)[2::4]

    def predict_ponly_reduced(table: dict[str, float]) -> np.ndarray:
        return predict_shared_reduced(table)[2::4]

    spec_full = full_fisher_with_plane(predict_shared_full, th)
    spec_red = fisher_rank(predict_shared_reduced, list(RED_NAMES), phi)
    spec_p_full = fisher_rank(predict_ponly_full, list(FULL_ORDER), th)
    spec_p_red = fisher_rank(predict_ponly_reduced, list(RED_NAMES), phi)

    def predict_readout_full(table: dict[str, float]) -> np.ndarray:
        rows = full_rows(table, null_forcing())
        if rows is None:
            raise SystemExit("readout prediction left the domain")
        return stack_readout(rows)

    spec_w = fisher_rank(predict_readout_full, list(FULL_ORDER), th)
    if spec_full["practical_rank"] != spec_red["practical_rank"]:
        raise SystemExit("shared practical rank was not preserved")
    if spec_full["practical_rank"] <= spec_p_full["practical_rank"]:
        raise SystemExit("four-channel schedule did not raise the practical rank")

    ledger = {
        "theta": {name: THETA[name] for name in THETA_NAMES},
        "forcing": {symbol: null_schedule(symbol) for symbol in FORCING_SYMBOLS},
        "calls": [],
    }
    baseline_sha = digest(ledger["theta"])
    for call in build_calls(doc):
        apply_call(doc, ledger, call)
    if digest(ledger["theta"]) != baseline_sha:
        raise SystemExit("Theta digest moved")
    admitted = [call for call in ledger["calls"] if call["status"] == "admitted"]
    if [call["admitted_symbol"] for call in admitted] != ["u_p", "u_q", "u_h", "u_w"]:
        raise SystemExit("admission set was not the four declared symbols")
    if any(not call["theta_unchanged"] for call in ledger["calls"]):
        raise SystemExit("a call changed Theta")

    scores = []
    profiles = {"null": full_rows(th, null_forcing())}
    for call in admitted:
        symbol = call["admitted_symbol"]
        forcing = one_forcing(symbol, amplitude)
        rows = full_rows(th, forcing)
        y_shared = stack_shared(rows)
        y_read = stack_readout(rows)
        red_hit = min_leakage_reduced(y_shared, phi)
        full_hit = min_leakage_full(y_read, th, "readout")
        d_reduced = float(red_hit["chi2"])
        d_full = float(full_hit["chi2"])
        outcome = classify(d_reduced, d_full)
        null_shared = chi2_distance(y_shared, stack_shared(full_rows(th, null_forcing())))
        null_read = chi2_distance(y_read, stack_readout(full_rows(th, null_forcing())))
        scores.append(
            {
                "call": call["id"],
                "row_id": call["admitted_row"],
                "symbol": symbol,
                "d_reduced": d_reduced,
                "d_full": d_full,
                "null_shared_chi2": null_shared,
                "null_readout_chi2": null_read,
                "reduced_coordinate": red_hit["coordinate"],
                "reduced_factor": red_hit["factor"],
                "full_coordinate": full_hit["coordinate"],
                "full_factor": full_hit["factor"],
                "reduced_runners": red_hit["runners"],
                "full_runners": full_hit["runners"],
                "outcome": outcome,
                "budget_profile": BUDGET_PROFILE,
                "budget_shared": BUDGET_SHARED,
                "budget_readout": BUDGET_READOUT,
                "reduced_inside_schedule_ball": d_reduced <= BUDGET_SHARED,
                "full_inside_schedule_ball": d_full <= BUDGET_READOUT,
            }
        )
        profiles[symbol] = rows
        if symbol == "u_w":
            moved = scale_table(th, full_hit["coordinate"], full_hit["factor"])
            profiles["leak_" + symbol] = full_rows(moved, null_forcing())
        else:
            profiles["leak_" + symbol] = leakage_profile_states(th, phi, symbol, red_hit)

    outcomes = {item["symbol"]: item["outcome"] for item in scores}
    if outcomes.get("u_h") != "confounded":
        raise SystemExit(f"gain shift was not confounded: {outcomes.get('u_h')}")
    if outcomes.get("u_w") != "artefact_unreduced":
        raise SystemExit(f"readout forcing was not an unreduced artefact: {outcomes.get('u_w')}")
    if outcomes.get("u_p") != "survives" or outcomes.get("u_q") != "survives":
        raise SystemExit(f"slow forcings did not survive: {outcomes}")

    witness = soft_prior_witness(doc, baseline_sha)
    gaps_int = {
        "null": integration_gap(th, null_forcing()),
        "u_p": integration_gap(th, one_forcing("u_p", amplitude)),
        "u_w": integration_gap(th, one_forcing("u_w", amplitude)),
        "u_h": integration_gap(th, one_forcing("u_h", amplitude)),
    }
    if max(gaps_int.values()) > 1e-6:
        raise SystemExit(f"trajectory left the steady map: {gaps_int}")

    elig = eligibility(doc)
    plot_spectra(spec_full, spec_red, FIG / "fim_spectra.png")
    plot_distances(scores, FIG / "chi_square_ratios.png")
    plot_signatures(profiles, FIG / "steady_signatures.png")
    plot_ledger(ledger["calls"], FIG / "ledger_digest.png")

    null_states = steady_table(th, null_forcing())
    results = {
        "seed": SEED,
        "rng_draws": 0,
        "sigma": SIGMA,
        "box": list(BOX),
        "practical_cut": PRACTICAL_CUT,
        "numerical_cut": NUMERICAL_CUT,
        "budget_profile": BUDGET_PROFILE,
        "budget_shared": BUDGET_SHARED,
        "budget_readout": BUDGET_READOUT,
        "n_shared": N_SHARED,
        "n_readout": N_READOUT,
        "amplitude": PROTOCOL_AMPLITUDE,
        "kappa": kappa_of(th, 0.0),
        "rho": rho_of(th),
        "factor_h_counterfeit": factor_h,
        "theta_sha256": baseline_sha,
        "theta": {name: THETA[name] for name in THETA_NAMES},
        "forcing_sha256_final": digest(ledger["forcing"]),
        "eligibility": elig,
        "calls": ledger["calls"],
        "n_refused": sum(call["status"] == "refused" for call in ledger["calls"]),
        "n_admitted": len(admitted),
        "ranks": {
            "shared_full": spec_full,
            "shared_reduced": spec_red,
            "product_only_full": spec_p_full,
            "product_only_reduced": spec_p_red,
            "readout_full": spec_w,
            "readout_reduced": None,
            "shared_practical_preserved": spec_full["practical_rank"] == spec_red["practical_rank"],
            "multi_channel_advantage_full": spec_full["practical_rank"] - spec_p_full["practical_rank"],
            "multi_channel_advantage_reduced": spec_red["practical_rank"] - spec_p_red["practical_rank"],
        },
        "scores": scores,
        "null_steady": null_states,
        "forced_steady": {
            symbol: steady_table(th, one_forcing(symbol, amplitude)) for symbol in ("u_p", "u_q", "u_h", "u_w")
        },
        "integration_gap": gaps_int,
        "algebraic_identity_gap": max(gaps),
        "soft_prior_witness": witness,
        "candidate_pin": CANDIDATE_PIN,
    }
    payload = jsonable(results)
    (ROOT / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"theta {baseline_sha}")
    print(
        "ranks shared",
        spec_full["practical_rank"],
        "/",
        spec_full["dimension"],
        "vs",
        spec_red["practical_rank"],
        "/",
        spec_red["dimension"],
        "product",
        spec_p_full["practical_rank"],
        "vs",
        spec_p_red["practical_rank"],
        "advantage",
        payload["ranks"]["multi_channel_advantage_full"],
        payload["ranks"]["multi_channel_advantage_reduced"],
    )
    print("kappa", payload["kappa"], "rho", payload["rho"], "budgets", BUDGET_SHARED, BUDGET_READOUT)
    for item in scores:
        print(
            item["symbol"],
            item["outcome"],
            "Dred",
            item["d_reduced"],
            item["reduced_coordinate"],
            item["reduced_factor"],
            "Dfull",
            item["d_full"],
            item["full_coordinate"],
            item["full_factor"],
            "null",
            item["null_shared_chi2"],
            item["null_readout_chi2"],
        )
    print("integration", gaps_int)
    print("witness", witness["posterior_h"], witness["fictional_theta_sha256"])
    print("refused", payload["n_refused"], "admitted", payload["n_admitted"])


if __name__ == "__main__":
    main()
