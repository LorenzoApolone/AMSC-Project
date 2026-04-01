#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re

COMM_RE = re.compile(r"^Comm times: Total ([0-9.]+) s \| Compute ([0-9.]+) s \| Bcast ([0-9.]+) s \| Allreduce ([0-9.]+) s \| Allgather ([0-9.]+) s \| Wait ([0-9.]+) s$")
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

RESULT_RE = re.compile(r"^([A-Za-z]+(?:[A-Za-z0-9_-]+)*)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.]+)s$")
TYPOLOGY_RE = re.compile(r"^\s*-\s*(.+?):\s*(\d+)/(\d+)\s*$")
CONVERGED_BY_TYPOLOGY_HEADER = "Converged by typology:"
NON_CONVERGED_HEADER = "Non-converged functions:"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse, aggregate and plot DPSO benchmark results.")
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--thresholds", default="0.1,0.01,0.001,0.0001,0.00001,0.000001")
    return parser.parse_args()


def parse_thresholds(raw: str) -> list[float]:
    return [float(token.strip()) for token in raw.split(",") if token.strip()]


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.fmean(values) if values else math.nan


def stdev(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.stdev(values) if len(values) > 1 else 0.0


def parse_float(raw: Any) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return math.nan


def fmt(value: float, digits: int = 9) -> str:
    return f"{value:.{digits}f}" if math.isfinite(value) else "nan"


def fmt_ratio(value: float) -> str:
    return f"{value:.6f}" if math.isfinite(value) else "nan"


def sanitize_name(raw: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("_")


def run_label(execution_mode: str, mpi_processes: int) -> str:
    return "serial" if execution_mode == "serial" else f"np={mpi_processes}"


def mode_sort_key(execution_mode: str, mpi_processes: int) -> tuple[int, int]:
    return (0 if execution_mode == "serial" else 1, mpi_processes)


def classify_function_outcome(status: str) -> str:
    if status == "True Convergence":
        return "success"
    if status.startswith("Failure"):
        return "failure"
    return "no_convergence"


def classify_run_outcome(return_code: int) -> str:
    return "completed" if return_code == 0 else "failure"


def config_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (row["battery"], row["suite"], row["family"], row["sub_swarm_size"])


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_stdout_log(log_path: Path) -> dict[str, Any]:
    parsed = {"functions": [], "typologies": [], "non_converged_lines": []}
    section = None
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = RESULT_RE.match(line)
        if match:
            current_function = match.group(1).strip()
            dim = int(match.group(2).strip())
            total_parts = int(match.group(3).strip())
            iters = int(match.group(4).strip())
            best_fitness = float(match.group(5).strip())
            fitness_error = float(match.group(6).strip())
            function_time_s = float(match.group(7).strip())
            status = "True Convergence" if fitness_error < 1e-4 else ("Failure" if iters >= 100 else "no_convergence")
            parsed["functions"].append({
                "function_name": current_function,
                "best_fitness": best_fitness,
                "fitness_error": fitness_error,
                "iterations": iters,
                "function_time_s": function_time_s,
                "comm_total_s": 0.0,
                "compute_total_s": 0.0,
                "comm_allgather_s": 0.0,
                "comm_bcast_s": 0.0,
                "comm_allreduce_s": 0.0,
                "comm_barrier_s": 0.0,
                "wait_total_s": 0.0,
                "status": status,
            })
            continue


        comm_match = COMM_RE.match(line)
        if comm_match and len(parsed["functions"]) > 0:
            last = parsed["functions"][-1]
            last["comm_total_s"] = float(comm_match.group(1))
            last["compute_total_s"] = float(comm_match.group(2))
            last["comm_bcast_s"] = float(comm_match.group(3))
            last["comm_allreduce_s"] = float(comm_match.group(4))
            last["comm_allgather_s"] = float(comm_match.group(5))
            last["wait_total_s"] = float(comm_match.group(6))
            continue
            
        stripped = line.strip()
        if "=== ALL BENCHMARKS COMPLETED ===" in stripped:
            section = "summary"
            continue
        if stripped == CONVERGED_BY_TYPOLOGY_HEADER:
            section = "typology"
            continue

        if section == "typology":
            if not stripped or stripped.startswith("- No typology metadata available."):
                continue
            match = TYPOLOGY_RE.match(line)
            if match:
                parsed["typologies"].append({"typology": match.group(1).strip(), "converged_count": int(match.group(2)), "total_count": int(match.group(3))})
            continue

        if section == "summary" and stripped.startswith("Function ") and stripped.endswith(" did not converge."):
            parsed["non_converged_lines"].append(stripped)

    return parsed

def aggregate_runs(run_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in run_rows:
        key = (row["battery"], row["suite"], row["family"], row["execution_mode"], row["mpi_processes"], row["dim"], row["sub_swarm_size"], row["total_particles"], row["max_iters"])
        grouped[key].append(row)
    serial_ref: dict[tuple[Any, ...], float] = {}
    parallel_baseline: dict[tuple[Any, ...], tuple[int, float]] = {}
    for key, rows in grouped.items():
        execution_mode, mpi_processes = key[3], key[4]
        cfg = (key[0], key[1], key[2], key[5], key[6], key[7], key[8])
        completed = [row for row in rows if row["run_outcome"] == "completed"]
        if not completed:
            continue
        current_mean = mean(row["suite_wall_time_s"] for row in completed)
        if execution_mode == "serial":
            serial_ref[cfg] = current_mean
        else:
            baseline = parallel_baseline.get(cfg)
            if baseline is None or mpi_processes < baseline[0]:
                parallel_baseline[cfg] = (mpi_processes, current_mean)
    summary_rows: list[dict[str, Any]] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1], item[2], item[5], item[3], item[4])):
        rows = grouped[key]
        execution_mode, mpi_processes = key[3], key[4]
        cfg = (key[0], key[1], key[2], key[5], key[6], key[7], key[8])
        completed = [row for row in rows if row["run_outcome"] == "completed"]
        failed_runs = len(rows) - len(completed)
        mean_wall = mean(row["suite_wall_time_s"] for row in completed)
        serial_time = serial_ref.get(cfg, math.nan)
        baseline = parallel_baseline.get(cfg)
        baseline_np = baseline[0] if baseline else None
        baseline_time = baseline[1] if baseline else math.nan
        if execution_mode == "parallel" and math.isfinite(mean_wall) and mean_wall > 0.0 and math.isfinite(baseline_time):
            speedup = baseline_time / mean_wall
            efficiency = speedup / mpi_processes if mpi_processes > 0 else math.nan
        else:
            speedup = math.nan
            efficiency = math.nan
        if execution_mode == "serial" and math.isfinite(mean_wall):
            rel_vs_serial = 1.0
        elif execution_mode == "parallel" and math.isfinite(serial_time) and mean_wall > 0.0:
            rel_vs_serial = serial_time / mean_wall
        else:
            rel_vs_serial = math.nan
        summary_rows.append({
            "battery": key[0], "suite": key[1], "family": key[2],
            "execution_mode": execution_mode, "run_label": run_label(execution_mode, mpi_processes), "mpi_processes": mpi_processes,
            "dim": key[5], "sub_swarm_size": key[6], "total_particles": key[7], "max_iters": key[8],
            "runs": len(rows), "completed_runs": len(completed), "failed_runs": failed_runs,
            "mean_suite_wall_time_s": fmt(mean_wall), "std_suite_wall_time_s": fmt(stdev(row["suite_wall_time_s"] for row in completed)),
            "mean_sum_function_time_s": fmt(mean(row["sum_function_time_s"] for row in completed)),
            "mean_sum_comm_total_s": fmt(mean(row["sum_comm_total_s"] for row in completed)),
            "mean_sum_compute_total_s": fmt(mean(row["sum_compute_total_s"] for row in completed)),
            "mean_sum_comm_allgather_s": fmt(mean(row["sum_comm_allgather_s"] for row in completed)),
            "mean_sum_comm_bcast_s": fmt(mean(row["sum_comm_bcast_s"] for row in completed)),
            "mean_sum_comm_allreduce_s": fmt(mean(row["sum_comm_allreduce_s"] for row in completed)),
            "mean_sum_comm_barrier_s": fmt(mean(row["sum_comm_barrier_s"] for row in completed)),
            "mean_sum_wait_total_s": fmt(mean(row["sum_wait_total_s"] for row in completed)),
            "mean_converged_functions": fmt_ratio(mean(row["converged_functions"] for row in completed)),
            "mean_convergence_rate": fmt_ratio(mean(row["convergence_rate"] for row in completed)),
            "mean_serial_reference_time_s": fmt(serial_time),
            "parallel_baseline_np": "" if baseline_np is None else baseline_np,
            "parallel_baseline_time_s": fmt(baseline_time),
            "speedup": fmt_ratio(speedup), "efficiency": fmt_ratio(efficiency),
            "relative_speedup_vs_serial": fmt_ratio(rel_vs_serial),
        })
    return summary_rows


def aggregate_typologies(typology_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in typology_rows:
        key = (row["battery"], row["suite"], row["family"], row["execution_mode"], row["mpi_processes"], row["dim"], row["typology"])
        grouped[key].append(row)
    out: list[dict[str, Any]] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1], item[2], item[5], item[3], item[4], item[6])):
        rows = grouped[key]
        total_mean = mean(row["total_count"] for row in rows)
        converged_mean = mean(row["converged_count"] for row in rows)
        ratio = converged_mean / total_mean if total_mean > 0 else math.nan
        out.append({
            "battery": key[0], "suite": key[1], "family": key[2],
            "execution_mode": key[3], "run_label": run_label(key[3], key[4]), "mpi_processes": key[4], "dim": key[5], "typology": key[6],
            "runs": len(rows), "mean_converged_count": fmt_ratio(converged_mean),
            "mean_total_count": fmt_ratio(total_mean), "mean_convergence_ratio": fmt_ratio(ratio),
        })
    return out


def build_precision_sweep(function_rows: list[dict[str, Any]], thresholds: list[float]) -> list[dict[str, Any]]:
    per_run: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in function_rows:
        if row["run_outcome"] != "completed":
            continue
        key = (row["battery"], row["suite"], row["family"], row["execution_mode"], row["mpi_processes"], row["dim"], row["sub_swarm_size"], row["total_particles"], row["max_iters"], row["seed"])
        per_run[key].append(row)
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for key, rows in per_run.items():
        total_functions = len(rows)
        for threshold in thresholds:
            converged = sum(1 for row in rows if row["fitness_error"] <= threshold)
            grouped[(key[0], key[1], key[2], key[3], key[4], key[5], key[6], key[7], key[8], threshold)].append({
                "converged_functions": converged,
                "convergence_ratio": converged / total_functions if total_functions else math.nan,
            })
    out: list[dict[str, Any]] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1], item[2], item[5], item[3], item[4], item[9])):
        rows = grouped[key]
        out.append({
            "battery": key[0], "suite": key[1], "family": key[2],
            "execution_mode": key[3], "run_label": run_label(key[3], key[4]), "mpi_processes": key[4], "dim": key[5],
            "sub_swarm_size": key[6], "total_particles": key[7], "max_iters": key[8], "threshold": f"{key[9]:.6g}",
            "runs": len(rows), "mean_converged_functions": fmt_ratio(mean(row["converged_functions"] for row in rows)),
            "std_converged_functions": fmt_ratio(stdev(row["converged_functions"] for row in rows)),
            "mean_convergence_ratio": fmt_ratio(mean(row["convergence_ratio"] for row in rows)),
        })
    return out


def aggregate_function_outcomes(function_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    per_run: dict[tuple[Any, ...], dict[str, float]] = defaultdict(lambda: {"success": 0.0, "no_convergence": 0.0, "failure": 0.0})
    for row in function_rows:
        key = (row["battery"], row["suite"], row["family"], row["execution_mode"], row["mpi_processes"], row["dim"], row["sub_swarm_size"], row["total_particles"], row["max_iters"], row["seed"])
        per_run[key][row["function_outcome"]] += 1.0
    grouped: dict[tuple[Any, ...], list[dict[str, float]]] = defaultdict(list)
    for key, counts in per_run.items():
        grouped[key[:-1]].append(counts)
    out: list[dict[str, Any]] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1], item[2], item[5], item[3], item[4])):
        rows = grouped[key]
        success_mean = mean(row["success"] for row in rows)
        no_conv_mean = mean(row["no_convergence"] for row in rows)
        failure_mean = mean(row["failure"] for row in rows)
        total_mean = success_mean + no_conv_mean + failure_mean
        out.append({
            "battery": key[0], "suite": key[1], "family": key[2], "execution_mode": key[3],
            "run_label": run_label(key[3], key[4]), "mpi_processes": key[4], "dim": key[5],
            "sub_swarm_size": key[6], "total_particles": key[7], "max_iters": key[8], "runs": len(rows),
            "mean_success_functions": fmt_ratio(success_mean), "mean_no_convergence_functions": fmt_ratio(no_conv_mean),
            "mean_failure_functions": fmt_ratio(failure_mean), "mean_total_functions": fmt_ratio(total_mean),
            "mean_success_ratio": fmt_ratio(success_mean / total_mean if total_mean > 0 else math.nan),
            "mean_no_convergence_ratio": fmt_ratio(no_conv_mean / total_mean if total_mean > 0 else math.nan),
            "mean_failure_ratio": fmt_ratio(failure_mean / total_mean if total_mean > 0 else math.nan),
        })
    return out


def build_serial_parallel_summary(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    serial_by_cfg = {config_key(row): row for row in summary_rows if row["execution_mode"] == "serial"}
    out: list[dict[str, Any]] = []
    parallel_rows = [row for row in summary_rows if row["execution_mode"] == "parallel"]
    for row in sorted(parallel_rows, key=lambda item: (item["battery"], item["suite"], item["family"], item["dim"], item["mpi_processes"])):
        serial_row = serial_by_cfg.get(config_key(row))
        out.append({
            "battery": row["battery"], "suite": row["suite"], "family": row["family"], "mpi_processes": row["mpi_processes"], "run_label": row["run_label"],
            "dim": row["dim"], "sub_swarm_size": row["sub_swarm_size"], "total_particles": row["total_particles"], "max_iters": row["max_iters"],
            "parallel_runs": row["runs"], "parallel_completed_runs": row["completed_runs"],
            "serial_runs": serial_row["runs"] if serial_row else 0, "serial_completed_runs": serial_row["completed_runs"] if serial_row else 0,
            "mean_parallel_wall_time_s": row["mean_suite_wall_time_s"], "mean_serial_wall_time_s": serial_row["mean_suite_wall_time_s"] if serial_row else "nan",
            "parallel_speedup": row["speedup"], "parallel_efficiency": row["efficiency"], "relative_speedup_vs_serial": row["relative_speedup_vs_serial"],
            "mean_parallel_convergence_rate": row["mean_convergence_rate"], "mean_serial_convergence_rate": serial_row["mean_convergence_rate"] if serial_row else "nan",
        })
    return out


def finite_values(values: Iterable[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def choose_time_unit(values: Iterable[float]) -> tuple[float, str]:
    valid = [abs(value) for value in finite_values(values) if not math.isclose(value, 0.0, abs_tol=1e-18)]
    if not valid:
        return 1.0, "s"
    vmax = max(valid)
    if vmax < 1.0e-3:
        return 1.0e6, "us"
    if vmax < 1.0:
        return 1.0e3, "ms"
    return 1.0, "s"


def scale_time_series(*series: list[float]) -> tuple[str, list[list[float]]]:
    all_values: list[float] = []
    for values in series:
        all_values.extend(finite_values(values))
    scale, unit = choose_time_unit(all_values)
    scaled = [
        [value * scale if math.isfinite(value) else value for value in values]
        for values in series
    ]
    return unit, scaled


def apply_focus_axis(axis: Any, values: Iterable[float], include_values: Iterable[float] = ()) -> None:
    valid = finite_values(list(values) + list(include_values))
    if not valid:
        return
    vmin = min(valid)
    vmax = max(valid)
    if math.isclose(vmin, vmax, rel_tol=1e-9, abs_tol=1e-12):
        pad = max(abs(vmax) * 0.15, 1e-12)
    else:
        pad = max((vmax - vmin) * 0.15, abs(vmax) * 0.03, 1e-12)
    lower = vmin - pad
    upper = vmax + pad
    if vmin > 0.0:
        lower = max(vmin * 0.85, lower)
    elif vmin >= 0.0:
        lower = max(0.0, lower)
    if not (math.isfinite(lower) and math.isfinite(upper)) or upper <= lower:
        lower, upper = vmin, vmax + max(abs(vmax) * 0.1, 1e-12)
    axis.set_ylim(lower, upper)


def apply_compact_time_axis(axis: Any, values: Iterable[float]) -> None:
    valid = finite_values(values)
    if not valid:
        return
    axis.ticklabel_format(style="plain", axis="y", useOffset=False)
    vmin = min(valid)
    vmax = max(valid)
    if math.isclose(vmin, vmax, rel_tol=1e-9, abs_tol=1e-12):
        pad = max(abs(vmax) * 0.2, 1e-12)
        lower = max(vmin * 0.85, vmin - pad) if vmin > 0.0 else max(0.0, vmin - pad)
        upper = vmax + pad
    else:
        span = vmax - vmin
        pad = max(span * 0.15, abs(vmax) * 0.05, 1e-12)
        lower = max(vmin * 0.85, vmin - pad) if vmin > 0.0 else max(0.0, vmin - pad)
        upper = vmax + pad
    if not (math.isfinite(lower) and math.isfinite(upper)) or upper <= lower:
        lower, upper = 0.0, max(1.0, vmax)
    axis.set_ylim(lower, upper)


def apply_positive_axis(axis: Any, values: Iterable[float]) -> None:
    valid = finite_values(values)
    if not valid:
        return
    apply_focus_axis(axis, valid)


def resolve_local_artifact_path(metadata_path: Path, raw_path: Any, fallback_name: str) -> Path:
    candidates: list[Path] = []
    if raw_path:
        candidates.append(Path(str(raw_path)))
    candidates.append(metadata_path.parent / fallback_name)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def safe_ratio(numerator: float, denominator: float) -> float:
    if math.isfinite(numerator) and math.isfinite(denominator) and denominator > 0.0:
        return numerator / denominator
    return math.nan


def propagated_ratio_error(numerator: float, numerator_std: float, denominator: float, denominator_std: float) -> float:
    ratio = safe_ratio(numerator, denominator)
    if not math.isfinite(ratio):
        return math.nan
    rel_terms: list[float] = []
    if math.isfinite(numerator_std) and numerator_std > 0.0 and numerator > 0.0:
        rel_terms.append((numerator_std / numerator) ** 2)
    if math.isfinite(denominator_std) and denominator_std > 0.0 and denominator > 0.0:
        rel_terms.append((denominator_std / denominator) ** 2)
    if not rel_terms:
        return 0.0
    return abs(ratio) * math.sqrt(sum(rel_terms))


def pretty_family_label(battery: str, suite: str, family: str) -> str:
    suite_labels = {
        ('comparable', 'strong_scaling_fixed_algorithm'): 'Reference strong scaling',
        ('dpso', 'weak_particles_per_process_constant'): 'Weak scaling with constant particles per process',
        ('dpso', 'weak_dimension_per_process_constant'): 'Weak scaling with constant dimension per process',
        ('dpso', 'weak_local_load_constant'): 'Weak scaling with constant local load',
        ('dpso', 'dimension_sweep_fixed_np8'): 'Dimension sweep with np=8',
        ('appendix', 'fixed_k28_owner_transition'): 'Fixed k=28 owner transition',
        ('validation', 'k28_consistency_check'): 'Fixed k=28 consistency check',
    }
    return suite_labels.get((battery, suite), family)


def maybe_plot(summary_rows: list[dict[str, Any]], serial_parallel_rows: list[dict[str, Any]], typology_summary_rows: list[dict[str, Any]], precision_rows: list[dict[str, Any]], function_outcome_rows: list[dict[str, Any]], output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('[warn] matplotlib not available; skipping plots.')
        return
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    for stale_png in plots_dir.glob('*.png'):
        stale_png.unlink()

    summary_lookup: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in summary_rows:
        key = (
            row['battery'], row['suite'], row['family'], row['execution_mode'], int(row['mpi_processes']),
            row['sub_swarm_size'],
        )
        summary_lookup[key] = row

    def lookup_summary_row(group_key: tuple[Any, ...], execution_mode: str, mpi_processes: int) -> dict[str, Any] | None:
        return summary_lookup.get((
            group_key[0], group_key[1], group_key[2], execution_mode, mpi_processes,
            group_key[3],
        ))

    def lookup_std(group_key: tuple[Any, ...], execution_mode: str, mpi_processes: int) -> float:
        row = lookup_summary_row(group_key, execution_mode, mpi_processes)
        return parse_float(row['std_suite_wall_time_s']) if row else math.nan

    parallel_groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        if row['execution_mode'] == 'parallel':
            parallel_groups[config_key(row)].append(row)
    for group_key, rows in parallel_groups.items():
        rows = sorted(rows, key=lambda item: int(item['mpi_processes']))
        rows = [row for row in rows if int(row.get('completed_runs', 0)) > 0 and math.isfinite(parse_float(row['mean_suite_wall_time_s']))]
        if not rows:
            continue
        np_values = [int(row['mpi_processes']) for row in rows]
        is_dim_sweep = (group_key[1] == "dimension_sweep_fixed_np8")
        x_label = 'Dimension' if is_dim_sweep else 'MPI processes'
        x_values = [int(row['dim']) for row in rows] if is_dim_sweep else np_values
        
        wall_times = [parse_float(row['mean_suite_wall_time_s']) for row in rows]
        wall_stds = [lookup_std(group_key, 'parallel', np_value) for np_value in np_values]
        speedups = [parse_float(row['speedup']) for row in rows]
        efficiencies = [parse_float(row['efficiency']) for row in rows]
        dim_disp = rows[0]["dim"] if group_key[1] != "weak_dimension_per_process_constant" and group_key[1] != "dimension_sweep_fixed_np8" else f"~{int(rows[0]['dim']) // max(1, int(rows[0]['mpi_processes']))}/rank"
        pps_disp = int(rows[0]["total_particles"]) // max(1, int(rows[0]["mpi_processes"]))
        slug = sanitize_name(f"{group_key[0]}__{group_key[1]}__{group_key[2]}__k{group_key[3]}")
        title_suffix = f"{pretty_family_label(group_key[0], group_key[1], group_key[2])} | dim={dim_disp} | k={group_key[3]} | pps={pps_disp}"

        wall_unit, scaled_wall = scale_time_series(wall_times, wall_stds)
        wall_times_plot, wall_stds_plot = scaled_wall
        figure, axis = plt.subplots(figsize=(7, 4))
        axis.errorbar(x_values, wall_times_plot, yerr=wall_stds_plot, marker='o', capsize=4)
        axis.set_title(f"Wall time vs {x_label}\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel(f'Mean suite wall time [{wall_unit}]')
        apply_compact_time_axis(axis, wall_times_plot)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__wall_time.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, speedups, marker='o', label='measured')
        baseline_np = max(1, int(rows[0]['mpi_processes']))
        ideal_speedups = [x_value / baseline_np for x_value in x_values] if not is_dim_sweep else [1.0 for _ in x_values]
        axis.plot(x_values, ideal_speedups, linestyle='--', color='tab:gray', label='ideal')
        axis.set_title(f"Speedup vs {x_label} (baseline: parallel np={baseline_np})\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel('parallel baseline time / parallel time')
        apply_focus_axis(axis, speedups, ideal_speedups)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__speedup_vs_parallel_baseline.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, efficiencies, marker='o', label='measured')
        ideal_efficiency = [1.0 for _ in x_values]
        axis.plot(x_values, ideal_efficiency, linestyle='--', color='tab:gray', label='ideal')
        axis.set_title(f"Efficiency vs {x_label} (baseline: parallel np={baseline_np})\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel('Speedup / MPI processes')
        apply_focus_axis(axis, efficiencies, ideal_efficiency)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__efficiency_vs_parallel_baseline.png", dpi=150)
        plt.close(figure)

        comm_totals = [parse_float(row['mean_sum_comm_total_s']) * float(row['mpi_processes']) for row in rows]
        wait_totals = [parse_float(row['mean_sum_wait_total_s']) * float(row['mpi_processes']) for row in rows]
        allgather_totals = [parse_float(row['mean_sum_comm_allgather_s']) * float(row['mpi_processes']) for row in rows]
        bcast_totals = [parse_float(row['mean_sum_comm_bcast_s']) * float(row['mpi_processes']) for row in rows]
        allreduce_totals = [parse_float(row['mean_sum_comm_allreduce_s']) * float(row['mpi_processes']) for row in rows]
        barrier_totals = [parse_float(row['mean_sum_comm_barrier_s']) * float(row['mpi_processes']) for row in rows]

        comm_unit, [comm_totals_plot] = scale_time_series(comm_totals)
        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, comm_totals_plot, marker='o')
        axis.set_title(f"Communication time vs {x_label}\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel(f'Aggregated communication time [{comm_unit}]')
        apply_compact_time_axis(axis, comm_totals_plot)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__comm_total.png", dpi=150)
        plt.close(figure)

        wait_unit, [wait_totals_plot] = scale_time_series(wait_totals)
        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, wait_totals_plot, marker='o')
        axis.set_title(f"Wait/sync time vs {x_label}\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel(f'Aggregated wait time [{wait_unit}]')
        apply_compact_time_axis(axis, wait_totals_plot)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__wait_total.png", dpi=150)
        plt.close(figure)

        breakdown_unit, breakdown_scaled = scale_time_series(allgather_totals, bcast_totals, allreduce_totals, barrier_totals)
        allgather_plot, bcast_plot, allreduce_plot, barrier_plot = breakdown_scaled
        figure, axis = plt.subplots(figsize=(8, 4.5))
        axis.plot(x_values, allgather_plot, marker='o', label='allgather')
        axis.plot(x_values, bcast_plot, marker='o', label='bcast')
        axis.plot(x_values, allreduce_plot, marker='o', label='allreduce')
        axis.plot(x_values, barrier_plot, marker='o', label='barrier')
        axis.set_title(f"Communication breakdown vs {x_label}\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel(f'Aggregated time [{breakdown_unit}]')
        apply_compact_time_axis(axis, allgather_plot + bcast_plot + allreduce_plot + barrier_plot)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__comm_breakdown.png", dpi=150)
        plt.close(figure)

    serial_groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in serial_parallel_rows:
        serial_groups[(row['battery'], row['suite'], row['family'], row['sub_swarm_size'])].append(row)
    for group_key, rows in serial_groups.items():
        rows = sorted(rows, key=lambda item: int(item['mpi_processes']))
        serial_time = parse_float(rows[0]['mean_serial_wall_time_s'])
        if not math.isfinite(serial_time):
            continue
        serial_std = lookup_std(group_key, 'serial', 0)
        serial_conv = parse_float(rows[0]['mean_serial_convergence_rate'])
        x_labels = ['serial'] + [row['run_label'] for row in rows]
        positions = list(range(len(x_labels)))
        np_values = [int(row['mpi_processes']) for row in rows]
        is_dim_sweep = (group_key[1] == "dimension_sweep_fixed_np8")
        x_label = 'Dimension' if is_dim_sweep else 'MPI processes'
        x_values = [int(row['dim']) for row in rows] if is_dim_sweep else np_values

        parallel_wall_times = [parse_float(row['mean_parallel_wall_time_s']) for row in rows]
        parallel_stds = [lookup_std(group_key, 'parallel', np_value) for np_value in np_values]
        wall_times = [serial_time] + parallel_wall_times
        wall_stds = [serial_std] + parallel_stds
        conv_rates = [serial_conv] + [parse_float(row['mean_parallel_convergence_rate']) for row in rows]
        speedups_vs_serial = [safe_ratio(serial_time, wall_time) for wall_time in parallel_wall_times]
        efficiencies_vs_serial = [safe_ratio(speedup, np_value) for speedup, np_value in zip(speedups_vs_serial, np_values)]
        wall_ratios_vs_serial = [safe_ratio(wall_time, serial_time) for wall_time in parallel_wall_times]
        wall_ratio_stds = [
            propagated_ratio_error(wall_time, wall_std, serial_time, serial_std)
            for wall_time, wall_std in zip(parallel_wall_times, parallel_stds)
        ]
        dim_disp = rows[0]["dim"] if group_key[1] != "weak_dimension_per_process_constant" and group_key[1] != "dimension_sweep_fixed_np8" else f"~{int(rows[0]['dim']) // max(1, int(rows[0]['mpi_processes']))}/rank"
        pps_disp = int(rows[0]["total_particles"]) // max(1, int(rows[0]["mpi_processes"]))
        slug = sanitize_name(f"{group_key[0]}__{group_key[1]}__{group_key[2]}__k{group_key[3]}")
        title_suffix = f"{pretty_family_label(group_key[0], group_key[1], group_key[2])} | dim={dim_disp} | k={group_key[3]} | pps={pps_disp}"

        serial_parallel_wall_unit, scaled_wall = scale_time_series(wall_times, wall_stds)
        wall_times_plot, wall_stds_plot = scaled_wall
        figure, axis = plt.subplots(figsize=(8, 4.5))
        axis.errorbar(positions, wall_times_plot, yerr=wall_stds_plot, marker='o', capsize=4)
        axis.set_xticks(positions, x_labels)
        axis.set_title(f"Serial and parallel wall time\n{title_suffix}")
        axis.set_xlabel('Execution mode')
        axis.set_ylabel(f'Mean suite wall time [{serial_parallel_wall_unit}]')
        apply_compact_time_axis(axis, wall_times_plot)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__serial_parallel_wall_time.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, speedups_vs_serial, marker='o', label='measured')
        ideal_speedups = [float(x_value) for x_value in x_values] if not is_dim_sweep else [1.0 for _ in x_values]
        axis.plot(x_values, ideal_speedups, linestyle='--', color='tab:gray', label='ideal')
        axis.set_title(f"Speedup vs serial reference\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel('serial time / parallel time')
        apply_focus_axis(axis, speedups_vs_serial, ideal_speedups)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__speedup_vs_serial.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(7, 4))
        axis.plot(x_values, efficiencies_vs_serial, marker='o', label='measured')
        ideal_efficiencies = [1.0 for _ in x_values]
        axis.plot(x_values, ideal_efficiencies, linestyle='--', color='tab:gray', label='ideal')
        axis.set_title(f"Efficiency vs serial reference\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel('(serial speedup) / MPI processes')
        apply_focus_axis(axis, efficiencies_vs_serial, ideal_efficiencies)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__efficiency_vs_serial.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(7, 4))
        axis.errorbar(x_values, wall_ratios_vs_serial, yerr=wall_ratio_stds, marker='o', capsize=4, label='measured')
        ideal_ratios = [1.0 / x_value for x_value in x_values] if not is_dim_sweep else [1.0] * len(x_values)
        axis.plot(x_values, ideal_ratios, linestyle='--', color='tab:gray', label='ideal')
        axis.axhline(1.0, linestyle=':', color='tab:red', linewidth=1.2, label='serial parity')
        axis.set_title(f"Parallel wall time normalized by serial\n{title_suffix}")
        axis.set_xlabel(x_label)
        axis.set_ylabel('parallel time / serial time')
        apply_focus_axis(axis, wall_ratios_vs_serial, ideal_ratios + [1.0])
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__parallel_to_serial_wall_ratio.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(8, 4.5))
        axis.plot(positions, conv_rates, marker='o')
        axis.set_xticks(positions, x_labels)
        axis.set_title(f"Serial vs parallel convergence ratio\n{title_suffix}")
        axis.set_xlabel('Execution mode')
        axis.set_ylabel('Mean convergence ratio')
        axis.set_ylim(0.0, 1.05)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        figure.savefig(plots_dir / f"{slug}__serial_vs_parallel_convergence.png", dpi=150)
        plt.close(figure)

    family_dim_summary: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        family_dim_summary[(row['battery'], row['suite'], row['family'])].append(row)
    for family_key, rows in family_dim_summary.items():
        dims = sorted({int(row['dim']) for row in rows})
        labels = sorted({(row['execution_mode'], int(row['mpi_processes'])) for row in rows}, key=lambda item: mode_sort_key(item[0], item[1]))
        if len(dims) < 2 or len(labels) < 2:
            continue
        slug = sanitize_name('__'.join(family_key))
        wall_by_dim_scale, wall_by_dim_unit = choose_time_unit(parse_float(row['mean_suite_wall_time_s']) for row in rows)
        pretty_label = pretty_family_label(*family_key)

        figure, axis = plt.subplots(figsize=(8, 4.5))
        has_data = False
        plotted_values: list[float] = []
        for execution_mode, mpi_processes in labels:
            label_rows = sorted([row for row in rows if row['execution_mode'] == execution_mode and int(row['mpi_processes']) == mpi_processes], key=lambda item: int(item['dim']))
            if not label_rows:
                continue
            y_values = [parse_float(row['mean_suite_wall_time_s']) * wall_by_dim_scale for row in label_rows]
            has_data = True
            plotted_values.extend(y_values)
            axis.plot([int(row['dim']) for row in label_rows], y_values, marker='o', label=run_label(execution_mode, mpi_processes))
        if has_data:
            axis.set_title(f"Wall time vs dimension\n{pretty_label}")
            axis.set_xlabel('Problem dimension')
            axis.set_ylabel(f'Mean suite wall time [{wall_by_dim_unit}]')
            apply_compact_time_axis(axis, plotted_values)
            axis.grid(True, alpha=0.3)
            axis.legend(loc='best', fontsize='small')
            figure.tight_layout()
            figure.savefig(plots_dir / f"{slug}__wall_time_by_dim.png", dpi=150)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(8, 4.5))
        has_data = False
        for execution_mode, mpi_processes in labels:
            label_rows = sorted([row for row in rows if row['execution_mode'] == execution_mode and int(row['mpi_processes']) == mpi_processes], key=lambda item: int(item['dim']))
            if not label_rows:
                continue
            has_data = True
            axis.plot([int(row['dim']) for row in label_rows], [parse_float(row['mean_convergence_rate']) for row in label_rows], marker='o', label=run_label(execution_mode, mpi_processes))
        if has_data:
            axis.set_title(f"Convergence ratio vs dimension\n{pretty_label}")
            axis.set_xlabel('Problem dimension')
            axis.set_ylabel('Mean convergence ratio')
            axis.set_ylim(0.0, 1.05)
            axis.grid(True, alpha=0.3)
            axis.legend(loc='best', fontsize='small')
            figure.tight_layout()
            figure.savefig(plots_dir / f"{slug}__convergence_by_dim.png", dpi=150)
        plt.close(figure)

    typology_groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in typology_summary_rows:
        typology_groups[(row['battery'], row['suite'], row['family'])].append(row)
    for family_key, rows in typology_groups.items():
        labels = sorted({(row['execution_mode'], int(row['mpi_processes'])) for row in rows}, key=lambda item: mode_sort_key(item[0], item[1]))
        by_typology: dict[str, dict[tuple[str, int], float]] = defaultdict(dict)
        for row in rows:
            by_typology[row['typology']][(row['execution_mode'], int(row['mpi_processes']))] = parse_float(row['mean_convergence_ratio'])
        if not by_typology:
            continue
        x_labels = [run_label(mode, np_value) for mode, np_value in labels]
        positions = list(range(len(x_labels)))
        figure, axis = plt.subplots(figsize=(8.5, 4.8))
        has_data = False
        for typology, values in sorted(by_typology.items()):
            ratios = [values.get(label, math.nan) for label in labels]
            if not any(math.isfinite(value) for value in ratios):
                continue
            has_data = True
            axis.plot(positions, ratios, marker='o', label=typology)
        if not has_data:
            plt.close(figure)
            continue
        dim_disp = rows[0]["dim"] if family_key[1] != "weak_dimension_per_process_constant" and family_key[1] != "dimension_sweep_fixed_np8" else f"~{int(rows[0]['dim']) // max(1, int(rows[0]['mpi_processes']))}/rank"
        axis.set_xticks(positions, x_labels)
        axis.set_title(f"Convergence ratio by typology\n{pretty_family_label(family_key[0], family_key[1], family_key[2])} | dim={dim_disp}")
        axis.set_xlabel('Execution mode')
        axis.set_ylabel('Mean convergence ratio')
        axis.set_ylim(0.0, 1.05)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        slug = sanitize_name(f"{family_key[0]}__{family_key[1]}__{family_key[2]}")
        figure.savefig(plots_dir / f"{slug}__typology_convergence_ratio.png", dpi=150)
        plt.close(figure)

    precision_groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in precision_rows:
        precision_groups[(row['battery'], row['suite'], row['family'])].append(row)
    for family_key, rows in precision_groups.items():
        labels = sorted({(row['execution_mode'], int(row['mpi_processes'])) for row in rows}, key=lambda item: mode_sort_key(item[0], item[1]))
        figure, axis = plt.subplots(figsize=(8.5, 4.8))
        has_data = False
        for execution_mode, mpi_processes in labels:
            label_rows = sorted([row for row in rows if row['execution_mode'] == execution_mode and int(row['mpi_processes']) == mpi_processes], key=lambda item: float(item['threshold']))
            if not label_rows:
                continue
            has_data = True
            axis.plot([float(row['threshold']) for row in label_rows], [parse_float(row['mean_convergence_ratio']) for row in label_rows], marker='o', label=run_label(execution_mode, mpi_processes))
        if not has_data:
            plt.close(figure)
            continue
        dim_disp = rows[0]["dim"] if family_key[1] != "weak_dimension_per_process_constant" and family_key[1] != "dimension_sweep_fixed_np8" else f"~{int(rows[0]['dim']) // max(1, int(rows[0]['mpi_processes']))}/rank"
        axis.set_xscale('log')
        axis.invert_xaxis()
        axis.set_title(f"Precision sweep convergence ratio\n{pretty_family_label(family_key[0], family_key[1], family_key[2])} | dim={dim_disp}")
        axis.set_xlabel('Error threshold')
        axis.set_ylabel('Mean convergence ratio')
        axis.set_ylim(0.0, 1.05)
        axis.grid(True, alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        slug = sanitize_name(f"{family_key[0]}__{family_key[1]}__{family_key[2]}")
        figure.savefig(plots_dir / f"{slug}__precision_sweep_ratio.png", dpi=150)
        plt.close(figure)

    outcome_groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in function_outcome_rows:
        outcome_groups[(row['battery'], row['suite'], row['family'])].append(row)
    for family_key, rows in outcome_groups.items():
        rows = sorted(rows, key=lambda item: mode_sort_key(item['execution_mode'], int(item['mpi_processes'])))
        x_labels = [row['run_label'] for row in rows]
        positions = list(range(len(x_labels)))
        success_counts = [parse_float(row['mean_success_functions']) for row in rows]
        no_conv_counts = [parse_float(row['mean_no_convergence_functions']) for row in rows]
        failure_counts = [parse_float(row['mean_failure_functions']) for row in rows]
        stacked_bottom = [success_counts[idx] + no_conv_counts[idx] for idx in range(len(x_labels))]
        figure, axis = plt.subplots(figsize=(8.5, 4.8))
        axis.bar(positions, success_counts, label='success')
        axis.bar(positions, no_conv_counts, bottom=success_counts, label='no_convergence')
        axis.bar(positions, failure_counts, bottom=stacked_bottom, label='failure')
        dim_disp = rows[0]["dim"] if family_key[1] != "weak_dimension_per_process_constant" and family_key[1] != "dimension_sweep_fixed_np8" else f"~{int(rows[0]['dim']) // max(1, int(rows[0]['mpi_processes']))}/rank"
        axis.set_xticks(positions, x_labels)
        axis.set_title(f"Function outcomes\n{pretty_family_label(family_key[0], family_key[1], family_key[2])} | dim={dim_disp}")
        axis.set_xlabel('Execution mode')
        axis.set_ylabel('Mean functions')
        axis.grid(True, axis='y', alpha=0.3)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        slug = sanitize_name(f"{family_key[0]}__{family_key[1]}__{family_key[2]}")
        figure.savefig(plots_dir / f"{slug}__function_outcomes.png", dpi=150)
        plt.close(figure)


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else results_root / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = parse_thresholds(args.thresholds)
    run_rows: list[dict[str, Any]] = []
    function_rows: list[dict[str, Any]] = []
    typology_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    for metadata_path in sorted(results_root.rglob("metadata.json")):
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        execution_mode = metadata.get("execution_mode", "parallel")
        mpi_processes = int(metadata.get("mpi_processes", 0))
        current_label = run_label(execution_mode, mpi_processes)
        stdout_log = resolve_local_artifact_path(metadata_path, metadata.get("stdout_log"), "stdout.log")
        stderr_log = resolve_local_artifact_path(metadata_path, metadata.get("stderr_log"), "stderr.log")
        parsed = {"functions": [], "typologies": [], "non_converged_lines": []}
        if stdout_log.exists():
            parsed = parse_stdout_log(stdout_log)
        converged = sum(1 for item in parsed["functions"] if item.get("status") == "True Convergence")
        total_functions = len(parsed["functions"])
        run_outcome = classify_run_outcome(metadata["return_code"])
        last_completed_function = parsed["functions"][-1]["function_name"] if parsed["functions"] else ""
        run_rows.append({
            "battery": metadata["battery"], "suite": metadata["suite"], "family": metadata["family"], "case_id": metadata["case_id"],
            "execution_mode": execution_mode, "run_label": current_label, "mpi_processes": mpi_processes,
            "dim": metadata["dim"], "sub_swarm_size": metadata["sub_swarm_size"], "total_particles": metadata["total_particles"], "max_iters": metadata["max_iters"], "seed": metadata["seed"],
            "suite_wall_time_s": metadata["suite_wall_time_s"], "sum_function_time_s": sum(item["function_time_s"] for item in parsed["functions"]),
            "sum_comm_total_s": sum(item.get("comm_total_s", 0.0) for item in parsed["functions"]), "sum_compute_total_s": sum(item.get("compute_total_s", 0.0) for item in parsed["functions"]),
            "sum_comm_allgather_s": sum(item.get("comm_allgather_s", 0.0) for item in parsed["functions"]), "sum_comm_bcast_s": sum(item.get("comm_bcast_s", 0.0) for item in parsed["functions"]),
            "sum_comm_allreduce_s": sum(item.get("comm_allreduce_s", 0.0) for item in parsed["functions"]), "sum_comm_barrier_s": sum(item.get("comm_barrier_s", 0.0) for item in parsed["functions"]),
            "sum_wait_total_s": sum(item.get("wait_total_s", 0.0) for item in parsed["functions"]),
            "functions_total": total_functions, "converged_functions": converged, "convergence_rate": converged / total_functions if total_functions else math.nan,
            "run_outcome": run_outcome, "last_completed_function": last_completed_function, "return_code": metadata["return_code"], "stdout_log": str(stdout_log), "stderr_log": str(stderr_log),
        })
        if run_outcome == "failure":
            failure_rows.append({
                "battery": metadata["battery"], "suite": metadata["suite"], "family": metadata["family"], "case_id": metadata["case_id"],
                "execution_mode": execution_mode, "run_label": current_label, "mpi_processes": mpi_processes,
                "dim": metadata["dim"], "sub_swarm_size": metadata["sub_swarm_size"], "total_particles": metadata["total_particles"], "max_iters": metadata["max_iters"], "seed": metadata["seed"],
                "suite_wall_time_s": metadata["suite_wall_time_s"], "completed_functions": total_functions, "last_completed_function": last_completed_function,
                "return_code": metadata["return_code"], "stdout_log": str(stdout_log), "stderr_log": str(stderr_log),
            })
        for item in parsed["functions"]:
            function_rows.append({
                "battery": metadata["battery"], "suite": metadata["suite"], "family": metadata["family"], "case_id": metadata["case_id"],
                "execution_mode": execution_mode, "run_label": current_label, "mpi_processes": mpi_processes,
                "dim": metadata["dim"], "sub_swarm_size": metadata["sub_swarm_size"], "total_particles": metadata["total_particles"], "max_iters": metadata["max_iters"], "seed": metadata["seed"],
                "run_outcome": run_outcome, "function_outcome": classify_function_outcome(item["status"]), **item,
            })
        for item in parsed["typologies"]:
            typology_rows.append({
                "battery": metadata["battery"], "suite": metadata["suite"], "family": metadata["family"],
                "execution_mode": execution_mode, "run_label": current_label, "mpi_processes": mpi_processes, "dim": metadata["dim"], "seed": metadata["seed"], **item,
            })
    summary_rows = aggregate_runs(run_rows)
    typology_summary_rows = aggregate_typologies(typology_rows)
    precision_rows = build_precision_sweep(function_rows, thresholds)
    function_outcome_rows = aggregate_function_outcomes(function_rows)
    serial_parallel_rows = build_serial_parallel_summary(summary_rows)
    function_fieldnames = ["battery", "suite", "family", "case_id", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "seed", "run_outcome", "function_outcome", "function_name", "best_fitness", "fitness_error", "iterations", "function_time_s", "comm_total_s", "compute_total_s", "comm_allgather_s", "comm_bcast_s", "comm_allreduce_s", "comm_barrier_s", "wait_total_s", "status"]
    write_csv(output_dir / "runs.csv", run_rows, ["battery", "suite", "family", "case_id", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "seed", "suite_wall_time_s", "sum_function_time_s", "sum_comm_total_s", "sum_compute_total_s", "sum_comm_allgather_s", "sum_comm_bcast_s", "sum_comm_allreduce_s", "sum_comm_barrier_s", "sum_wait_total_s", "functions_total", "converged_functions", "convergence_rate", "run_outcome", "last_completed_function", "return_code", "stdout_log", "stderr_log"])
    write_csv(output_dir / "functions.csv", function_rows, function_fieldnames)
    write_csv(output_dir / "functions_converged.csv", [row for row in function_rows if row["function_outcome"] == "success"], function_fieldnames)
    write_csv(output_dir / "functions_no_convergence.csv", [row for row in function_rows if row["function_outcome"] == "no_convergence"], function_fieldnames)
    write_csv(output_dir / "functions_failure.csv", [row for row in function_rows if row["function_outcome"] == "failure"], function_fieldnames)
    write_csv(output_dir / "summary.csv", summary_rows, ["battery", "suite", "family", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "runs", "completed_runs", "failed_runs", "mean_suite_wall_time_s", "std_suite_wall_time_s", "mean_sum_function_time_s", "mean_sum_comm_total_s", "mean_sum_compute_total_s", "mean_sum_comm_allgather_s", "mean_sum_comm_bcast_s", "mean_sum_comm_allreduce_s", "mean_sum_comm_barrier_s", "mean_sum_wait_total_s", "mean_converged_functions", "mean_convergence_rate", "mean_serial_reference_time_s", "parallel_baseline_np", "parallel_baseline_time_s", "speedup", "efficiency", "relative_speedup_vs_serial"])
    write_csv(output_dir / "serial_parallel_summary.csv", serial_parallel_rows, ["battery", "suite", "family", "mpi_processes", "run_label", "dim", "sub_swarm_size", "total_particles", "max_iters", "parallel_runs", "parallel_completed_runs", "serial_runs", "serial_completed_runs", "mean_parallel_wall_time_s", "mean_serial_wall_time_s", "parallel_speedup", "parallel_efficiency", "relative_speedup_vs_serial", "mean_parallel_convergence_rate", "mean_serial_convergence_rate"])
    write_csv(output_dir / "typology_summary.csv", typology_summary_rows, ["battery", "suite", "family", "execution_mode", "run_label", "mpi_processes", "dim", "typology", "runs", "mean_converged_count", "mean_total_count", "mean_convergence_ratio"])
    write_csv(output_dir / "precision_sweep.csv", precision_rows, ["battery", "suite", "family", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "threshold", "runs", "mean_converged_functions", "std_converged_functions", "mean_convergence_ratio"])
    write_csv(output_dir / "function_outcome_summary.csv", function_outcome_rows, ["battery", "suite", "family", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "runs", "mean_success_functions", "mean_no_convergence_functions", "mean_failure_functions", "mean_total_functions", "mean_success_ratio", "mean_no_convergence_ratio", "mean_failure_ratio"])
    write_csv(output_dir / "failures.csv", failure_rows, ["battery", "suite", "family", "case_id", "execution_mode", "run_label", "mpi_processes", "dim", "sub_swarm_size", "total_particles", "max_iters", "seed", "suite_wall_time_s", "completed_functions", "last_completed_function", "return_code", "stdout_log", "stderr_log"])
    maybe_plot(summary_rows, serial_parallel_rows, typology_summary_rows, precision_rows, function_outcome_rows, output_dir)
    print(f"[done] Wrote aggregated outputs to {output_dir}")


if __name__ == "__main__":
    main()
