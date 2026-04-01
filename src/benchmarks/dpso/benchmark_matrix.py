#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence

DEFAULT_PROCESSES = (1, 2, 4, 8, 16, 28)
DEFAULT_SEEDS = (123, 456, 789, 2024, 4242)
DEFAULT_SHUFFLE_FREQ = 50
DEFAULT_STAGNATION_PATIENCE = 50


@dataclass(frozen=True)
class BenchmarkCase:
    battery: str
    suite: str
    family: str
    case_id: str
    execution_mode: str
    mpi_processes: int
    dim: int
    total_particles: int
    max_iters: int
    seed: int
    w: float = 0.729
    c1: float = 1.49445
    c2: float = 1.49445
    regrouping_period: int = 5
    sub_swarm_size: int = 5
    hmcr: float = 0.98
    par_min: float = 0.01
    par_max: float = 0.99
    description: str = ""


def parse_seed_list(raw: str | None) -> list[int]:
    if not raw:
        return list(DEFAULT_SEEDS)
    return [int(token.strip()) for token in raw.split(",") if token.strip()]


def _make_case(
    *,
    battery: str,
    suite: str,
    family: str,
    execution_mode: str,
    mpi_processes: int,
    dim: int,
    total_particles: int,
    max_iters: int,
    seed: int,
    w: float = 0.729,
    c1: float = 1.49445,
    c2: float = 1.49445,
    regrouping_period: int = 5,
    sub_swarm_size: int = 5,
    hmcr: float = 0.98,
    par_min: float = 0.01,
    par_max: float = 0.99,
    description: str = "",
) -> BenchmarkCase:
    if execution_mode == "serial":
        case_id = f"{family}__serial__dim{dim}__tp{total_particles}__seed{seed}"
    else:
        case_id = f"{family}__np{mpi_processes:02d}__dim{dim}__tp{total_particles}__seed{seed}"
    return BenchmarkCase(
        battery=battery,
        suite=suite,
        family=family,
        case_id=case_id,
        execution_mode=execution_mode,
        mpi_processes=mpi_processes,
        dim=dim,
        total_particles=total_particles,
        max_iters=max_iters,
        seed=seed,
        w=w,
        c1=c1,
        c2=c2,
        regrouping_period=regrouping_period,
        sub_swarm_size=sub_swarm_size,
        hmcr=hmcr,
        par_min=par_min,
        par_max=par_max,
        description=description,
    )


def _make_parallel_case(**kwargs) -> BenchmarkCase:
    return _make_case(execution_mode="parallel", **kwargs)


def _make_serial_case(**kwargs) -> BenchmarkCase:
    return _make_case(execution_mode="serial", mpi_processes=0, **kwargs)


def build_comparable_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    base_definitions = (
        {
            "family": "cmp_strong_dim32_k32_pps8",
            "suite": "strong_scaling_fixed_algorithm",
            "dim": 32,
            "total_particles": 256,
            "max_iters": 10000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione DPSO fissa, varia solo il numero di processi."
            ),
        },
        {
            "family": "cmp_strong_dim64_k32_pps8",
            "suite": "strong_scaling_fixed_algorithm",
            "dim": 64,
            "total_particles": 256,
            "max_iters": 10000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione DPSO fissa, varia solo il numero di processi."
            ),
        },
        {
            "family": "cmp_strong_dim128_k32_pps8",
            "suite": "strong_scaling_fixed_algorithm",
            "dim": 128,
            "total_particles": 256,
            "max_iters": 20000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione DPSO fissa, varia solo il numero di processi."
            ),
        },
    )

    cases: list[BenchmarkCase] = []
    for definition in base_definitions:
        for seed in seeds:
            cases.append(
                _make_serial_case(
                    battery="comparable",
                    suite=definition["suite"],
                    family=definition["family"],
                    dim=definition["dim"],
                    total_particles=definition["total_particles"],
                    max_iters=definition["max_iters"],
                    seed=seed,
                    description=(
                        definition["description"]
                        + " Riferimento seriale per confronti implementativi seriale-vs-parallelo."
                    ),
                )
            )
            for mpi_processes in DEFAULT_PROCESSES:
                if mpi_processes > 32:
                    continue
                cases.append(
                    _make_parallel_case(
                        battery="comparable",
                        suite=definition["suite"],
                        family=definition["family"],
                        mpi_processes=mpi_processes,
                        dim=definition["dim"],
                        total_particles=definition["total_particles"],
                        max_iters=definition["max_iters"],
                        seed=seed,
                        description=definition["description"],
                    )
                )
    return cases


def build_dpso_only_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    cases: list[BenchmarkCase] = []

    weak_particles_definitions = (
        {
            "family": "dpso_weak_particles_dim64_k_equals_p_pps32",
            "suite": "weak_particles_per_process_constant",
            "dim": 64,
            "particles_per_swarm": 32,
            "max_iters": 10000,
            "description": (
                "Benchmark DPSO weak scaling: particelle per processo costanti, "
                "dimensione globale fissa, k uguale al numero di processi."
            ),
        },
    )

    for definition in weak_particles_definitions:
        for seed in seeds:
            for mpi_processes in DEFAULT_PROCESSES:
                if mpi_processes > definition["dim"]:
                    continue
                common = dict(
                    battery="dpso",
                    suite=definition["suite"],
                    family=definition["family"],
                    dim=definition["dim"],
                    total_particles=mpi_processes * definition["particles_per_swarm"],
                    max_iters=definition["max_iters"],
                    seed=seed,
                )
                cases.append(
                    _make_serial_case(
                        **common,
                        description=(
                            definition["description"]
                            + " Riferimento seriale con stessa decomposizione del caso parallelo."
                        ),
                    )
                )
                cases.append(
                    _make_parallel_case(
                        **common,
                        mpi_processes=mpi_processes,
                        description=definition["description"],
                    )
                )

    weak_dimension_local_dim = 8
    weak_dimension_total_particles = 224
    for seed in seeds:
        for mpi_processes in DEFAULT_PROCESSES:
            dim = weak_dimension_local_dim * mpi_processes
            particles_per_swarm = weak_dimension_total_particles // mpi_processes
            max_iters = 10000 if dim <= 64 else 20000
            common = dict(
                battery="dpso",
                suite="weak_dimension_per_process_constant",
                family="dpso_weak_dimensions_localdim8_totalparticles224",
                dim=dim,
                total_particles=mpi_processes * particles_per_swarm,
                max_iters=max_iters,
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        "Benchmark DPSO weak scaling: dimensioni per processo costanti, "
                        "numero totale di particelle fisso, k uguale al numero di processi. "
                        "Riferimento seriale con stessa decomposizione del caso parallelo."
                    ),
                )
            )
            cases.append(
                _make_parallel_case(
                    **common,
                    mpi_processes=mpi_processes,
                    description=(
                        "Benchmark DPSO weak scaling: dimensioni per processo costanti, "
                        "numero totale di particelle fisso, k uguale al numero di processi."
                    ),
                )
            )

    for seed in seeds:
        for mpi_processes in DEFAULT_PROCESSES:
            dim = 4 * mpi_processes
            max_iters = 10000 if dim <= 64 else 20000
            cases.append(
                _make_parallel_case(
                    battery="dpso",
                    suite="weak_local_load_constant",
                    family="dpso_weak_local_particles32_dim_per_swarm4",
                    mpi_processes=mpi_processes,
                    dim=dim,
                    total_particles=mpi_processes * 32,
                    max_iters=max_iters,
                    seed=seed,
                    description=(
                        "Benchmark DPSO weak scaling: carico locale costante, "
                        "particelle per sottoswarm costanti e dimensioni per sottoswarm costanti."
                    ),
                )
            )

    dimension_sweep = (
        (16, 10000),
        (32, 10000),
        (64, 15000),
        (128, 20000),
        (256, 40000),
    )
    for seed in seeds:
        for dim, max_iters in dimension_sweep:
            common = dict(
                battery="dpso",
                suite="dimension_sweep_fixed_np8",
                family="dpso_dim_sweep_np8_k8_pps16",
                dim=dim,
                total_particles=128,
                max_iters=max_iters,
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        "Benchmark specifico DPSO: sweep della dimensione con confronto seriale-vs-parallelo, "
                        "mantenendo costante la stessa decomposizione algoritmica del caso parallelo."
                    ),
                )
            )
            cases.append(
                _make_parallel_case(
                    **common,
                    mpi_processes=8,
                    description=(
                        "Benchmark specifico DPSO: sweep della dimensione con numero di processi e "
                        "decomposizione fissi per studiare convergenza e costo al crescere del problema."
                    ),
                )
            )

    return cases


def build_cases(battery: str, seeds: Sequence[int]) -> list[BenchmarkCase]:
    selected = battery.lower()
    cases: list[BenchmarkCase] = []

    if selected in {"all", "comparable"}:
        cases.extend(build_comparable_cases(seeds))
    if selected in {"all", "dpso"}:
        cases.extend(build_dpso_only_cases(seeds))

    # Finer splitting to avoid wall time timeouts
    if selected == "dpso_weak_particles":
        all_dpso = build_dpso_only_cases(seeds)
        cases.extend([c for c in all_dpso if c.suite == "weak_particles_per_process_constant"])
    if selected == "dpso_weak_dimension":
        all_dpso = build_dpso_only_cases(seeds)
        cases.extend([c for c in all_dpso if c.suite == "weak_dimension_per_process_constant"])
    if selected == "dpso_weak_local":
        all_dpso = build_dpso_only_cases(seeds)
        cases.extend([c for c in all_dpso if c.suite == "weak_local_load_constant"])
    if selected == "dpso_sweep_dim":
        all_dpso = build_dpso_only_cases(seeds)
        cases.extend([c for c in all_dpso if c.suite == "dimension_sweep_fixed_np8"])

    if not cases:
        raise ValueError(f"Unsupported battery selection: {battery}")

    return cases


def emit_tsv(cases: Iterable[BenchmarkCase], include_header: bool) -> None:
    if include_header:
        print(
            "battery\tsuite\tfamily\tcase_id\texecution_mode\tmpi_processes\tdim\ttotal_particles"
            "\tmax_iters\tseed\tw\tc1\tc2\tregrouping_period\tsub_swarm_size\thmcr\tpar_min\tpar_max\tdescription"
        )
    for case in cases:
        print(
            "	".join(
                [
                    case.battery,
                    case.suite,
                    case.family,
                    case.case_id,
                    case.execution_mode,
                    str(case.mpi_processes),
                    str(case.dim),
                    str(case.total_particles),
                    str(case.max_iters),
                    str(case.seed),
                    str(case.w),
                    str(case.c1),
                    str(case.c2),
                    str(case.regrouping_period),
                    str(case.sub_swarm_size),
                    str(case.hmcr),
                    str(case.par_min),
                    str(case.par_max),
                    case.description,
                ]
            )
        )


def emit_jsonl(cases: Iterable[BenchmarkCase]) -> None:
    for case in cases:
        print(json.dumps(asdict(case), ensure_ascii=True))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate benchmark matrices for comparable and DPSO-only campaigns."
    )
    parser.add_argument(
        "--battery",
        default="all",
        choices=("all", "comparable", "dpso", "dpso_weak_particles", "dpso_weak_dimension", "dpso_weak_local", "dpso_sweep_dim"),
        help="Which battery to emit.",
    )
    parser.add_argument(
        "--format",
        default="tsv",
        choices=("tsv", "jsonl"),
        help="Output format.",
    )
    parser.add_argument(
        "--seeds",
        default=None,
        help="Comma-separated list of integer seeds.",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Suppress TSV header output.",
    )
    args = parser.parse_args()

    cases = build_cases(args.battery, parse_seed_list(args.seeds))
    if args.format == "jsonl":
        emit_jsonl(cases)
        return
    emit_tsv(cases, include_header=not args.no_header)


if __name__ == "__main__":
    main()
