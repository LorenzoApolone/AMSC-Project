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
    k_subswarms: int
    particles_per_swarm: int
    max_iters: int
    shuffle_freq: int
    stagnation_patience: int
    seed: int
    runtime_profile: str
    env_overrides: str
    description: str


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
    k_subswarms: int,
    particles_per_swarm: int,
    max_iters: int,
    seed: int,
    description: str,
    runtime_profile: str = "default",
    env_overrides: str = "-",
    shuffle_freq: int = DEFAULT_SHUFFLE_FREQ,
    stagnation_patience: int = DEFAULT_STAGNATION_PATIENCE,
) -> BenchmarkCase:
    if execution_mode == "serial":
        case_id = (
            f"{family}__serial__dim{dim}__k{k_subswarms}"
            f"__pps{particles_per_swarm}__seed{seed}"
        )
    else:
        case_id = (
            f"{family}__np{mpi_processes:02d}__dim{dim}__k{k_subswarms}"
            f"__pps{particles_per_swarm}__seed{seed}"
        )
    return BenchmarkCase(
        battery=battery,
        suite=suite,
        family=family,
        case_id=case_id,
        execution_mode=execution_mode,
        mpi_processes=mpi_processes,
        dim=dim,
        k_subswarms=k_subswarms,
        particles_per_swarm=particles_per_swarm,
        max_iters=max_iters,
        shuffle_freq=shuffle_freq,
        stagnation_patience=stagnation_patience,
        seed=seed,
        runtime_profile=runtime_profile,
        env_overrides=env_overrides,
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
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 10000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione CPSO fissa, varia solo il numero di processi."
            ),
        },
        {
            "family": "cmp_strong_dim64_k32_pps8",
            "suite": "strong_scaling_fixed_algorithm",
            "dim": 64,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 10000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione CPSO fissa, varia solo il numero di processi."
            ),
        },
        {
            "family": "cmp_strong_dim128_k32_pps8",
            "suite": "strong_scaling_fixed_algorithm",
            "dim": 128,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 20000,
            "description": (
                "Benchmark comparabile inter-metodo: problema totale fisso, "
                "decomposizione CPSO fissa, varia solo il numero di processi."
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
                    k_subswarms=definition["k"],
                    particles_per_swarm=definition["particles_per_swarm"],
                    max_iters=definition["max_iters"],
                    seed=seed,
                    description=(
                        definition["description"]
                        + " Riferimento seriale per confronti implementativi seriale-vs-parallelo."
                    ),
                )
            )
            for mpi_processes in DEFAULT_PROCESSES:
                if mpi_processes > definition["k"]:
                    continue
                cases.append(
                    _make_parallel_case(
                        battery="comparable",
                        suite=definition["suite"],
                        family=definition["family"],
                        mpi_processes=mpi_processes,
                        dim=definition["dim"],
                        k_subswarms=definition["k"],
                        particles_per_swarm=definition["particles_per_swarm"],
                        max_iters=definition["max_iters"],
                        seed=seed,
                        description=definition["description"],
                    )
                )
    return cases


def build_appendix_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    appendix_definitions = (
        {
            "family": "appendix_fixedk28_dim64_pps16",
            "suite": "fixed_k28_owner_transition",
            "dim": 64,
            "k": 28,
            "particles_per_swarm": 16,
            "max_iters": 10000,
            "description": (
                "Appendice CPSO: problema totale fisso con k=28, "
                "varia solo il numero di processi per osservare la transizione "
                "verso il caso np=28 in cui ogni rank possiede un solo sottoswarm."
            ),
        },
        {
            "family": "appendix_fixedk28_dim128_pps8",
            "suite": "fixed_k28_owner_transition",
            "dim": 128,
            "k": 28,
            "particles_per_swarm": 8,
            "max_iters": 20000,
            "description": (
                "Appendice CPSO: problema totale fisso con k=28, "
                "varia solo il numero di processi per osservare la transizione "
                "verso il caso np=28 in cui ogni rank possiede un solo sottoswarm."
            ),
        },
    )

    cases: list[BenchmarkCase] = []
    for definition in appendix_definitions:
        for seed in seeds:
            common = dict(
                battery="appendix",
                suite=definition["suite"],
                family=definition["family"],
                dim=definition["dim"],
                k_subswarms=definition["k"],
                particles_per_swarm=definition["particles_per_swarm"],
                max_iters=definition["max_iters"],
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        definition["description"]
                        + " Riferimento seriale con la stessa decomposizione algoritmica."
                    ),
                )
            )
            for mpi_processes in DEFAULT_PROCESSES:
                if mpi_processes > definition["k"]:
                    continue
                cases.append(
                    _make_parallel_case(
                        **common,
                        mpi_processes=mpi_processes,
                        description=definition["description"],
                    )
                )
    return cases


def build_validation_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    validation_definitions = (
        {
            "family": "validation_fixedk28_dim64_pps16",
            "suite": "k28_consistency_check",
            "dim": 64,
            "k": 28,
            "particles_per_swarm": 16,
            "max_iters": 10000,
            "description": (
                "Micro-benchmark CPSO di validazione: controlla la stabilita' dei risultati "
                "nel caso k=28 confrontando seriale, np=1, np=4 e np=28 sul problema dim=64."
            ),
        },
        {
            "family": "validation_fixedk28_dim128_pps8",
            "suite": "k28_consistency_check",
            "dim": 128,
            "k": 28,
            "particles_per_swarm": 8,
            "max_iters": 20000,
            "description": (
                "Micro-benchmark CPSO di validazione: controlla la stabilita' dei risultati "
                "nel caso k=28 confrontando seriale, np=1, np=4 e np=28 sul problema dim=128."
            ),
        },
    )
    validation_processes = (1, 4, 28)

    cases: list[BenchmarkCase] = []
    for definition in validation_definitions:
        for seed in seeds:
            common = dict(
                battery="validation",
                suite=definition["suite"],
                family=definition["family"],
                dim=definition["dim"],
                k_subswarms=definition["k"],
                particles_per_swarm=definition["particles_per_swarm"],
                max_iters=definition["max_iters"],
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        definition["description"]
                        + " Riferimento seriale per il controllo del relativo speedup vs serial."
                    ),
                )
            )
            for mpi_processes in validation_processes:
                cases.append(
                    _make_parallel_case(
                        **common,
                        mpi_processes=mpi_processes,
                        description=definition["description"],
                    )
                )
    return cases


def build_communication_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    communication_processes = (8, 16, 28)
    communication_definitions = (
        {
            "family": "comm_baseline_dim64_k32_pps8",
            "suite": "greedy_merge_ablation",
            "dim": 64,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 10000,
            "runtime_profile": "baseline",
            "env_overrides": "-",
            "description": (
                "Ablation CPSO sulla cooperazione MPI: configurazione baseline "
                "con fallback greedy completo per misurare wall time, costo di comunicazione "
                "e convergenza ai processi medio-alti."
            ),
        },
        {
            "family": "comm_no_greedy_dim64_k32_pps8",
            "suite": "greedy_merge_ablation",
            "dim": 64,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 10000,
            "runtime_profile": "no_greedy_merge_fallback",
            "env_overrides": "CPSO_MPI_DISABLE_GREEDY_MERGE=1",
            "description": (
                "Ablation CPSO sulla cooperazione MPI: stessa configurazione baseline "
                "ma senza il fallback greedy nel loop principale, per stimare quanto "
                "pesa la cooperazione globale extra."
            ),
        },
        {
            "family": "comm_baseline_dim128_k32_pps8",
            "suite": "greedy_merge_ablation",
            "dim": 128,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 20000,
            "runtime_profile": "baseline",
            "env_overrides": "-",
            "description": (
                "Ablation CPSO sulla cooperazione MPI: configurazione baseline "
                "con fallback greedy completo per misurare wall time, costo di comunicazione "
                "e convergenza ai processi medio-alti."
            ),
        },
        {
            "family": "comm_no_greedy_dim128_k32_pps8",
            "suite": "greedy_merge_ablation",
            "dim": 128,
            "k": 32,
            "particles_per_swarm": 8,
            "max_iters": 20000,
            "runtime_profile": "no_greedy_merge_fallback",
            "env_overrides": "CPSO_MPI_DISABLE_GREEDY_MERGE=1",
            "description": (
                "Ablation CPSO sulla cooperazione MPI: stessa configurazione baseline "
                "ma senza il fallback greedy nel loop principale, per stimare quanto "
                "pesa la cooperazione globale extra."
            ),
        },
    )

    cases: list[BenchmarkCase] = []
    for definition in communication_definitions:
        for seed in seeds:
            for mpi_processes in communication_processes:
                if mpi_processes > definition["k"]:
                    continue
                cases.append(
                    _make_parallel_case(
                        battery="communication",
                        suite=definition["suite"],
                        family=definition["family"],
                        mpi_processes=mpi_processes,
                        dim=definition["dim"],
                        k_subswarms=definition["k"],
                        particles_per_swarm=definition["particles_per_swarm"],
                        max_iters=definition["max_iters"],
                        seed=seed,
                        runtime_profile=definition["runtime_profile"],
                        env_overrides=definition["env_overrides"],
                        description=definition["description"],
                    )
                )
    return cases


def build_cpso_only_cases(seeds: Sequence[int]) -> list[BenchmarkCase]:
    cases: list[BenchmarkCase] = []

    weak_particles_definitions = (
        {
            "family": "cpso_weak_particles_dim64_k_equals_p_pps32",
            "suite": "weak_particles_per_process_constant",
            "dim": 64,
            "particles_per_swarm": 32,
            "max_iters": 10000,
            "description": (
                "Benchmark CPSO weak scaling: particelle per processo costanti, "
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
                    battery="cpso",
                    suite=definition["suite"],
                    family=definition["family"],
                    dim=definition["dim"],
                    k_subswarms=mpi_processes,
                    particles_per_swarm=definition["particles_per_swarm"],
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
            max_iters = 20000
            common = dict(
                battery="cpso",
                suite="weak_dimension_per_process_constant",
                family="cpso_weak_dimensions_localdim8_totalparticles224",
                dim=dim,
                k_subswarms=mpi_processes,
                particles_per_swarm=particles_per_swarm,
                max_iters=max_iters,
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        "Benchmark CPSO weak scaling: dimensioni per processo costanti, "
                        "numero totale di particelle fisso, k uguale al numero di processi, "
                        "budget iterativo costante. Riferimento seriale con stessa decomposizione del caso parallelo."
                    ),
                )
            )
            cases.append(
                _make_parallel_case(
                    **common,
                    mpi_processes=mpi_processes,
                    description=(
                        "Benchmark CPSO weak scaling: dimensioni per processo costanti, "
                        "numero totale di particelle fisso, k uguale al numero di processi, "
                        "budget iterativo costante."
                    ),
                )
            )

    for seed in seeds:
        for mpi_processes in DEFAULT_PROCESSES:
            dim = 4 * mpi_processes
            max_iters = 20000
            cases.append(
                _make_parallel_case(
                    battery="cpso",
                    suite="weak_local_load_constant",
                    family="cpso_weak_local_particles32_dim_per_swarm4",
                    mpi_processes=mpi_processes,
                    dim=dim,
                    k_subswarms=mpi_processes,
                    particles_per_swarm=32,
                    max_iters=max_iters,
                    seed=seed,
                    description=(
                        "Benchmark CPSO weak scaling: carico locale costante, "
                        "particelle per sottoswarm costanti, dimensioni per sottoswarm costanti "
                        "e budget iterativo costante."
                    ),
                )
            )

    dimension_sweep = (16, 32, 64, 128, 256)
    dimension_sweep_max_iters = 20000
    for seed in seeds:
        for dim in dimension_sweep:
            common = dict(
                battery="cpso",
                suite="dimension_sweep_fixed_np8",
                family="cpso_dim_sweep_np8_k8_pps16",
                dim=dim,
                k_subswarms=8,
                particles_per_swarm=16,
                max_iters=dimension_sweep_max_iters,
                seed=seed,
            )
            cases.append(
                _make_serial_case(
                    **common,
                    description=(
                        "Benchmark specifico CPSO: sweep della dimensione con confronto seriale-vs-parallelo, "
                        "decomposizione fissa e budget iterativo costante."
                    ),
                )
            )
            cases.append(
                _make_parallel_case(
                    **common,
                    mpi_processes=8,
                    description=(
                        "Benchmark specifico CPSO: sweep della dimensione con numero di processi, "
                        "decomposizione e budget iterativo fissi per studiare convergenza e costo al crescere del problema."
                    ),
                )
            )

    return cases


def build_cases(battery: str, seeds: Sequence[int]) -> list[BenchmarkCase]:
    selected = battery.lower()
    cases: list[BenchmarkCase] = []

    if selected in {"all", "comparable"}:
        cases.extend(build_comparable_cases(seeds))
    if selected == "appendix":
        cases.extend(build_appendix_cases(seeds))
    if selected == "validation":
        cases.extend(build_validation_cases(seeds))
    if selected == "communication":
        cases.extend(build_communication_cases(seeds))
    if selected in {"all", "cpso"}:
        cases.extend(build_cpso_only_cases(seeds))
    if not cases:
        raise ValueError(f"Unsupported battery selection: {battery}")

    return cases


def emit_tsv(cases: Iterable[BenchmarkCase], include_header: bool) -> None:
    if include_header:
        print(
            "battery	suite	family	case_id	execution_mode	mpi_processes	dim	k_subswarms"
            "	particles_per_swarm	max_iters	shuffle_freq	stagnation_patience"
            "	seed	runtime_profile	env_overrides	description"
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
                    str(case.k_subswarms),
                    str(case.particles_per_swarm),
                    str(case.max_iters),
                    str(case.shuffle_freq),
                    str(case.stagnation_patience),
                    str(case.seed),
                    case.runtime_profile,
                    case.env_overrides,
                    case.description,
                ]
            )
        )


def emit_jsonl(cases: Iterable[BenchmarkCase]) -> None:
    for case in cases:
        print(json.dumps(asdict(case), ensure_ascii=True))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate benchmark matrices for comparable and CPSO-only campaigns."
    )
    parser.add_argument(
        "--battery",
        default="all",
        choices=("all", "comparable", "cpso", "appendix", "validation", "communication"),
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
