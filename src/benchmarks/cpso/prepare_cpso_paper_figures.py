#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent / 'results'
ANALYSIS = ROOT / 'analysis'
PAPER_ROOT = ROOT / 'paper_figures'
FIGURES_DIR = PAPER_ROOT / 'figures'
README_PATH = PAPER_ROOT / 'PLOTS_FOR_PAPER.md'

FIGURES = [
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim64_k32_pps8__dim64__k32__pps8__iters10000__wall_time.png',
        'target': '01_reference_strong_scaling_dim64_wall_time.png',
        'title': 'Reference strong scaling, dim=64, wall time',
        'purpose': 'Baseline wall-time reduction for the main CPSO campaign with the reference decomposition k=32, pps=8.',
        'use': 'Use this as the clean runtime baseline before discussing any k=28 sensitivity study.',
    },
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim64_k32_pps8__dim64__k32__pps8__iters10000__speedup_vs_serial.png',
        'target': '02_reference_strong_scaling_dim64_speedup_vs_serial.png',
        'title': 'Reference strong scaling, dim=64, speedup vs serial',
        'purpose': 'Serial-referenced speedup for the baseline campaign, suitable for paper presentation because it avoids the np=1 parallel baseline ambiguity.',
        'use': 'Use this instead of the parallel-baseline speedup when you want a direct serial-to-parallel comparison.',
    },
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim128_k32_pps8__dim128__k32__pps8__iters20000__wall_time.png',
        'target': '03_reference_strong_scaling_dim128_wall_time.png',
        'title': 'Reference strong scaling, dim=128, wall time',
        'purpose': 'Baseline wall-time reduction for the heavier reference case.',
        'use': 'Pair this with the dim=64 figure to show the same trend on a more expensive workload.',
    },
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim128_k32_pps8__dim128__k32__pps8__iters20000__speedup_vs_serial.png',
        'target': '04_reference_strong_scaling_dim128_speedup_vs_serial.png',
        'title': 'Reference strong scaling, dim=128, speedup vs serial',
        'purpose': 'Serial-referenced speedup for the heavier baseline case.',
        'use': 'Use this to keep the paper narrative consistent across dimensions.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__parallel_to_serial_wall_ratio.png',
        'target': '05_appendix_k28_dim64_parallel_to_serial_wall_ratio.png',
        'title': 'Appendix k=28, dim=64, normalized wall time',
        'purpose': 'Shows directly whether each parallel configuration is slower or faster than serial, making the np=1 penalty visible.',
        'use': 'This is the best figure to explain why superlinear-looking speedup appears when k=28 is compared against the np=1 parallel baseline.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__speedup_vs_serial.png',
        'target': '06_appendix_k28_dim64_speedup_vs_serial.png',
        'title': 'Appendix k=28, dim=64, speedup vs serial',
        'purpose': 'Serial-referenced speedup for the fixed-k=28 owner-transition study.',
        'use': 'Use this instead of the parallel-baseline speedup if you want a paper-safe performance figure.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__wait_total.png',
        'target': '07_appendix_k28_dim64_wait_time.png',
        'title': 'Appendix k=28, dim=64, wait/sync time',
        'purpose': 'Highlights the synchronization overhead trend across MPI process counts in the owner-transition study.',
        'use': 'Use this as the implementation-level explanation for the performance change when moving toward k=np.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__parallel_to_serial_wall_ratio.png',
        'target': '08_appendix_k28_dim128_parallel_to_serial_wall_ratio.png',
        'title': 'Appendix k=28, dim=128, normalized wall time',
        'purpose': 'Same normalized comparison as the dim=64 case, but on the heavier owner-transition workload.',
        'use': 'Use this to confirm that the k=28 effect persists on a more expensive problem.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__speedup_vs_serial.png',
        'target': '09_appendix_k28_dim128_speedup_vs_serial.png',
        'title': 'Appendix k=28, dim=128, speedup vs serial',
        'purpose': 'Serial-referenced speedup for the heavier appendix case.',
        'use': 'Use this if you want a compact performance figure for the k=28 appendix without exposing the superlinear parallel-baseline plot.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__wait_total.png',
        'target': '10_appendix_k28_dim128_wait_time.png',
        'title': 'Appendix k=28, dim=128, wait/sync time',
        'purpose': 'Synchronization trend for the heavier appendix case.',
        'use': 'Use this to connect the runtime behavior to communication/synchronization overhead.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__parallel_to_serial_wall_ratio.png',
        'target': '11_validation_k28_dim64_parallel_to_serial_wall_ratio.png',
        'title': 'Validation k=28, dim=64, normalized wall time',
        'purpose': 'Short validation battery focused on serial, np=1, np=4 and np=28, normalized by serial.',
        'use': 'Use this as the sanity-check plot to show that np=1 parallel is the real outlier, not serial.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__speedup_vs_serial.png',
        'target': '12_validation_k28_dim64_speedup_vs_serial.png',
        'title': 'Validation k=28, dim=64, speedup vs serial',
        'purpose': 'Serial-referenced validation speedup for the dim=64 sanity-check battery.',
        'use': 'Use this to confirm that the strong effect remains even when the validation battery is reduced and focused.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__serial_vs_parallel_convergence.png',
        'target': '13_validation_k28_dim64_convergence.png',
        'title': 'Validation k=28, dim=64, convergence ratio',
        'purpose': 'Checks that the validation speedup story is not hiding a clear convergence collapse.',
        'use': 'Use this as the numerical-quality companion to the dim=64 validation runtime plots.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__parallel_to_serial_wall_ratio.png',
        'target': '14_validation_k28_dim128_parallel_to_serial_wall_ratio.png',
        'title': 'Validation k=28, dim=128, normalized wall time',
        'purpose': 'Validation version of the serial-normalized wall-time plot on the heavier workload.',
        'use': 'Use this to show that the same interpretation survives on dim=128.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__speedup_vs_serial.png',
        'target': '15_validation_k28_dim128_speedup_vs_serial.png',
        'title': 'Validation k=28, dim=128, speedup vs serial',
        'purpose': 'Serial-referenced validation speedup for the heavier sanity-check battery.',
        'use': 'Use this to present the cleaned-up k=28 speedup story on dim=128.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__serial_vs_parallel_convergence.png',
        'target': '16_validation_k28_dim128_convergence.png',
        'title': 'Validation k=28, dim=128, convergence ratio',
        'purpose': 'Numerical-quality companion plot for the heavier validation case.',
        'use': 'Use this to reassure the reader that the runtime effect is not obviously bought by a convergence collapse.',
    },
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim64_k32_pps8__dim64__k32__pps8__iters10000__efficiency_vs_serial.png',
        'target': '17_reference_strong_scaling_dim64_efficiency_vs_serial.png',
        'title': 'Reference strong scaling, dim=64, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the main CPSO baseline on the lighter strong-scaling case.',
        'use': 'Use this only as a companion to the baseline speedup figure if you want an efficiency metric in the paper.',
    },
    {
        'source': ANALYSIS / 'cpso_all_scale_free_final' / 'plots' / 'comparable__strong_scaling_fixed_algorithm__cmp_strong_dim128_k32_pps8__dim128__k32__pps8__iters20000__efficiency_vs_serial.png',
        'target': '18_reference_strong_scaling_dim128_efficiency_vs_serial.png',
        'title': 'Reference strong scaling, dim=128, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the heavier baseline case.',
        'use': 'Use this as an optional supporting figure if you want to report efficiency together with the baseline speedup.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__efficiency_vs_serial.png',
        'target': '19_appendix_k28_dim64_efficiency_vs_serial.png',
        'title': 'Appendix k=28, dim=64, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the lighter owner-transition appendix case.',
        'use': 'Use this only as a supporting metric after the normalized wall-time and speedup-vs-serial figures.',
    },
    {
        'source': ANALYSIS / 'cpso_appendix_k28_final' / 'plots' / 'appendix__fixed_k28_owner_transition__appendix_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__efficiency_vs_serial.png',
        'target': '20_appendix_k28_dim128_efficiency_vs_serial.png',
        'title': 'Appendix k=28, dim=128, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the heavier owner-transition appendix case.',
        'use': 'Use this as optional support if you want the same efficiency metric on the heavier k=28 appendix workload.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim64_pps16__dim64__k28__pps16__iters10000__efficiency_vs_serial.png',
        'target': '21_validation_k28_dim64_efficiency_vs_serial.png',
        'title': 'Validation k=28, dim=64, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the lighter validation battery.',
        'use': 'Use this only if you want a compact efficiency confirmation for the dim=64 validation case.',
    },
    {
        'source': ANALYSIS / 'cpso_validation_k28_final' / 'plots' / 'validation__k28_consistency_check__validation_fixedk28_dim128_pps8__dim128__k28__pps8__iters20000__efficiency_vs_serial.png',
        'target': '22_validation_k28_dim128_efficiency_vs_serial.png',
        'title': 'Validation k=28, dim=128, efficiency vs serial',
        'purpose': 'Serial-referenced efficiency for the heavier validation battery.',
        'use': 'Use this as optional support if you want an efficiency check alongside the dim=128 validation speedup.',
    },
]


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for stale_png in FIGURES_DIR.glob('*.png'):
        stale_png.unlink()

    copied = []
    for figure in FIGURES:
        source = figure['source']
        if not source.exists():
            raise FileNotFoundError(f'Missing source figure: {source}')
        target = FIGURES_DIR / figure['target']
        shutil.copy2(source, target)
        copied.append((target.name, figure))

    lines = [
        '# CPSO Paper Figures',
        '',
        'This folder collects the CPSO figures that are safest and most informative to use in the paper.',
        '',
        'General note:',
        '- Prefer the `*_speedup_vs_serial` and `*_parallel_to_serial_wall_ratio` figures in the paper narrative.',
        '- Use `*_efficiency_vs_serial` only as a supporting metric, not as the main plot for the fixed-k=28 story.',
        '- Keep the `*_vs_parallel_baseline` figures only for internal diagnostics, because the fixed `k=28` studies can make the `np=1` parallel baseline artificially weak.',
        '- The `parallel_to_serial_wall_ratio` plots are the clearest way to show that `np=1` parallel is the outlier, not serial.',
        '',
        '## Included figures',
        '',
    ]
    for filename, figure in copied:
        lines.extend([
            f'### `{filename}`',
            '',
            f"- Title: {figure['title']}",
            f"- Source: `{figure['source']}`",
            f"- Purpose: {figure['purpose']}",
            f"- Recommended use: {figure['use']}",
            '',
        ])

    README_PATH.write_text('\n'.join(lines), encoding='utf-8', newline='\n')
    print(f'[done] Copied {len(copied)} paper figures into {FIGURES_DIR}')
    print(f'[done] Wrote guide to {README_PATH}')


if __name__ == '__main__':
    main()
