#!/usr/bin/env python3
"""
debug_dwell_comparison.py -- Diagnostic script for comparing Python vs MATLAB
CMOST simulation dwell-time statistics.

Runs a 10,000-patient CMOST13 simulation (no screening, no surveillance)
and outputs per-year statistics to CSV for comparison with MATLAB.

Usage:
    python debug_dwell_comparison.py

    # With per-event trace logging (writes debug_cancer_trace.csv):
    set CMOST_DEBUG_TRACE=1
    python debug_dwell_comparison.py

Outputs:
    debug_dwell_yearly_stats.csv   -- per-year cancer counts by pathway
    debug_dwell_summary.txt        -- overall summary statistics
    debug_cancer_trace.csv         -- per-event trace (if CMOST_DEBUG_TRACE=1)
"""

import os
import sys
import copy
import csv
import numpy as np

# Ensure we can import sibling modules
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calculate_sub import calculate_sub


def load_cmost13_settings():
    """Load CMOST13 settings and override for diagnostic run."""
    import importlib.util
    settings_path = os.path.join(_this_dir, 'settings', 'CMOST13.py')
    spec = importlib.util.spec_from_file_location('_cmost13', settings_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    variables = copy.deepcopy(mod.settings)

    # Override for diagnostic: 10K patients, no screening, no surveillance
    variables['Number_patients'] = 10000
    variables['Screening']['Mode'] = 'off'
    variables['Polyp_Surveillance'] = 'off'
    variables['Cancer_Surveillance'] = 'off'

    return variables


def collect_yearly_stats(data):
    """Extract per-year cancer statistics from simulation data."""
    rows = []
    for y in range(100):  # years 1-100
        prog = int(data['ProgressedCancer'][y]) if y < len(data['ProgressedCancer']) else 0
        prog_r = int(data['ProgressedCancerR'][y]) if y < len(data['ProgressedCancerR']) else 0

        # Count direct cancers (DirectCancer2 = de-novo cancers per year)
        direct = int(data['DirectCancer2'][y]) if y < len(data['DirectCancer2']) else 0
        direct_r = int(data['DirectCancer2R'][y]) if y < len(data['DirectCancer2R']) else 0

        # Collect dwell times recorded this year
        dwell_prog_row = data['DwellTimeProgression'][y, :]
        nz_prog = dwell_prog_row[dwell_prog_row != 0]
        dwell_fast_row = data['DwellTimeFastCancer'][y, :]
        nz_fast = dwell_fast_row[dwell_fast_row != 0]

        rows.append({
            'year': y + 1,
            'progressed_cancer': prog,
            'progressed_cancer_R': prog_r,
            'direct_cancer': direct,
            'direct_cancer_R': direct_r,
            'n_dwell_prog': len(nz_prog),
            'mean_dwell_prog': float(np.mean(nz_prog)) if len(nz_prog) > 0 else 0,
            'median_dwell_prog': float(np.median(nz_prog)) if len(nz_prog) > 0 else 0,
            'n_dwell_fast': len(nz_fast),
            'mean_dwell_fast': float(np.mean(nz_fast)) if len(nz_fast) > 0 else 0,
            'median_dwell_fast': float(np.median(nz_fast)) if len(nz_fast) > 0 else 0,
        })
    return rows


def compute_summary(data):
    """Compute overall summary statistics matching Evaluation.py output."""
    # Gather all dwell times (same approach as Evaluation.py)
    all_dwell_prog = []
    all_dwell_fast = []

    for y in range(100):
        row_p = data['DwellTimeProgression'][y, :]
        nz_p = row_p[row_p != 0]
        all_dwell_prog.extend(nz_p.tolist())

        row_f = data['DwellTimeFastCancer'][y, :]
        nz_f = row_f[row_f != 0]
        all_dwell_fast.extend(nz_f.tolist())

    # Replicate the Evaluation.py doubling (AllDwellCa appends each value twice)
    doubled_prog = []
    for v in all_dwell_prog:
        doubled_prog.append(v)
        doubled_prog.append(v)
    doubled_fast = []
    for v in all_dwell_fast:
        doubled_fast.append(v)
        doubled_fast.append(v)

    combined = doubled_prog + doubled_fast

    summary = {}
    summary['n_patients'] = data['n']
    summary['total_progressed_cancer'] = int(np.sum(data['ProgressedCancer']))
    summary['total_direct_cancer'] = int(np.sum(data['DirectCancer2']))
    summary['n_dwell_prog_raw'] = len(all_dwell_prog)
    summary['n_dwell_fast_raw'] = len(all_dwell_fast)

    if all_dwell_prog:
        summary['median_dwell_prog_raw'] = float(np.median(all_dwell_prog))
        summary['mean_dwell_prog_raw'] = float(np.mean(all_dwell_prog))
        summary['q25_dwell_prog_raw'] = float(np.quantile(all_dwell_prog, 0.25))
        summary['q75_dwell_prog_raw'] = float(np.quantile(all_dwell_prog, 0.75))
    else:
        summary['median_dwell_prog_raw'] = 0
        summary['mean_dwell_prog_raw'] = 0
        summary['q25_dwell_prog_raw'] = 0
        summary['q75_dwell_prog_raw'] = 0

    if all_dwell_fast:
        summary['median_dwell_fast_raw'] = float(np.median(all_dwell_fast))
        summary['mean_dwell_fast_raw'] = float(np.mean(all_dwell_fast))
    else:
        summary['median_dwell_fast_raw'] = 0
        summary['mean_dwell_fast_raw'] = 0

    # These match the Evaluation.py "DwellTimeProgressedCa" (using doubled values)
    if doubled_prog:
        summary['DwellTimeProgressedCa_doubled'] = float(np.median(doubled_prog))
    else:
        summary['DwellTimeProgressedCa_doubled'] = 0

    if doubled_fast:
        summary['DwellTimeFastCa_doubled'] = float(np.median(doubled_fast))
    else:
        summary['DwellTimeFastCa_doubled'] = 0

    if combined:
        summary['DwellTimeAllCa_doubled'] = float(np.median(combined))
    else:
        summary['DwellTimeAllCa_doubled'] = 0

    return summary


def main():
    print("=" * 70)
    print("CMOST Dwell Time Diagnostic")
    print("=" * 70)
    print()

    # Set a fixed seed for reproducibility
    np.random.seed(42)

    variables = load_cmost13_settings()
    print(f"Settings loaded: DwellSpeed = {variables.get('DwellSpeed', 'NOT SET')}")
    print(f"  Number_patients = {variables['Number_patients']}")
    print(f"  Screening.Mode = {variables['Screening']['Mode']}")
    print(f"  Polyp_Surveillance = {variables['Polyp_Surveillance']}")
    print(f"  Cancer_Surveillance = {variables['Cancer_Surveillance']}")
    print()

    handles = {'Variables': variables}
    handles, bm = calculate_sub(handles)

    if 'data' not in handles:
        print("ERROR: simulation did not produce data.")
        return

    data = handles['data']

    # Yearly stats CSV
    yearly = collect_yearly_stats(data)
    csv_path = os.path.join(_this_dir, 'debug_dwell_yearly_stats.csv')
    fieldnames = list(yearly[0].keys())
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(yearly)
    print(f"Wrote per-year stats to: {csv_path}")

    # Summary
    summary = compute_summary(data)
    summary_path = os.path.join(_this_dir, 'debug_dwell_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CMOST Dwell Time Diagnostic Summary\n")
        f.write("=" * 50 + "\n\n")
        for k, v in summary.items():
            line = f"  {k}: {v}"
            f.write(line + "\n")
            print(line)

    print(f"\nWrote summary to: {summary_path}")

    # Key comparison values
    print("\n" + "=" * 70)
    print("KEY COMPARISON VALUES (compare with MATLAB CMOST13 output):")
    print("=" * 70)
    print(f"  DwellTimeProgressedCa (median, doubled): {summary['DwellTimeProgressedCa_doubled']:.2f}")
    print(f"  DwellTimeFastCa (median, doubled):       {summary['DwellTimeFastCa_doubled']:.2f}")
    print(f"  DwellTimeAllCa (median, doubled):        {summary['DwellTimeAllCa_doubled']:.2f}")
    print(f"  Total progressed cancers:                {summary['total_progressed_cancer']}")
    print(f"  Total direct cancers:                    {summary['total_direct_cancer']}")
    print(f"  Raw progressed dwell events:             {summary['n_dwell_prog_raw']}")
    print(f"  Raw fast cancer dwell events:            {summary['n_dwell_fast_raw']}")
    print()
    print("  MATLAB CMOST13 reference: DwellTimeProgressedCa ~ 13 years")
    if summary['DwellTimeProgressedCa_doubled'] > 30:
        print("  ** STILL DIVERGENT ** -- further investigation needed")
    elif summary['DwellTimeProgressedCa_doubled'] > 0:
        print("  Result looks closer to MATLAB -- verify with full comparison")
    print()


if __name__ == '__main__':
    main()
