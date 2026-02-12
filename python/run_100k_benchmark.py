#!/usr/bin/env python3
"""
run_100k_benchmark.py -- Run full 100K patient CMOST13 simulation and
print key metrics for comparison with MATLAB reference.

Usage:
    python run_100k_benchmark.py
"""

import os
import sys
import copy
import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calculate_sub import calculate_sub


def load_cmost13():
    import importlib.util
    settings_path = os.path.join(_this_dir, 'settings', 'CMOST13.py')
    spec = importlib.util.spec_from_file_location('_cmost13', settings_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return copy.deepcopy(mod.settings)


def main():
    np.random.seed(42)

    variables = load_cmost13()
    # CMOST13 defaults: 100K patients, screening off, surveillance on
    variables['Number_patients'] = 100000
    # Keep surveillance and screening as they are in CMOST13
    print("=" * 72)
    print("CMOST 100K Patient Simulation - CMOST13 Settings")
    print("=" * 72)
    print(f"  Number_patients       = {variables['Number_patients']}")
    print(f"  DwellSpeed            = {variables.get('DwellSpeed', 'NOT SET')}")
    print(f"  Screening.Mode        = {variables['Screening']['Mode']}")
    print(f"  Polyp_Surveillance    = {variables.get('Polyp_Surveillance')}")
    print(f"  Cancer_Surveillance   = {variables.get('Cancer_Surveillance')}")
    print()

    handles = {'Variables': variables}
    handles, bm = calculate_sub(handles)

    if 'data' not in handles:
        print("ERROR: simulation did not produce data.")
        return

    data = handles['data']
    n = data['n']

    # ================================================================
    # Collect summary metrics (matching Evaluation.py SummaryVariable)
    # ================================================================
    female_mask = data['Gender'] == 2
    male_mask = data['Gender'] == 1
    n_female = np.sum(female_mask)
    n_male = np.sum(male_mask)

    avg_age = np.sum(data['DeathYear']) / n - 1
    avg_age_m = np.sum(data['DeathYear'][male_mask]) / n_male - 1 if n_male > 0 else 0
    avg_age_f = np.sum(data['DeathYear'][female_mask]) / n_female - 1 if n_female > 0 else 0

    screening_colo = np.sum(data['Number']['Screening_Colonoscopy'])
    symptom_colo = np.sum(data['Number']['Symptoms_Colonoscopy'])
    followup_colo = np.sum(data['Number']['Follow_Up_Colonoscopy'])

    ca_deaths = int(np.sum(data['DeathCause'] == 2))
    ca_death_mask = data['DeathCause'] == 2
    yrs_lost_ca = np.sum(data['NaturalDeathYear'][ca_death_mask] - data['DeathYear'][ca_death_mask])

    colo_deaths = int(np.sum(data['DeathCause'] == 3))

    # Dwell times (same as Evaluation.py)
    all_dwell_prog = []
    all_dwell_fast = []
    for y in range(100):
        row_p = data['DwellTimeProgression'][y, :]
        nz = row_p[row_p != 0]
        all_dwell_prog.extend(nz.tolist())
        row_f = data['DwellTimeFastCancer'][y, :]
        nz = row_f[row_f != 0]
        all_dwell_fast.extend(nz.tolist())

    # Evaluation doubles each value
    doubled_prog = [v for v in all_dwell_prog for _ in range(2)]
    doubled_fast = [v for v in all_dwell_fast for _ in range(2)]
    combined = doubled_prog + doubled_fast

    dwell_prog_med = np.median(doubled_prog) if doubled_prog else 0
    dwell_fast_med = np.median(doubled_fast) if doubled_fast else 0
    dwell_all_med = np.median(combined) if combined else 0

    # Sojourn time
    all_sojourn = []
    for y in range(100):
        row = data['TumorRecord']['Sojourn'][y, :]
        nz = row[row != 0]
        all_sojourn.extend(nz.tolist())
    sojourn_med = np.median(all_sojourn) if all_sojourn else 0

    # Cancer stage distribution from TumorRecord
    all_stages = []
    for y in range(100):
        row = data['TumorRecord']['Stage'][y, :]
        nz = row[row != 0]
        all_stages.extend(nz.tolist())
    total_detected = len(all_stages)
    stage_counts = {7: 0, 8: 0, 9: 0, 10: 0}
    for s in all_stages:
        si = int(s)
        if si in stage_counts:
            stage_counts[si] += 1
    stage_pct = {k: v / total_detected * 100 if total_detected > 0 else 0 for k, v in stage_counts.items()}

    # Detection mode distribution
    all_detect = []
    for y in range(100):
        row = data['TumorRecord']['Detection'][y, :]
        nz_idx = np.nonzero(data['TumorRecord']['Stage'][y, :])[0]
        for i in nz_idx:
            all_detect.append(int(data['TumorRecord']['Detection'][y, i]))
    detect_counts = {1: 0, 2: 0, 3: 0, 4: 0}
    for d in all_detect:
        if d in detect_counts:
            detect_counts[d] += 1
    detect_labels = {1: 'screening', 2: 'symptoms', 3: 'surveillance', 4: 'baseline'}

    # Cancer incidence
    total_prog = int(np.sum(data['ProgressedCancer']))
    total_direct = int(np.sum(data['DirectCancer2']))
    total_fast = len(all_dwell_fast)

    # ================================================================
    # Print results
    # ================================================================
    print()
    print("=" * 72)
    print("SIMULATION RESULTS")
    print("=" * 72)

    print(f"\n--- Demographics ---")
    print(f"  Number of patients:         {n}")
    print(f"  Average age at death:       {avg_age:.2f}")
    print(f"  Average age male:           {avg_age_m:.2f}")
    print(f"  Average age female:         {avg_age_f:.2f}")

    print(f"\n--- Colonoscopies ---")
    print(f"  Screening colonoscopies:    {screening_colo:.0f}")
    print(f"  Symptom colonoscopies:      {symptom_colo:.0f}")
    print(f"  Follow-up colonoscopies:    {followup_colo:.0f}")

    print(f"\n--- Cancer Deaths ---")
    print(f"  Colon cancer deaths:        {ca_deaths}")
    print(f"  Years lost to colon cancer: {yrs_lost_ca:.2f}")
    print(f"  Colonoscopy deaths:         {colo_deaths}")

    print(f"\n--- Cancer Incidence (per {n} patients) ---")
    print(f"  Progressed cancers:         {total_prog}")
    print(f"  Fast cancers:               {total_fast}")
    print(f"  Direct cancers:             {total_direct}")
    print(f"  Total detected cancers:     {total_detected}")

    print(f"\n--- Cancer Stage Distribution ---")
    for stage in [7, 8, 9, 10]:
        label = ['I', 'II', 'III', 'IV'][stage - 7]
        print(f"  Stage {label}:  {stage_counts[stage]:>5}  ({stage_pct[stage]:.1f}%)")

    print(f"\n--- Detection Mode ---")
    for mode, label in detect_labels.items():
        cnt = detect_counts.get(mode, 0)
        pct = cnt / total_detected * 100 if total_detected > 0 else 0
        print(f"  {label:>15}: {cnt:>5}  ({pct:.1f}%)")

    print(f"\n--- Dwell Time (median, Evaluation doubling) ---")
    print(f"  DwellTimeProgressedCa:      {dwell_prog_med:.2f}")
    print(f"  DwellTimeFastCa:            {dwell_fast_med:.2f}")
    print(f"  DwellTimeAllCa:             {dwell_all_med:.2f}")

    print(f"\n--- Sojourn Time ---")
    print(f"  Sojourn time (median):      {sojourn_med:.2f}")

    # ================================================================
    # MATLAB reference comparison
    # ================================================================
    print()
    print("=" * 72)
    print("COMPARISON WITH MATLAB CMOST13 REFERENCE")
    print("=" * 72)
    print(f"{'Metric':<35} {'Python':>10} {'MATLAB ref':>12} {'Status':>8}")
    print("-" * 72)

    def check(name, py_val, matlab_val, tol_pct=30):
        if matlab_val == 0:
            status = "N/A"
        elif abs(py_val - matlab_val) / abs(matlab_val) * 100 < tol_pct:
            status = "OK"
        else:
            status = "DIFF"
        print(f"  {name:<33} {py_val:>10.2f} {matlab_val:>12.2f} {status:>8}")
        return status

    # MATLAB CMOST13 reference values (from prior analysis):
    check("Avg age at death", avg_age, 74.5)
    check("DwellTimeProgressedCa", dwell_prog_med, 13.0, tol_pct=40)
    check("DwellTimeFastCa", dwell_fast_med, 18.0, tol_pct=40)
    check("DwellTimeAllCa", dwell_all_med, 13.0, tol_pct=40)
    check("Sojourn time (median)", sojourn_med, 3.0, tol_pct=40)
    check("Cancer deaths", ca_deaths, 2790, tol_pct=20)
    check("Yrs lost to cancer", yrs_lost_ca, 35652, tol_pct=30)

    # Per 10K normalization for comparison with 10K run
    print()
    print(f"  --- Per 10K patients ---")
    factor = 10000 / n
    check("Progressed Ca /10K", total_prog * factor, 83, tol_pct=50)
    check("Fast Ca /10K", total_fast * factor, 12, tol_pct=50)
    check("Direct Ca /10K", total_direct * factor, 15, tol_pct=50)

    print()
    print("(MATLAB ref values are approximate; tol_pct allows stochastic variation)")
    print()


if __name__ == '__main__':
    main()
