"""Preserve failed current-run pixel checkpoints before corrected replay."""
from common import *
def main():
    dest=OUT/'failed_weighted_roundoff_checkpoints';assert not dest.exists();dest.mkdir()
    paths=[]
    for pattern in ['framewise_local_background_qc_*.csv.gz','position_integrals_*.csv.gz','geometry_support_D*.csv.gz','shared_pixel_support_D*.csv.gz']:
        paths.extend(OUT.glob(pattern))
    if (OUT/'localization').exists():paths.append(OUT/'localization')
    moves=[]
    for p in paths:
        # Both resolved source and destination stay inside this run's output.
        p.resolve().relative_to(OUT.resolve());target=dest/p.name;target.resolve().relative_to(OUT.resolve());assert not target.exists()
        p.rename(target);moves.append(dict(source=p.name,preserved_as=str(target.relative_to(OUT))))
    js('weighted_numerical_fix.json',dict(failure='NumPy sequential cumulative mass versus Python 3.12 compensated sum disagreed near half weight',
       evidence='weighted_discrepancy_diagnosis.csv',fix='Exact binary rational integer arithmetic when floating error bound overlaps half boundary; independent Fraction reference',
       definition_changed=False,threshold_relaxed=False,formal_signal_changed=False,preserved_failed_checkpoints=moves,rerun_scope='pixel statistics and independent 45-frame checks; table statistics retained',status='corrected_replay_pending'))
    print('Preserved',len(paths),'current-run failed checkpoint paths; no prior v1 or input changes')
if __name__=='__main__':main()
