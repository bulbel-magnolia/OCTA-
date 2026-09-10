"""Run once before scientific statistics; preserve contract and input identities."""
from common import *
import subprocess, scipy, shutil
def main():
    assert not (OUT/'analysis_plan.json').exists(),'Existing analysis contract: do not overwrite'
    allowed={'common.py','initialize.py'}
    assert all(p.name in allowed for p in OUT.iterdir()),'Conflicting output files'
    text=Path('C:/Users/zby/Downloads/SV_ROI_v2_Codex_实施提示词.md').read_text(encoding='utf-8')
    (OUT/'analysis_specification.md').write_text(text,encoding='utf-8')
    plan=dict(version='sv_fixed_diameter_selective_v2',frozen_utc=now(),informed_by_historical_results=True,independent_validation=False,
        specification_file='analysis_specification.md',specification_sha256=sha(OUT/'analysis_specification.md'),
        core_questions=['matched-flow absolute source/tail and relative tail/source','within-volume depth association and zero-lag specificity','matched local background, spatial controls and localization sensitivity'],
        cohort=dict(diameters_um=[128,235,285],flows_mm_s=[1,3,5,7,10],volumes=15,nominal_frames=7500,valid_frames=7018,frame_index=[0,499]),
        excluded=dict(diameter_um=500,volumes=8,timing='post-hoc after historical image/QC review; no arrays read'),
        signal='SV_raw=var(abs(IMG),1,3); denominator N; linear; mean=Q/area; Q=sum(SV*w)*dx*dz',dx_um=DX,dz_um=DZ,
        geometry='frozen X4/z_top; physical D circle; lower boundary z_top+D/dz is prior; 16x16 source lattice; 6 bins; tail 20x25um; guard=0',
        primary_scalar='RI500 = per-frame tail500 mean / source mean, followed by volume median; zero denominator NaN',
        other_scalars=['RI100','source mean','tail100 mean','tail500 mean','20-band tail mean and per-frame RI profile'],
        C16=C,C16_rule='median of all 16 signed lag0 rhos; NaN if any undefined',representatives=REP,
        stats=dict(correlation='pairwise average ranks followed by Pearson; >=3 nonconstant pairs; no p-values',detrend_windows=[51,31,101],
                   trend='centered moving median on original coordinates; min_periods=1; ignore neighborhood NaN, preserve missing center; endpoints truncate',
                   lag=list(range(-50,51)),lag_definition='Source(i) versus Tail(i+lag), non-circular',far_lags=LAGS[FAR],
                   specificity='rho0 minus median of fixed 42 far-lag rhos; all required',peak='maximum signed rho; tie min abs(lag), then smaller signed lag; all 101 required',
                   rank='1+count(rho_lag>rho0); all 101 required',matched_specificity='median of 42 coordinate-matched rho0-rho_lag deltas; all required',
                   differences='raw only, original adjacent i and i+1 both present'),
        controls=dict(offset_D=1.5,combined='both complete required; concatenate pixels and weights',common_support='real/L/R intersect before separate detrending',
                      Z='(weighted_median_real-weighted_median_BG)/(1.4826*weighted_MAD_BG); no epsilon, no thresholds',
                      weighted_median='cumulative mass reaches half; exact half averages adjacent values',
                      excess='cellwise real_matched-(rho_L+rho_R)/2; strict 16-cell median; separate signed real-L and real-R'),
        position_offsets_px=OFFSETS,position_support='baseline plus all four offsets intersect before detrending',
        sensitivity_scope='31/101 representatives full lag and C16 lag0 only; position source 6 bins, tail T1-T4, tail100/500 only',
        conditional_extensions='A-E exactly as specification; log trigger and fixed scope before execution; do not change primary contract',
        endpoint_drift='indices 0..50 and 449..499 inclusive, 51 positions each; signed difference last-first; absolute magnitude also stored; relative difference divided by whole-curve median',
        validation=dict(integration_relative_tolerance=1e-10,zero_reference='exact absolute zero check',rho_absolute_tolerance=1e-12),
        settings_consistency='unverified; raw units instrument units',slow_axis_step='unknown; frame units only',independent_unit='scan volume',
        final_scope='local outputs for user review; no commit/push/publication')
    js('analysis_plan.json',plan)
    mapping={r:{m:col(r,m) for m in ['mean','q','area']} for r in REGIONS if r!='s36'}
    mapping['s36']={'mean':'sum(source_s3..s6_q_raw)/sum(source_s3..s6_area_um2)','q':'sum(source_s3..s6_q_raw)','area':'sum(source_s3..s6_area_um2)'}
    js('column_mapping.json',dict(input=str(INPUT),join_keys=KEY,geometry_backbone='geometry_manifest.csv.gz',regions=mapping))
    before=[]
    for folder in [INPUT,PREV,DEP,QC,ROOT/'analysis/sv_diameter_stage2_scientific_v1']:
        for p in sorted(folder.rglob('*')):
            if p.is_file() and '__pycache__' not in str(p):before.append(dict(input_path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,role='read_only_input_or_history'))
    p=ROOT/'analysis/ACTIVE_ANALYSIS_POLICY.json';before.append(dict(input_path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,role='active_policy'))
    for r in read(INPUT/'output_sha256.csv').itertuples():assert sha(INPUT/r.file)==r.sha256,r.file
    csv('input_manifest.csv',before)
    fw=frame_table();assert fw.groupby('diameter_um').valid_geometry.sum().to_dict()=={128:2422,235:2225,285:2371}
    spots=read(INPUT/'independent_integration_spotchecks.csv')[KEY].drop_duplicates();assert len(spots)==45
    spots['selection']='recovered existing geometry-v2 fixed first/nearest-original-249/last IDs; no result selection'
    csv('localization_fixed_frames.csv',spots)
    identities=read(INPUT/'array_identity_audit.csv.gz').fillna('');assert len(identities)==7018
    absent=[p for p in identities.input_path.unique() if not Path(p).exists()]
    js('preflight.json',dict(status='passed' if not absent else 'pixel_inputs_blocked',array_paths_checked=len(identities.input_path.unique()),missing_paths=absent,
                            input_output_hashes_passed=True,nominal=7500,valid=7018,output_preexisted_as_results=False))
    js('provenance.json',dict(starting_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
        starting_worktree_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        contract_sha256=sha(OUT/'analysis_plan.json'),software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
        input_manifest_sha256=sha(OUT/'input_manifest.csv'),analysis_guard='D:/research-skills/skills/analysis-guard/SKILL.md',
        authorization='current user requests complete, definition-consistent, traceable v2 answers; attachment supplies analysis specification, not independent authority'))
    (OUT/'selection_log.md').write_text('# Selection log\n\nThe current user requests complete v2 answers to all three questions. The attached draft supplies the scientific contract. All parameters were frozen before new scientific statistics.\n\nHistorical same-data observations informed C16, the allowed sensitivity set and the D500 exclusion. This is not preregistration or independent validation. Every flow, signed cell, deep band and undefined result remains represented. No old background or coupling is relabeled as v2.\n\nThe existing 45 geometry-v2 verification frame IDs are recovered verbatim (first, nearest original 249, last); see localization_fixed_frames.csv. Endpoint drift uses original 0..50 and 449..499 inclusive, matching the prior explicit coordinate convention; relative drift divides by the full-curve median.\n',encoding='utf-8')
    status('contract_frozen',overall='running',pixel_inputs_available=not absent)
    print('Contract frozen; input hashes and 7500-position identity passed; missing retained paths:',len(absent),flush=True)
if __name__=='__main__':main()
