"""Check delivery identities and write a complete output SHA256 inventory."""
import json,subprocess
from derive_qc import OUT,ROOT,digest,csv,read
import pandas as pd
def main():
    p=json.loads((OUT/'provenance.json').read_text());v=json.loads((OUT/'validation.json').read_text())
    assert v['status']=='passed' and not v['unfinished_computation'] and not v['missing_required_data']
    assert digest(OUT/'input_manifest.csv')==p['input_manifest_sha256']
    for name,h in p['delivered_script_sha256'].items():assert digest(OUT/name)==h
    before=json.loads((OUT/'frozen_files_before.json').read_text())
    assert all(digest(ROOT/name)==h for name,h in before.items())
    changes=subprocess.check_output(['git','diff','--cached','--name-status'],cwd=ROOT,text=True).splitlines()
    assert all(line.startswith('A\tanalysis/sv_data_quality_qc_v1/') for line in changes)
    required=['README.md','input_manifest.csv','provenance.json','validation.json','volume_qc_summary.csv',
      'volume_background_stability.csv','volume_geometry_qc.csv','background_roi_coverage.csv',
      'source_depth_detectability_summary.csv','tail_depth_detectability_summary.csv','proximal_tail_detectability_summary.csv',
      'background_real_tail_correlation.csv','pseudo_source_tail_coupling_summary.csv','control_excess_coupling_summary.csv',
      'negative_control_volume_summary.csv','qc_vs_coupling_summary.csv','d500_quality_vs_other_diameters_summary.csv','structural_qc_summary.csv']
    assert all((OUT/x).is_file() for x in required)
    assert len(list(OUT.glob('figure*.png')))==10 and len(list(OUT.glob('figure*.pdf')))==10
    files=[f for f in OUT.iterdir() if f.is_file() and not f.name.startswith('checkpoint_') and f.name not in ['extraction_progress.json','output_sha256.csv']]
    assert all(f.stat().st_size<100_000_000 for f in files)
    csv('output_sha256.csv',pd.DataFrame([dict(file=f.name,sha256=digest(f),bytes=f.stat().st_size) for f in sorted(files)]))
    print('Delivery passed:',len(files),'hashed files; total',sum(f.stat().st_size for f in files),'bytes; only QC additions staged')
if __name__=='__main__':main()
