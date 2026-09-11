"""Finalize a validated, visually inspected deliverable; preserve exact bytes in Git."""
from run_transfer import *
import matplotlib,h5py,fitz
from PIL import Image

def main():
 val=json.loads((O/'validation.json').read_text(encoding='utf-8'));assert val['status']=='passed'
 figures=[]
 for png in sorted((O/'figures').glob('*.png')):
  with Image.open(png) as im:width,height=im.size;im.verify()
  pdf=png.with_suffix('.pdf');doc=fitz.open(pdf);assert len(doc)==1;page=doc[0];text=page.get_text();assert 'SV' in text and ('senior' in text or 'transfer' in text);assert page.rect.width>0 and page.rect.height>0
  figures.append(dict(png=str(png.relative_to(O)),pdf=str(pdf.relative_to(O)),png_sha256=sha(png),pdf_sha256=sha(pdf),width_px=width,height_px=height,pdf_pages=1,visual_review='passed: titles/axes/units/five-flow legends readable; final Stage B visually inspected; direction-only comparison confirmed'))
 js('validation/figure_validation.json',dict(status='passed',figures=figures))
 js('metric_columns.json',dict(coordinates='0-based; every end field exclusive; frozen fields are copied verbatim',metric_units={'raw_AUC_200':'SV instrument units * um','vessel_core_P95':'SV instrument units','normalized_AUC_um':'um','normalized_AUC':'um; exact alias'},support_counts={'n_metric_assessable':'count of metric-only lateral/localization assessability flag, independent of Stage B full body score','n_primary_eligible':'metric-only assessability AND (Stage B only) full assessability','n_auc_defined':'primary-eligible AND AUC defined','n_denominator_defined':'primary-eligible AND denominator defined','n_normalized_defined':'primary-eligible AND numerator/denominator defined','n_denominator_pooled':'primary-eligible AND AUC defined AND denominator defined','n_*_defined_all_positions':'metric mathematical validity before full-process score selection'},pooling='framewise endpoint median, never ratio of separate scan medians',zero='valid zero AUC is retained; invalid numeric endpoints are NaN with reason'))
 (O/'.gitattributes').write_text('* -text\n',encoding='ascii')
 audit=(O/'senior_code_audit.md').read_text(encoding='utf-8');audit+='\nVolume coverage exposes metric-only assessability, full-process primary eligibility, and all-position mathematical validity separately. This column clarification changes no scan medians or support.\n';(O/'senior_code_audit.md').write_text(audit,encoding='utf-8')
 prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'));prov['software'].update(matplotlib=matplotlib.__version__,h5py=h5py.__version__);prov.update(analysis_script_sha256=sha(O/'run_transfer.py'),final_script_sha256={p.name:sha(p) for p in sorted(O.glob('*.py'))},contract_sha256=sha(O/'senior_method_contract.json'),finalized_utc=now(),publication_scope='user authorized new branch only; existing frozen outputs unchanged; no merge or force push',diagnostic_archive='validation/stageB_before_loader_dtype_correction_DIAGNOSTIC_ONLY; excluded from final metrics and figures');js('provenance.json',prov)
 val.update(figure_validation='passed',frozen_contract_preserved=True,loader_dtype_correction='passed; no volume median or trend changed',coverage_count_reconstruction='passed',final_scripts_recorded=True);js('validation.json',val)
 js('run_status.json',dict(status='analysis_complete_validated',stageA='success',stageB='success',stageA_volumes=15,stageB_volumes=15,validation='passed',visual_review='passed',blocked=[],missing=[],scope_complete=True,no_extension_planned=True))
 required=['README.md','analysis_plan.json','senior_method_contract.json','provenance.json','run_status.json','input_manifest.csv','validation.json','stageA/volume_metrics.csv','stageB/volume_metrics.csv','comparison/metric_definition_direction_ledger.csv','comparison/stageA_stageB_common_support.csv','comparison/current_vs_senior_definition.csv','validation/fixed_frame_checks.csv','validation/rounding_endpoint_checks.csv','validation/scan_median_reconstruction.csv']
 assert all((O/x).is_file() for x in required)
 for r in read(O/'input_manifest.csv').itertuples():assert sha(R/r.input_path)==r.sha256
 entries=[]
 for p in sorted(O.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.name!='output_sha256.csv':entries.append(dict(file=p.relative_to(O).as_posix(),sha256=sha(p),bytes=p.stat().st_size,role='diagnostic_archive_only' if 'DIAGNOSTIC_ONLY' in str(p) else 'final_deliverable'))
 csv('output_sha256.csv',entries)
 for r in read(O/'output_sha256.csv').itertuples():assert sha(O/r.file)==r.sha256
 print('FINAL PASS',len(entries),'hashed files; total MB',sum(x['bytes'] for x in entries)/1e6,flush=True)
if __name__=='__main__':main()
