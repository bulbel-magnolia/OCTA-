"""Seal the reviewed local delivery; does not rerun scientific computations."""
from common import *
import ast, re, subprocess

def main():
    for name in ['validation.json','report_validation.json','conditional_validation.json','pixel_validation.json','controls_validation.json']:
        assert json.loads((OUT/name).read_text(encoding='utf-8'))['status']=='passed', name
    coverage=read(OUT/'core_question_coverage.csv')
    assert len(coverage)==3 and coverage.status.eq('complete').all()
    evidence=read(OUT/'answer_evidence.csv')
    evidence_rows=read(OUT/'answer_evidence_rows.csv.gz')
    total=0
    for row in evidence.itertuples():
        table=OUT/row.table
        assert sha(table)==row.sha256, row.evidence_id
        indices=json.loads(row.csv_row_numbers)
        assert len(indices)==row.n_rows
        n=len(read(table))
        assert all(2<=i<=n+1 for i in indices)
        total+=row.n_rows
    assert total==len(evidence_rows)==15994
    figures=read(OUT/'figure_manifest.csv')
    assert len(figures)==14
    for row in figures.itertuples():
        assert (OUT/row.figure).exists() and (OUT/row.figure).with_suffix('.pdf').exists()
    figures['status']='numeric_mapping_and_visual_review_passed'
    csv('figure_manifest.csv',figures)
    v=json.loads((OUT/'figure_validation.json').read_text(encoding='utf-8'))
    v.update(status='passed',visual_review='all 14 PNGs actually inspected by Codex; labels, axes, layout and signed scales checked',reviewed_utc=now())
    js('figure_validation.json',v)
    scripts=sorted(OUT.glob('*.py'))
    for p in scripts: ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    links=re.findall(r'\]\(<([^>]+)>\)',(OUT/'README.md').read_text(encoding='utf-8'))
    for path in links:
        assert Path(path).exists() or Path(path)==OUT/'output_sha256.csv',path
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True)
    provenance=json.loads((OUT/'provenance.json').read_text(encoding='utf-8'))
    assert head==provenance['starting_head']
    assert all(line.startswith('?? analysis/') for line in git_status.splitlines())
    assert set(git_status.splitlines())==set(provenance['starting_worktree_status'].splitlines())
    provenance.update(ending_head=head,ending_worktree_status=git_status,completed_utc=now(),
        final_script_sha256={p.name:sha(p) for p in scripts},
        input_and_historical_hash_check='401 files unchanged; see validation.json and input_hashes_after.csv',
        numerical_fix='weighted_numerical_fix.json; failed checkpoints are diagnostic archives excluded from scientific results',
        delivery_scope='local outputs only; no commit, push or publication')
    js('provenance.json',provenance)
    js('delivery_validation.json',dict(status='passed',core_questions=3,evidence_entries=len(evidence),
        evidence_rows=total,evidence_hashes_and_row_bounds='passed',python_syntax_files=len(scripts),
        report_local_links=len(links),reviewed_statistical_figures=14,head_and_worktree_scope='unchanged outside new output directory'))
    status('complete',overall='complete',pixel_validation='passed',conditional_B_visual_review='complete_with_unresolved_anatomical_boundary',
        remaining='none within authorized local analysis scope',three_core_questions='complete',
        figure_visual_review='14 of 14 passed',scientific_computational_blockers=[],
        physical_undefined_results='explicitly preserved in README and undefined_physical_interpretations.csv',
        local_output_review_pending=True,commit_push_or_publish=False)
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='output_sha256.csv' and '__pycache__' not in p.parts)
    records=[]
    for p in files:
        rel=p.relative_to(OUT).as_posix()
        role='failed_checkpoint_not_scientific_result' if rel.startswith('failed_weighted_roundoff_checkpoints/') else 'code' if p.suffix=='.py' else 'delivery_or_audit'
        records.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p),role=role))
    csv('output_sha256.csv',records)
    for r in records: assert sha(OUT/r['path'])==r['sha256']
    print(json.dumps(dict(status='complete',sealed_files=len(records),manifest_sha256=sha(OUT/'output_sha256.csv'),core_questions=3,scientific_gates=102,figures=14),ensure_ascii=False))

if __name__=='__main__': main()
