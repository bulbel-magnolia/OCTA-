"""Record actual assistant visual assessment of all 45 fixed panels."""
from common import *
def main():
    review=read(OUT/'localization_review.csv');logs=read(OUT/'conditional_extension_log.csv').to_dict('records')
    review['figure']=review.scan_id.map(lambda s:f'localization_review_display/{s}.png')
    review['review_status']='图像不足以判断'
    review['review_basis']='Retained structural and linear raw-SV panels inspected; incomplete structural contour prevents validation of full physical-D circle or lower wall. SV signal location is descriptive, not truth.'
    review['structural_gradient']=True;review['clear_misalignment']=False;review['single_side_background_mismatch']=False
    review['note']=review.diameter_um.map({128:'Narrow central SV column lies within source and proximal tail; structural upper-to-lower brightness gradient and incomplete lower contour.',235:'Central low-structural-signal region overlaps source; SV distribution within circle changes from first to later frames; lower structural boundary not independently resolved.',285:'Central dark structural region overlaps source, with an axial gradient; SV brighter in upper part in later frames. Full lower contour remains unresolved.'})
    hit=(review.diameter_um==285)&review.flow_mm_s.isin([1,3])&review.frame_index.eq(0)
    review.loc[hit,'note']='Frame 0: bright structural feature at left source/lower-boundary vicinity extends toward proximal tail; fixed panels insufficient to distinguish interface crossing from local feature. Trigger B; no frame exclusion or geometry adjustment.'
    review['reviewer']='Codex visual inspection of rendered retained arrays';review['reviewed_utc']=now()
    csv('localization_review.csv',review)
    logs=[r for r in logs if r.get('extension') not in ['B','D','E']]
    for scan,g in review.groupby('scan_id'):
        meta={k:g[k].iloc[0] for k in ID};logs.append(meta|dict(extension='D',trigger='fixed panels show axial structural brightness gradient intersecting source/tail whose relation changes with retained ROI position',evidence_table='localization_review.csv',evidence='fixed frames '+','.join(map(str,g.frame_index)),scope='whole-source/tail100/tail500 structural means real/L/R, same-ROI structural-SV and structural source-tail raw/51 correlations only',status='triggered_before_execution'))
        b=scan in ['D285_F01_V01','D285_F03_V01'];logs.append(meta|dict(extension='B',trigger='frame 0 structural feature near source lower-left/tail origin; fixed three panels insufficient' if b else 'no clear misalignment or qualifying interface crossing beyond gradient',evidence_table='localization_review.csv',evidence='frame 0' if b else 'three fixed panels',scope='maximum adjacent |delta X4| and |delta z_top| frame pairs, earliest ties; no geometry changes',status='triggered_before_execution' if b else 'not_triggered'))
        logs.append(meta|dict(extension='C',trigger='neither signed C16 side-dependence nor clear single-side structural mismatch observed',evidence_table='fixed_region_control_excess.csv; localization_review.csv',evidence='real-minus-L and real-minus-R both positive; same-depth field gradient recorded under D',scope='no +/-2.5D extension; retain primary +/-1.5D',status='not_triggered'))
        logs.append(meta|dict(extension='E',trigger='report describes fixed broad region and all six rows; does not rely on exact single-bin location',evidence_table='fixed_region_summary.csv; position_sensitivity_cells.csv',evidence='no exact-bin localization claim required',scope='no 6-to-3 bin extension',status='not_triggered'))
    csv('conditional_extension_log.csv',logs)
    fixed=[];j=read(OUT/'geometry_jitter_adjacent.csv.gz')
    for scan in ['D285_F01_V01','D285_F03_V01']:
        for parameter in ['X4','z_top']:
            a=j[(j.scan_id==scan)&(j.parameter==parameter)].sort_values(['absolute_delta_px','frame_i'],ascending=[False,True]).iloc[0]
            fixed.append(dict(scan_id=scan,parameter=parameter,frame_i=int(a.frame_i),frame_next=int(a.frame_next),absolute_delta_px=a.absolute_delta_px,selection='maximum actual adjacent absolute change; earliest tie'))
    csv('conditional_B_fixed_pairs.csv',fixed)
    status('visual_review_complete',conditional_extensions='D all 15; B two scans; A complete; C and E not triggered')
if __name__=='__main__':main()
