"""Verify new cohort/geometry tables and publish the local active analysis policy."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def table(d):return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None)])
def main():
    val=json.loads((OUT/'validation.json').read_text(encoding='utf-8'));assert val['status']=='passed'
    frames=pd.concat([pd.read_csv(p,float_precision='round_trip') for p in sorted(OUT.glob('framewise_fixed_diameter_D*.csv.gz'))],ignore_index=True)
    assert set(frames.diameter_um)=={128,235,285} and frames.scan_id.nunique()==15 and len(frames)==7500
    valid=frames[frames.valid_geometry];invalid=frames[~frames.valid_geometry]
    assert len(valid)==7018 and len(invalid)==482
    assert invalid[['X4','z_top','source_mean_raw','tail_mean_raw_100um','tail_mean_raw_500um']].isna().all().all()
    assert np.allclose((valid.x_right_edge_px-valid.x_left_edge_px)*12.7,valid.diameter_um,rtol=0,atol=1e-9)
    assert np.allclose((valid.z_bottom_edge_px-valid.z_top)*6.7,valid.diameter_um,rtol=0,atol=1e-9)
    original=pd.read_csv(ROOT/'analysis/sv_source_tail_framewise_coupling_v1/framewise_source_tail_metrics.csv',float_precision='round_trip')
    joined=valid.merge(original,on=['scan_id','frame_index'],validate='one_to_one',suffixes=('','_old'))
    assert np.array_equal(joined.X4,joined.X4_old) and np.array_equal(joined.z_top,joined.z_top_old)
    assert np.array_equal(joined.original_apparent_width_um,joined.X1)
    exclusion=pd.read_csv(OUT/'d500_exclusion_record.csv');assert len(exclusion)==8 and exclusion.diameter_um.eq(500).all()
    cover=pd.read_csv(OUT/'coverage.csv');changes=pd.read_csv(OUT/'volume_geometry_change_summary.csv')
    width=valid.groupby('diameter_um').original_apparent_width_um.agg(['median','min','max']).reset_index()
    width['new_width_um']=width.diameter_um
    policy=json.loads((OUT/'analysis_policy.json').read_text(encoding='utf-8'))
    policy.update(status='validated_geometry_and_raw_integration_inputs',active_input_directory='analysis/sv_physical_diameter_geometry_v2',
      primary_volumes=15,primary_valid_frames=7018,primary_nominal_frames=7500,
      provenance='user-approved change on 2026-09-10; prior versions retained',
      excluded_volume_ids=exclusion.scan_id.tolist(),excluded_valid_frames=3913,
      downstream_rule='New scientific analyses must use this fixed-physical-D geometry and this 15-volume cohort; never combine old-X1 and new-D metrics under one method label.')
    js(ROOT/'analysis/ACTIVE_ANALYSIS_POLICY.json',policy)
    text='''# 固定真实直径 ROI v2：已批准的几何与积分输入更新

2026-09-10 用户批准将 source 改为固定真实直径，并明确选择“D500 不再进入后续科学分析，但保留数据及排除记录”。本目录已落实该决定；分析入口政策为 `analysis/ACTIVE_ANALYSIS_POLICY.json`。

## 新版范围

D128、D235、D285，各共同流速1/3/5/7/10 mm/s，共15卷、7,018个原有效帧。7,500个名义位置全部保留，482个原无效geometry位置仍为缺失。没有依据新结果恢复或新增筛除任何帧。

D500全部8卷（含2/9/12 mm/s扩展）不进入后续科学分析；原始数据、3,913个有效帧的已有结果、QC、截图及排除记录全部保留。排除决定发生在查看已有图像、QC及coupling之后，属于事后分析范围调整，不是预注册排除，也不表示仪器采集失败已被证实。后续结论的直径范围限于128–285 μm。

## 已执行的几何更新

- 保留每个原有效帧的冻结X4和z_top。
- Source物理宽度=真实D、物理高度=真实D，圆形截面先验；像素坐标中的椭圆宽D/12.7、高D/6.7，不取整。
- 圆心z=z_top+D/(2×6.7)，左右边界X4±D/(2×12.7)。下界仍为z_top+D/6.7。
- Tail横向宽度同步为真实D，从同一物理下界开始，guard=0，无扩散锥；保留0–100/0–500 μm整体窗口和20个25 μm bands。
- Source保留6个归一化轴向bins和16×16 fractional weighting。
- X1原值改以original_apparent_width_um保存，作为原算法的表观宽度记录；不将其冒充真实直径。
- 正式SV仍为var(abs(IMG),1,3)，分母N，linear raw SV，无gain/log/normalization/background subtraction。

## 完成范围与后续边界

已从同身份raw arrays完成15卷新版source/tail mean/Q/area与6×20深度输入的重新积分，另保存逐卷新旧差值和有效像素支持。旧版两个coupling目录、旧QC及原几何结果均未改写。

本次没有重跑新的coupling或新的matched-background QC。旧coupling与旧QC属于X1几何，不能贴上v2标签继续使用；后续科学分析应读取本目录输入并按v2重算。旧版本脚本为复现旧结果而保留，不应直接作为新版科学分析入口。固定D也不自动修正上边界或中心的潜在误差。

## 为什么记录D500为范围调整

此前旧几何QC中，D500五个共同流速的S1 robust Z中位数为6.428，S6为0.612，whole-source为1.516，proximal tail为0.258；上部与深部可检测性存在明显差异。该事实与上部亮信号集中的观察相容，不等同于只有上部有血流，不能证明整组数据无效。用户决定将其排除于后续分析；排除原因和时点在cohort_manifest及d500_exclusion_record中公开记录。

## 覆盖与宽度

'''+table(cover)+'\n\n'+table(width)+'''\n
## 验证与复现

所有7,018个原有效帧的数组SHA256与正式输入逐一一致，D128所用25个release ZIP验证通过。先用同一数组重放原source/tail mean/Q核对输入，再积分新版ROI。新旧X4、z_top、物理bottom一致，source/tail均完整落在FOV内，新增无效帧=0。45个预先固定帧（每卷首、中、末有效帧）的全图source权重独立对照通过；各项数值误差见validation.json。

```powershell
python analysis/sv_physical_diameter_geometry_v2/build_fixed_diameter_inputs.py
python analysis/sv_physical_diameter_geometry_v2/finalize_geometry_migration.py
```

主要文件：analysis_policy.json与全局ACTIVE_ANALYSIS_POLICY.json；cohort_manifest.csv；d500_exclusion_record.csv；geometry_manifest.csv.gz；三个framewise_fixed_diameter_D*.csv.gz；volume_geometry_change_summary.csv；array_identity_audit.csv.gz；input_manifest.csv；validation.json；provenance.json；output_sha256.csv。
'''
    (OUT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    val.update(independent_table_geometry_check='passed',independent_frozen_center_top_check='passed',invalid_metric_missingness_check='passed',active_policy_written=True)
    js(OUT/'validation.json',val)
    prov=json.loads((OUT/'provenance.json').read_text(encoding='utf-8'));prov['delivered_scripts_sha256']={p.name:sha(p) for p in OUT.glob('*.py')};prov['active_policy_sha256']=sha(ROOT/'analysis/ACTIVE_ANALYSIS_POLICY.json');js(OUT/'provenance.json',prov)
    files=[p for p in OUT.iterdir() if p.is_file() and p.name!='output_sha256.csv']
    pd.DataFrame([dict(file=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)]).to_csv(OUT/'output_sha256.csv',index=False,lineterminator='\n')
    print('Validated active policy: 15 volumes / 7018 valid frames; D500 excluded and retained.')
    print(changes.groupby(['diameter_um','metric']).paired_frame_change_percent_median.median().to_string())
if __name__=='__main__':main()
