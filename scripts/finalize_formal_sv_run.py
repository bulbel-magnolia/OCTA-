"""Finalize the approved full2500 run after all checkpoint workers finish."""
from pathlib import Path
import sys,json,subprocess,html
import pandas as pd
import numpy as np
sys.path.insert(0,'src')
from svrecttail.collection import sha256,json_write,TRACKING
p=Path('results/formal_sv_d128_v21_full2500_run001');bridge=p.parent/'formal_sv_d128_v21_bridge15_run001'
t=pd.read_csv(p/'frame_results.csv',float_precision='round_trip');summary=pd.read_csv(p/'scan_summary.csv');sel=pd.read_csv(p/'preselected_qc_manifest.csv')
assert len(t)==2500 and not t.duplicated(['scan_id','frame_index_0based']).any()
assert t.input_mapping_valid.all()
for r in sel.itertuples():assert (p/r.mapping_image).is_file()
comparisons=[]
for scan,g in t.groupby('scan_id'):
 tr=pd.read_csv(TRACKING/f'tracking/{scan}/{scan}_mentor_tracking.csv')
 joined=g.merge(tr,left_on=['scan_id','frame_index_0based'],right_on=['scan_id','frame_index'],suffixes=('_run','_frozen'),validate='one_to_one')
 for col in ['z_upper_px','seed_z_upper_px','x4_centroid_isolated_jump_corrected_px','x1_local_geometry_px','local_body_run_width_px']:
  np.testing.assert_array_equal(joined[col+'_run'],joined[col+'_frozen'])
 comparisons.append({'scan_id':scan,'rows':len(g),'frozen_coordinates_identical':True})
small=pd.read_csv(bridge/'frame_results.csv');match=small.merge(t,on=['scan_id','frame_index_0based'],suffixes=('_bridge','_full'),validate='one_to_one')
for col in ['q_vessel','q_tail','ratio_tail_to_vessel','z_upper_px']:
 np.testing.assert_allclose(match[col+'_bridge'],match[col+'_full'],rtol=1e-12,equal_nan=True)
json_write(p/'frozen_mapping_validation.json',{'scans':comparisons,'bridge15_metrics_match_full_run':True,'no_localization_rerun':True})
status=json.loads((p/'execution_status.json').read_text());status['coordinate_mapping_review']='stage_A_all_15_passed; full_preselected_diagnostics_and_review_images_available';status['fixed_qc_positions']=len(sel);json_write(p/'execution_status.json',status)
env=json.loads((bridge/'software_versions.json').read_text());env['python_tests']='59 passed; one unrelated requests dependency-version warning';env['frame_processing_code_commits']=sorted(t.code_commit.unique());env['final_verification_code_commit']=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();json_write(p/'software_versions.json',env)
files=list((p/'tracking').rglob('*'))+[p/'run_config.json',p/'tracking_config.json']
pd.DataFrame([{'file':f.relative_to(p).as_posix(),'sha256':sha256(f),'bytes':f.stat().st_size} for f in files if f.is_file()]).to_csv(p/'frozen_inputs_sha256.csv',index=False)
excluded=t.loc[~t.valid].copy();excluded['detail']=excluded.invalid_reason+';'+excluded.localization_invalid_reason.fillna('')
excluded.groupby(['scan_id','detail']).size().reset_index(name='frames').to_csv(p/'exclusion_detail_counts.csv',index=False)
view=['scan_id','planned_frames','z_present','geometry_qc_valid','quantification_valid','quantification_excluded','direct_supported_included','assisted_included']
header='| 扫描卷 | 计划 | 有 z | 几何通过 | 定量有效 | 排除 | 纳入直接支持 | 纳入辅助定位 |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
rows=''.join('| '+' | '.join(map(str,r))+' |\n' for r in summary[view].itertuples(index=False,name=None))
intro=f'''# 128 μm v2.1：五卷完整 SV 数据

已完成 5×500=2500 帧的真实 OCT 重建、固定区域定量和二维归档。定量有效 {int(t.valid.sum())} 帧，排除 {int((~t.valid).sum())} 帧；全部计划帧均在主表和数组包中。检出深度未运行，原因是没有匹配空白。

{header}{rows}
直接候选支持与短缺口补全的总数、纳入数分别见 `scan_summary.csv`。定位总覆盖率与定量有效率分别计算，分母始终为计划的 2500 帧。`exclusions.csv` 保留全部排除帧；`exclusion_detail_counts.csv` 将 source QC 失败展开到原定位原因。

'''
base=(bridge/'README.md').read_text(encoding='utf-8')
method=base[base.index('定位、几何、背景'):base.index('所有 15 张')]
review='''原 15 帧接入判据已通过，见同级 bridge15 结果目录。全卷与首包的原 15 帧指标逐项一致，2500 行冻结坐标与原定位表一对一核对通过，见 `frozen_mapping_validation.json`。`archive_validation.json` 保存独立从 NPZ 二维矩阵重算积分和背景的检查结果。`preselected_qc_manifest.csv` 列出全部 91 个固定 QC 位置，包括原 15 帧、均匀位置、flow03 原错误分支与长缺口/重锁定邻域。

打开 `qc_gallery.html` 查看固定图及每帧纳入状态。有几何但未通过 QC 的帧仍画模板用于诊断，画出椭圆并不表示纳入分析。已复看 flow03 第 420、423、432、434、436、447、469–471 帧的同帧叠加；没有新增坐标映射错误。第 432、447 帧因冻结 assessability 不通过被排除。此轮没有重估定位准确率或调整核心追踪。

所有输入信号均来自该帧自身的三次重复，辅助定位只影响坐标。无匹配空白，`detected` 和长度保留 NA；不存在“长度为零”的未运行结果。实验身份、独立重复、时间及慢轴物理标定未知处保留 NA；每卷 500 个空间位置不作为 500 次独立实验。该批支持空间描述，未进行显著性检验或参数筛选。

大数组按每卷每 100 帧一个 ZIP，共 25 包。每包含该段的完整二维数据、相应 profile/QC 和共同表格/冻结配置；全批表格为全 2500 行，数组只含包名对应段。Release 入口为 [formal-sv-d128-v21-run001](https://github.com/bulbel-magnolia/OCTA-/releases/tag/formal-sv-d128-v21-run001)。逐包名称、大小和 SHA-256 在 `download_packages.csv`，逐帧在 `arrays_sha256.csv`。首包 bridge15 是另一个自包含包，不替换全卷数据。

`output_sha256.csv` 和 `frozen_inputs_sha256.csv` 对应 ZIP 内实际字节；Git 可能规范化文本换行，Git 浏览文件的字节哈希因此可能不同。冻结算法同时记录实际字节和 Git LF 哈希。`collection_run.json` 是准备阶段的冻结快照，最终执行状态以 `execution_status.json` 和验证报告为准。原 OCT/DICOM/MAT 留本地，不上传原始采集文件。
'''
(p/'README.md').write_text(intro+method+review,encoding='utf-8')
body=['<!doctype html><meta charset="utf-8"><title>SV v2.1 fixed QC</title><style>body{font:16px sans-serif;margin:28px}img{width:100%;max-width:1500px}.excluded{color:#b11}article{margin-bottom:32px}</style><h1>SV v2.1 fixed QC</h1><p>Shapes show frozen geometry. Only valid=true frames enter formal summaries. Pink: both background strips; cyan: full source; yellow: 500 μm tail.</p>']
indexed=t.set_index(['scan_id','frame_index_0based'])
for r in sel.itertuples():
 rec=indexed.loc[(r.scan_id,r.frame_index_0based)];flag='INCLUDED' if rec.valid else 'EXCLUDED'
 body.append(f'<article><h2 class="{"" if rec.valid else "excluded"}">{r.scan_id} / {r.frame_index_0based}: {flag}</h2><p>{html.escape(str(rec.invalid_reason))}; {html.escape(str(rec.localization_source))}; z={rec.z_upper_px}</p><img loading="lazy" src="{r.mapping_image}"></article>')
(p/'qc_gallery.html').write_text('\n'.join(body),encoding='utf-8')
files=[f for f in p.rglob('*') if f.is_file() and f.name not in ['output_sha256.csv','download_packages.csv']]
pd.DataFrame([{'file':f.relative_to(p).as_posix(),'bytes':f.stat().st_size,'sha256':sha256(f)} for f in sorted(files)]).to_csv(p/'output_sha256.csv',index=False)
print(summary[view].to_string(index=False))

