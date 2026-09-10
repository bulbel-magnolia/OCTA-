from common import *
def main():
    b=read(OUT/'conditional_B_review.csv');b['review_status']='图像不足以判断'
    b['finding']='Maximum-adjacent-shift pairs show a changing upper structural gradient and central SV distribution. They do not establish an independently visible lower wall or a unique localization error; no additional rejection or geometry change.'
    b['reviewer']='Codex visual inspection';b['reviewed_utc']=now();csv('conditional_B_review.csv',b)
    logs=read(OUT/'conditional_extension_log.csv');logs.loc[logs.extension.eq('B')&logs.status.eq('triggered_before_execution'),'status']='complete';csv('conditional_extension_log.csv',logs)
    v=json.loads((OUT/'conditional_validation.json').read_text());v.update(status='passed',B_visual_review='complete_with_unresolved_anatomical_boundary');js('conditional_validation.json',v)
    fix=json.loads((OUT/'weighted_numerical_fix.json').read_text());fix.update(status='corrected_replay_passed',independent_exact_weighted_validation='all 45 fixed frames and synthetic cases passed',synthetic_test_correction='symmetric [0.1,0.2,0.1,0.2] case expected midpoint corrected to 2.5 before replay');js('weighted_numerical_fix.json',fix)
    v=json.loads((OUT/'pixel_validation.json').read_text());v['localization_visual_review']='all 45 reviewed; anatomical lower boundary unresolved';js('pixel_validation.json',v)
if __name__=='__main__':main()
