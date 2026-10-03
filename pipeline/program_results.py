"""Generate membership from programs.py and canonical reporting tables from R outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from programs import PROGRAMS


def export_membership(destination):
    destination.write_text(json.dumps(PROGRAMS,indent=2)+'\n')


def export_results(output):
    analysis=output/'program_analysis'
    results=pd.read_csv(analysis/'combined_results.csv',float_precision='round_trip')
    scores=pd.read_csv(output/'module_scores.csv',index_col=0,float_precision='round_trip')
    primary=[];sensitivity=[]
    for model,frame in results.groupby('Model',sort=False):
        for row in frame.to_dict('records'):
            name=row['Program']
            common={'program':name,'n_genes':row['NGenes'],'model':'~ W_1 + genotype' if model=='RUV_adjusted' else '~ genotype',
                    'contrast':'dcKO minus CON','n_libraries_per_group':3,
                    'residual_df':3 if model=='RUV_adjusted' else 4}
            seed=20261003 if model=='RUV_adjusted' else 20261004
            if model=='RUV_adjusted':
                con=scores.loc[scores.group=='CON',name].mean();ko=scores.loc[scores.group=='dcKO',name].mean()
                primary.append({**common,'mean_CON':con,'mean_dcKO':ko,'delta_score':ko-con,
                    'score_interpretation':'Descriptive mean gene-wise z score; not the modeled ROAST statistic',
                    'direction':row['ROAST_direction'],'ROAST_P':row['ROAST_P'],'ROAST_BH_adjusted_P':row['ROAST_FDR'],
                    'method':'limma::mroast mean; two-sided directional PValue','rotations':99999,'seed':seed,
                    'adjustment_family':'ROAST_directional_RUV_adjusted_four_programs','midp':False})
            for method,prefix in [('ROAST','ROAST'),('ROAST_mixed','ROAST_mixed')]:
                if model=='RUV_adjusted' and method=='ROAST':continue
                sensitivity.append({**common,'method':method,'direction':row['ROAST_direction'] if method=='ROAST' else 'Mixed',
                    'raw_P':row[prefix+'_P'],'BH_adjusted_P':row[prefix+'_FDR'],
                    'adjustment_family':f'{method}_{model}_four_programs',
                    'rotations':99999,'seed':seed,
                    'role':'supplementary sensitivity' if method!='ROAST_mixed' else 'secondary mixed-direction endpoint'})
    for name,rows in [('Supplementary_Data_S3_program_tests.csv',primary),('Supplementary_Data_S3_program_test_sensitivity.csv',sensitivity)]:
        pd.DataFrame(rows).to_csv(output/name,index=False,float_format='%.17g')
    provenance=json.loads((analysis/'input_provenance.json').read_text())
    provenance['membership_authority_sha256']=hashlib.sha256(Path(__file__).with_name('programs.py').read_bytes()).hexdigest()
    (analysis/'input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--membership',type=Path);ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    if args.membership:export_membership(args.membership)
    if args.output:export_results(args.output)
