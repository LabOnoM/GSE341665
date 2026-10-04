#!/usr/bin/env python3
"""Draw the author-confirmed September experimental allocation (Figure S1).

Authority: September Methods and WeChat notes, ingested as
.wiki/sources/september_2026_experimental_update.md. The superseded July generator
is preserved in _revision_history/2026-09-24/before_update.zip.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.path import Path as MplPath

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'scripts/JDR_reanalysis/figures/FigureS1_breeding_scheme'

# Explicit design metadata audited independently of drawing placement.
DESIGN={
 'constitutive':{'driver':'K14-cre','age':'E16.5','sex':'male','4OHT':False,'figures':['1']},
 'inducible':{'driver':'K14-creER','age':'P0','sex':'male','4OHT':True,'figures':['2','3','4','5','6','7','S7','S4','S2']},
 'transplant':{'recovery_days':14,'recipient_age_weeks':8,'paired_recipients':3},
}

def main():
 assert DESIGN['constitutive']['figures']==['1'] and not DESIGN['constitutive']['4OHT']
 assert '3' in DESIGN['inducible']['figures'] and DESIGN['inducible']['4OHT']
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
 # Match the author's displayed figure; editorial notes belong outside the artwork.
 y_min,y_max=2.440073,7.332069
 fig,ax=plt.subplots(figsize=(7.2,y_max-y_min));fig.subplots_adjust(0,0,1,1);ax.set_xlim(0,7.2);ax.set_ylim(y_min,y_max);ax.axis('off')
 def box(x,y,w,h,txt,color='#EEF2F5',size=10):
  ax.add_patch(Rectangle((x,y),w,h,facecolor=color,edgecolor='#708090',lw=.8))
  ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=size,linespacing=1.3)
 def arrow(x,y0,y1):ax.add_patch(FancyArrowPatch((x,y0),(x,y1),arrowstyle='-|>',mutation_scale=12,color='#555555',lw=.9))
 box(.25,6.57,3.15,.68,'Constitutive K14-cre × Smad4 fl/fl\nCre-negative littermate controls','#DCE6F4')
 box(3.8,6.57,3.15,.68,'Inducible K14-creER × Smad4 fl/fl\nCre-negative littermate controls','#FAEBDD')
 arrow(1.825,6.57,6.31);arrow(5.375,6.57,6.31)
 box(.25,5.48,3.15,.83,'E16.5 male first-molar donors\n24-h culture; no 4-OHT\nConstitutive model: Figure 1')
 box(3.8,5.48,3.15,.83,'P0 male first-molar donors\n90 µmol/L 4-OHT for 24 h\nBoth genotypes treated identically')
 arrow(1.825,5.48,5.20);arrow(5.375,5.48,5.20)
 box(.25,3.88,3.15,1.32,'Subrenal-capsule transplantation\n8-week-old male recipients\nOne mutant + one control per kidney\nRecovery after 14 days\nThree recipients; six grafts')
 box(3.8,3.88,3.15,1.32,'Mechanical epithelial separation\nRT-qPCR / immunoblot: Figures 2–3\nRNA sequencing: Figures 4–7\nInhibitor-response assays: Figure 3\nTissue-region image: Figure S2',size=9.5)
 arrow(1.825,3.88,3.62)
 box(.25,2.90,3.15,.72,'Figure 1\nRadiography, H&E, AMBN staining','#DCE6F4')
 # Branch directly from treated tooth germs, bypassing epithelial separation.
 box(3.8,2.53,3.15,1.03,'Inducible transplants: Figure S7\nSame graft allocation\nRecovery after 14 days',color='#FAEBDD',size=9.5)
 branch=MplPath([(6.95,5.86),(7.15,5.86),(7.15,3.05),(6.95,3.05)],
                [MplPath.MOVETO,MplPath.LINETO,MplPath.LINETO,MplPath.LINETO])
 ax.add_patch(FancyArrowPatch(path=branch,arrowstyle='-|>',mutation_scale=12,
                             color='#555555',lw=.9,joinstyle='miter'))
 OUT.parent.mkdir(parents=True,exist_ok=True)
 for ext in ['pdf','svg','png']:fig.savefig(OUT.with_suffix('.'+ext),dpi=300)
 plt.close(fig)
 print('Updated Figure S1 from September experimental allocation')

if __name__=='__main__':main()
