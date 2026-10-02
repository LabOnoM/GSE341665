#!/usr/bin/env python3
"""Draw the author-confirmed September experimental allocation (Figure S6).

Authority: September Methods and WeChat notes, ingested as
.wiki/sources/september_2026_experimental_update.md. The superseded July generator
is preserved in _revision_history/2026-09-24/before_update.zip.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'scripts/JDR_reanalysis/figures/FigureS6_breeding_scheme'

# Explicit design metadata audited independently of drawing placement.
DESIGN={
 'constitutive':{'driver':'K14-cre','age':'E16.5','4OHT':False,'figures':['1']},
 'inducible':{'driver':'K14-creER','age':'P0','4OHT':True,'figures':['2','3','4','5','6','S1','S2','S7']},
 'transplant':{'recovery_days':14,'recipient_age_weeks':8,'paired_recipients':3},
}

def main():
 assert DESIGN['constitutive']['figures']==['1'] and not DESIGN['constitutive']['4OHT']
 assert '3' in DESIGN['inducible']['figures'] and DESIGN['inducible']['4OHT']
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
 fig,ax=plt.subplots(figsize=(7.2,7.9));fig.subplots_adjust(0,0,1,1);ax.set_xlim(0,7.2);ax.set_ylim(0,7.9);ax.axis('off')
 ax.text(.25,7.58,'Experimental models and sample collection',fontweight='bold',fontsize=14)
 def box(x,y,w,h,txt,color='#EEF2F5',size=10):
  ax.add_patch(Rectangle((x,y),w,h,facecolor=color,edgecolor='#708090',lw=.8))
  ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=size,linespacing=1.3)
 def arrow(x,y0,y1):ax.add_patch(FancyArrowPatch((x,y0),(x,y1),arrowstyle='-|>',mutation_scale=12,color='#555555',lw=.9))
 box(.25,6.57,3.15,.68,'Constitutive K14-cre × Smad4 fl/fl\nCre-negative littermate controls','#DCE6F4')
 box(3.8,6.57,3.15,.68,'Inducible K14-creER × Smad4 fl/fl\nCre-negative littermate controls','#FAEBDD')
 arrow(1.825,6.57,6.31);arrow(5.375,6.57,6.31)
 box(.25,5.48,3.15,.83,'E16.5 first-molar germs\n24-h culture; no 4-OHT\nConstitutive model: Figure 1')
 box(3.8,5.48,3.15,.83,'P0 first-molar germs\n90 µmol/L 4-OHT for 24 h\nBoth genotypes treated identically')
 arrow(1.825,5.48,5.20);arrow(5.375,5.48,5.20)
 box(.25,3.88,3.15,1.32,'Subrenal-capsule transplantation\n8-week-old male recipients\nOne mutant + one control per kidney\nRecovery after 14 days\nn = 3 paired recipients')
 box(3.8,3.88,3.15,1.32,'Mechanical epithelial separation\nRT-qPCR / immunoblot: Figures 2–3\nRNA sequencing: Figures 4–7\nInhibitor-response assays: Figure 3\nTissue-region image: Figure S7',size=9.5)
 arrow(1.825,3.88,3.62)
 box(.25,2.90,3.15,.72,'Figure 1\nRadiography, H&E, AMBN staining','#DCE6F4')
 # The transplant and dose-series notes are separate protocol branches, not
 # downstream steps of molecular tissue collection.
 box(3.8,2.53,3.15,1.03,'Inducible transplants: Figure S1\nSame paired-recipient design\nRecovery after 14 days',color='#FAEBDD',size=9.5)
 ax.plot([6.98,6.98,6.98],[5.86,3.52,3.05],color='#777777',lw=.8)
 ax.annotate('',xy=(6.94,3.05),xytext=(6.98,3.05),arrowprops={'arrowstyle':'->','lw':.8})
 ax.text(.3,2.26,'Inhibitors (24 h)',fontsize=10,fontweight='bold')
 ax.text(.3,2.02,'SB203580 10 µmol/L; dorsomorphin 10 µmol/L;\nSB525334 1 µmol/L; U0126 10 µmol/L; matched DMSO vehicle.',fontsize=9.5,va='top')
 ax.text(.3,1.36,'Scope and record checks',fontsize=10,fontweight='bold')
 ax.text(.3,1.16,'The five inhibitor groups belong to the inducible / 4-OHT-treated arm; individual\nCre assignments are not encoded in the supplied tables. Figure S3 is an external\nsingle-cell reference. Figure S5 examines the sequenced allele. Stock 016230,\nrecorded for constitutive K14-cre, conflicts with the JAX catalog; stock identity\nremains unresolved. Figure S2 shows the separate dose series (50/90/120 µmol/L),\nwhich preceded selection of 90 µmol/L.',fontsize=9,va='top',linespacing=1.35)
 for ext in ['pdf','png']:fig.savefig(OUT.with_suffix('.'+ext),dpi=300)
 plt.close(fig)
 print('Updated Figure S6 from September experimental allocation')

if __name__=='__main__':main()
