#!/usr/bin/env python3
"""Rebuild nonredundant expression figures and their complete numerical source data.

This is a presentation/reporting revision of archived DESeq2 results, not a new
count-model analysis. Gene programs are imported unchanged from programs.py.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, itertools, json, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy import stats

ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,required=True)
ap.add_argument('--destination',type=Path,required=True)
ap.add_argument('--claim-source',type=Path,required=True,help='Captured manuscript paragraph evidence JSON from the data deposit')
ap.add_argument('--scratch',type=Path,required=True)
args=ap.parse_args()
ROOT=args.root; DEST=args.destination; OUT=ROOT/'scripts/JDR_reanalysis/output'
sys.path.insert(0,str(ROOT/'scripts/JDR_reanalysis'))
from programs import PROGRAMS
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor'))
import figstyle as fs
fs.setup(base_pt=8.7,style='neutral')
plt.rcParams.update({'font.sans-serif':['DejaVu Sans'],'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':True})
BLUE='#0072B2'; ORANGE='#D55E00'; GREY='#6C737A'
CMAP=LinearSegmentedColormap.from_list('expression',['#2166AC','#FFFFFF','#D55E00'])
NORM=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2)
SAMPLES=['TPM_CON_1','TPM_CON_2','TPM_CON_3','TPM_dcKO_1','TPM_dcKO_2','TPM_dcKO_3']
F=DEST/'Figures'; S=DEST/'Supplementary'
for d in [F,S,args.scratch]: d.mkdir(parents=True,exist_ok=True)
deg=pd.read_csv(OUT/'DEG_table.csv',skiprows=1,float_precision='round_trip',low_memory=False)
tpm=pd.read_csv(OUT/'TPM_matrix.csv',skiprows=1,float_precision='round_trip',low_memory=False)
assert deg.ensembl_gene_id.is_unique and tpm.ensembl_gene_id.is_unique
mod=pd.read_csv(OUT/'module_scores.csv',index_col=0,float_precision='round_trip')
dec=pd.read_csv(OUT/'deconv_proportions.csv',index_col=0,float_precision='round_trip')
specificity=pd.read_csv(OUT/'program_celltype_specificity.csv',index_col=0,float_precision='round_trip')
ai=pd.read_csv(OUT/'AI_gene_directions.csv',float_precision='round_trip')

DISPLAY={
 'ERK-associated':['Ereg','Areg','Hras','Dusp1','Dusp4','Dusp5','Dusp6','Spry2','Spry4','Fosl1','Egr1','Fos','Fosb','Dusp8','Dusp10','Dusp16'],
 'ER stress':['Atf3','Atf4','Ddit3','Hspa5','Trib3'],
 'p38 and TGF-beta context':['Gadd45a','Gadd45b','Gadd45g','Map2k6','Map2k3','Mapk14','Map3k7','Traf6','Tgfb2','Tgfb3','Tgfbr2','Tgfbr3','Grem1','Bmpr1a','Smad1','Smad5','Smad9','Smad3','Smad4'],
 'Additional epithelial markers':['Krt76','Gli1','Igfbp5','Krt15','Sox2','Lgr5','Msx2','Tuft1','Odam'],
 'Sex markers':['Xist','Ddx3y','Eif2s3y','Uty','Kdm5d'],
 'Secretory and maturation markers':['Amelx','Ambn','Klk4','Slc24a4']}
PANEL={'Identity':'Figure 5A','Wound':'Figure 5C','Sterol':'Figure 6A','IFN':'Figure 6B',
 'ERK-associated':'Figure 7A','ER stress':'Figure 7C','p38 and TGF-beta context':'Figure 7B',
 'Additional epithelial markers':'Figure 5E','Sex markers':'Figure 4D',
 'Secretory and maturation markers':'Figure 5F'}
all_groups={**PROGRAMS,**DISPLAY,'Human AI orthologs':ai.external_gene_name.tolist()}
PANEL['Human AI orthologs']='Figure 5B'
display_genes=set(sum(all_groups.values(),[]))
selected=deg[deg.external_gene_name.isin(display_genes)]
assert not selected.external_gene_name.duplicated().any(), 'Ambiguous gene symbols require explicit mapping'
bygene=selected.set_index('external_gene_name')
T=tpm.groupby('external_gene_name')[SAMPLES].mean()

def zscore_rows(x):
    x=np.log2(x+1)
    return x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1,ddof=1).replace(0,np.nan),axis=0)

program_tests=[]
for name,genes in PROGRAMS.items():
    assert set(genes)<=set(bygene.index)
    calc=zscore_rows(T.loc[genes]).mean(axis=0)
    assert np.allclose(calc.values,mod.loc[SAMPLES,name].values,rtol=1e-12,atol=1e-12)
    vals=calc.values; diff=vals[3:].mean()-vals[:3].mean()
    diffs=[abs(vals[list(c)].mean()-np.delete(vals,list(c)).mean()) for c in itertools.combinations(range(6),3)]
    p=sum(x>=abs(diff)-1e-12 for x in diffs)/len(diffs)
    program_tests.append({'program':name,'n_genes':len(genes),'mean_CON':vals[:3].mean(),'mean_dcKO':vals[3:].mean(),'delta_score':diff,'exact_permutation_P':p,'n_libraries_per_group':3,'n_assignments':len(diffs)})
assert sum((deg.padj<.05)&(deg.log2FoldChange.abs()>1))==563
assert sum(deg.padj==0)==5
assert all(x['exact_permutation_P']==.1 for x in program_tests)

def pformat(p):
    if pd.isna(p): return 'NA'
    if p==0: return '0*'
    if p<.001: return f'{p:.1e}'.replace('e-0','e-')
    return f'{p:.3f}'

def newfig(w,h):
    fig=plt.figure(figsize=(w/25.4,h/25.4),facecolor='white');fig._mm=(w,h);return fig

def axes(fig,x,y,w,h):
    fw,fh=fig._mm
    return fig.add_axes([x/fw,1-(y+h)/fh,w/fw,h/fh])

def title(fig,x,y,s):
    w,h=fig._mm;fig.text(x/w,1-y/h,s,ha='left',va='top',weight='bold',fontsize=10.2)

def heat_legend(fig,y):
    w,h=fig._mm
    ax=axes(fig,12,y,23,2.2)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=ax,orientation='horizontal',ticks=[-2,0,2])
    cb.ax.tick_params(labelsize=8,length=2,pad=1);cb.outline.set_visible(False)
    fig.text(39/w,1-(y+1)/h,'Row z of log₂(TPM + 1)',va='center',fontsize=8.6)
    fig.text(100/w,1-(y+1)/h,'log₂FC: dcKO / CON',va='center',fontsize=8.6)
    fig.text(100/w,1-(y+6)/h,'adj. P: DESeq2 BH-adjusted',va='center',fontsize=8.6)

def gene_panel(fig,rect,genes,label,heat=True,lim=(-4.5,4.5),rowfontsize=8.5):
    x,y,w,h=rect;title(fig,x,y,label)
    ax=axes(fig,x,y+9,w,h-9)
    n=len(genes);ax.set_xlim(0,100);ax.set_ylim(n+.55,-2.1);ax.axis('off')
    hx=(23,43) if heat else None
    fx=(49,67) if heat else (27,65)
    qx=91;nx=75
    def xx(v):return fx[0]+(v-lim[0])/(lim[1]-lim[0])*(fx[1]-fx[0])
    if heat:
        z=zscore_rows(T.reindex(genes));zval=z.to_numpy()
        ax.imshow(np.ma.masked_invalid(zval),extent=[hx[0],hx[1],n-.5,-.5],aspect='auto',cmap=CMAP,norm=NORM,interpolation='nearest')
        ax.axvline((hx[0]+hx[1])/2,ymin=(1.05)/(n+2.65),ymax=(n+1.05)/(n+2.65),color='white',lw=1.6)
        for center,txt,col in [(28,'CON',BLUE),(38,'dcKO',ORANGE)]:
            ax.text(center,-1.48,txt,ha='center',va='center',color=col,fontsize=8)
        for j in range(6): ax.text(hx[0]+(j+.5)*(hx[1]-hx[0])/6,-.81,str(j%3+1),ha='center',va='center',fontsize=7.8)
    for loc,txt in [(0,'Gene'),((fx[0]+fx[1])/2,'log₂FC'),(nx,'Value'),(qx,'adj. P')]:
        ax.text(loc,-1.2,txt,ha='left' if loc==0 else 'center',va='center',fontsize=8.3)
    ax.plot([xx(0)]*2,[-.5,n-.5],c='#B6BCC2',lw=.65,zorder=0)
    for i,g in enumerate(genes):
        ax.text(0,i,g,va='center',ha='left',style='italic',fontsize=rowfontsize)
        if g not in bygene.index:
            for loc in [nx,qx]:ax.text(loc,i,'NA',ha='center',va='center',fontsize=rowfontsize,color=GREY)
            if heat:ax.text((hx[0]+hx[1])/2,i,'not retained',ha='center',va='center',fontsize=6.6,color=GREY)
            continue
        r=bygene.loc[g];lfc=float(r.log2FoldChange);p=float(r.padj)
        col=ORANGE if lfc>0 else BLUE
        if not np.isfinite(p) or p>=.05:col=GREY
        ax.plot([xx(0),xx(lfc)],[i,i],color=col,lw=1.0)
        ax.scatter(xx(lfc),i,s=12,edgecolor=col,facecolor=col if p<.05 else 'white',linewidth=.7,zorder=3)
        fc_label='0.00' if abs(lfc)<.005 else f'{lfc:+.2f}'.replace('-','−')
        ax.text(nx,i,fc_label,ha='center',va='center',fontsize=rowfontsize)
        ax.text(qx,i,pformat(p),ha='center',va='center',fontsize=rowfontsize)
    ax.plot(fx,[n-.20,n-.20],color='#9DA5AD',lw=.55)
    ticks=([lim[0],0,lim[1]] if lim!=(-4.5,4.5) else [-4,0,4])
    for v in ticks:
        ax.text(xx(v),n+.22,f'{v:g}'.replace('-','−'),ha='center',va='center',fontsize=7.7)
    return ax

figure_meta=[]
def save(fig,path):
    w,h=fig._mm
    for ext in ['png','pdf','svg']:
        fig.savefig(path.with_suffix('.'+ext),dpi=300,facecolor='white')
    minfont=min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text() and t.get_visible())
    figure_meta.append({'file':path.name,'width_mm':w,'height_mm':h,'minimum_font_pt':minfont,'png_dpi':300,
                        'formats':['png','pdf','svg'],'svg_text':'editable text'})
    plt.close(fig);print('WROTE',path.name,flush=True)

# Figure 4: retain a single program overview with all six libraries.
fig=newfig(180,230)
title(fig,3,3,'A  Transcriptome separation')
ax=axes(fig,15,17,61,75)
lg=np.log2(T+1);top=lg.var(axis=1).sort_values(ascending=False).head(2000).index
M=lg.loc[top].T.to_numpy();M-=M.mean(0);U,s,V=np.linalg.svd(M,full_matrices=False)
pc=U[:,:2]*s[:2];ev=(s*s/(s*s).sum())[:2]*100
for sl,c,m,g in [(slice(0,3),BLUE,'o','CON'),(slice(3,6),ORANGE,'s','dcKO')]:
    ax.scatter(pc[sl,0],pc[sl,1],s=31,c=c,marker=m,label=g,edgecolors='white',lw=.5)
ax.set(xlabel=f'PC1 ({ev[0]:.1f}%)',ylabel=f'PC2 ({ev[1]:.1f}%)');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,fontsize=8.5)
title(fig,91,3,'B  Differential expression')
ax=axes(fig,103,17,72,75)
d=deg[deg.padj.notna()].copy();sig=(d.padj<.05)&(d.log2FoldChange.abs()>1)
fin=d.padj>0;nlp=-np.log10(d.loc[fin,'padj']);yy=np.sqrt(nlp)
xx=d.loc[fin,'log2FoldChange'];colors=np.where(~sig[fin],'#B6BCC2',np.where(xx>0,ORANGE,BLUE))
ax.scatter(np.clip(xx,-6.2,6.2),yy,s=np.where(sig[fin],4.5,1.5),c=colors,alpha=.65,lw=0,rasterized=False)
ymax=float(np.sqrt(nlp.max()));zero_y=ymax+1.05
zeros=d[d.padj==0];ax.scatter(np.clip(zeros.log2FoldChange,-6.2,6.2),np.repeat(zero_y,len(zeros)),c=np.where(zeros.log2FoldChange>0,ORANGE,BLUE),marker='^',s=18)
ax.axhline(zero_y,color='#999999',ls=':',lw=.7,zorder=0)
ax.text(-6.3,zero_y+.33,f'0* ({len(zeros)} genes)',fontsize=7.5,ha='left')
ax.axhline(np.sqrt(-np.log10(.05)),color='#666666',ls='--',lw=.6)
for val in [-1,1]:ax.axvline(val,color='#777777',ls='--',lw=.6)
ticks=[0,5,20,50,100,200];ax.set_yticks(np.sqrt(ticks),[str(x) for x in ticks])
ax.set(xlim=(-6.4,6.4),ylim=(0,zero_y+1),xlabel='log₂ fold change (dcKO / CON)',ylabel='−log₁₀ adjusted P (square-root scale)')
ax.text(.02,.90,f'{sum(sig & (d.log2FoldChange<0))} down',transform=ax.transAxes,color=BLUE,weight='bold',fontsize=8)
ax.text(.98,.90,f'{sum(sig & (d.log2FoldChange>0))} up',transform=ax.transAxes,color=ORANGE,ha='right',weight='bold',fontsize=8)
from volcano_panel import LabelPlacer
lp=LabelPlacer(fig);lp.add_point_cloud(ax,np.clip(xx,-6.2,6.2),yy)
items=[]
for gene in ['Amelx','Ambn','Krt76','Ereg','Ldlr','Isg15','Mx1','Hmgcr']:
    r=bygene.loc[gene]; y=np.sqrt(-np.log10(r.padj));lp.add_anchor(gene,ax,r.log2FoldChange,y)
    items.append({'text':gene,'ax':ax,'x':float(r.log2FoldChange),'y':float(y)})
# LabelPlacer's public input contract is defined in volcano_panel.py.
try: lp.place(items,fontsize=7.6)
except (TypeError,ValueError):
    # No fallback positions or numbers: the plot remains valid without selected labels.
    raise
title(fig,3,109,'C  Program scores across individual libraries')
ax=axes(fig,17,122,156,29)
order=['Identity','IFN','Sterol','Wound']
for i,name in enumerate(order):
    for sl,col,dx,marker in [(slice(0,3),BLUE,-.17,'o'),(slice(3,6),ORANGE,.17,'s')]:
        v=mod.loc[SAMPLES,name].to_numpy()[sl]
        ax.scatter(i+dx+np.array([-.045,0,.045]),v,s=19,c=col,marker=marker,edgecolors='white',lw=.3,zorder=3)
        ax.plot([i+dx-.11,i+dx+.11],[v.mean()]*2,color=col,lw=1.4)
    p=next(x['exact_permutation_P'] for x in program_tests if x['program']==name)
    ax.text(i,1.30,f'P = {p:.2f}',ha='center',fontsize=8.4)
ax.axhline(0,color='#999999',lw=.6);ax.set_ylim(-1.15,1.6);ax.set_xlim(-.5,3.5)
ax.set_xticks(range(4),['Ameloblast identity','Interferon','Sterol biosynthesis','Wound / stress']);ax.set_ylabel('Program z-score');ax.set_yticks([-1,0,1]);ax.tick_params(labelsize=8.4)
fig.text(.5,1-161/230,'Points: libraries; short lines: means; exact two-sided permutation tests.',ha='center',fontsize=8.3)
gene_panel(fig,(3,171,174,42),DISPLAY['Sex markers'],'D  Sex-associated expression across the same libraries')
heat_legend(fig,218)
save(fig,F/'Figure4_transcriptome')

# Figure 5: complete identity and wound programs, AI orthologs, and composition.
fig=newfig(180,255);heat_legend(fig,6)
gene_panel(fig,(3,22,84,82),PROGRAMS['Identity'],'A  Ameloblast identity')
gene_panel(fig,(3,111,84,72),PROGRAMS['Wound'],'C  Wound / stress')
ag=ai.sort_values('log2FoldChange').external_gene_name.tolist()
gene_panel(fig,(93,22,84,125),ag,'B  Human AI gene orthologs',heat=False,lim=(-1.5,2.1),rowfontsize=8.3)
sens=pd.read_csv(OUT/'AI_enrichment_sensitivity.csv',float_precision='round_trip')
r=sens[(sens.resolution=='CORRECTED')&(sens.lfc_cutoff==0)].iloc[0]
fig.text(94/180,1-151/255,f'Down-enrichment: OR {r.odds_ratio:.2f}; P {r.p_one_sided:.1e}',fontsize=8.0)
title(fig,93,158,'D  Estimated tissue composition')
ax=axes(fig,105,169,68,27)
keep=[c for c in dec if dec[c].mean()>.005];pvals=[]
for i,name in enumerate(keep):
    a=dec.loc[SAMPLES[:3],name].values*100;b=dec.loc[SAMPLES[3:],name].values*100
    p=stats.ttest_ind(a,b,equal_var=True).pvalue;pvals.append({'cell_type':name,'mean_CON_pct':a.mean(),'mean_dcKO_pct':b.mean(),'P_unpaired_t':p})
    for v,c,dx,m in [(a,BLUE,-.17,'o'),(b,ORANGE,.17,'s')]:
        ax.scatter(i+dx+np.array([-.04,0,.04]),v,c=c,s=17,marker=m,zorder=3)
        ax.plot([i+dx-.1,i+dx+.1],[v.mean()]*2,color=c,lw=1.2)
    ax.text(i,max(a.max(),b.max())+6,f'P={p:.3f}',ha='center',fontsize=7.8)
ax.set(ylim=(0,78),ylabel='Estimated fraction (%)');ax.set_xticks(range(3),['AB\nmaturation','AB\nsecretory','Mesenchyme']);ax.tick_params(labelsize=8)
gene_panel(fig,(3,189,84,52),DISPLAY['Additional epithelial markers'],'E  Additional epithelial genes',lim=(-5,5))
title(fig,93,209,'F  Secretory / maturation markers')
ax=axes(fig,107,220,66,20)
for i,g in enumerate(DISPLAY['Secretory and maturation markers']):
    for sl,c,dx,m in [(slice(0,3),BLUE,-.15,'o'),(slice(3,6),ORANGE,.15,'s')]:
        v=T.loc[g].values[sl];ax.scatter(i+dx+np.array([-.045,0,.045]),np.log10(v+1),s=16,c=c,marker=m)
ax.set_xticks(range(4),DISPLAY['Secretory and maturation markers']);ax.set_yticks([0,2,4],['0','99','9,999']);ax.set_ylabel('TPM (log scale)');ax.tick_params(labelsize=8)
for l in ax.get_xticklabels():l.set_fontstyle('italic')
fig.text(.02,.015,'NA: not retained for DESeq2. 0*: numerical underflow. Full source values: Data S3.',fontsize=8.3)
save(fig,F/'Figure5_identity_switch')

# Figure 6: gene-level sterol / IFN evidence replaces both redundant score views.
fig=newfig(180,205);heat_legend(fig,6)
gene_panel(fig,(3,22,84,105),PROGRAMS['Sterol'],'A  Sterol biosynthesis')
gene_panel(fig,(93,22,84,105),PROGRAMS['IFN'],'B  Interferon response')
title(fig,3,137,'C  Gene-program specificity in the single-cell reference')
ax=axes(fig,25,150,119,27)
sp=specificity.loc[['Identity','Sterol','Wound','IFN']]
im=ax.imshow(sp.values,aspect='auto',cmap=CMAP,norm=TwoSlopeNorm(vmin=-1.6,vcenter=0,vmax=1.6))
ax.set_yticks(range(4),['Identity','Sterol','Wound','IFN']);ax.set_xticks(range(7),['AB mat.','AB sec.','Endo.','IEE/SI','Erythr.','Immune','Mes.'],rotation=60,ha='right');ax.tick_params(length=0,labelsize=7.9)
for i in range(4):
    for j in range(7): ax.text(j,i,f'{sp.iloc[i,j]:.1f}'.replace('-0.0','0.0'),ha='center',va='center',fontsize=7.5,color='white' if abs(sp.iloc[i,j])>1 else '#182029')
for s in ax.spines.values():s.set_visible(False)
cax=axes(fig,154,151,2,25);cb=fig.colorbar(im,cax=cax,orientation='vertical',ticks=[-1.6,0,1.6]);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=2,pad=1);cb.set_label('Reference z',fontsize=8.5)
fig.text(.02,.015,'0*: numerical underflow in source. Expression programs do not measure sterol flux.',fontsize=8.4)
save(fig,F/'Figure6_metabolic_convergence')

# Figure 7: all additional signaling genes belong in the main figures, per author.
fig=newfig(180,199);heat_legend(fig,6)
gene_panel(fig,(3,22,84,112),DISPLAY['ERK-associated'],'A  ERK-associated / MAPK feedback')
gene_panel(fig,(93,22,84,112),DISPLAY['p38 and TGF-beta context'],'B  p38 / TGF-β / BMP pathway genes')
gene_panel(fig,(3,145,174,42),DISPLAY['ER stress'],'C  ER stress')
fig.text(.02,.013,'All listed genes are shown, including nonsignificant results. Full numerical records: Data S3.',fontsize=8.4)
save(fig,F/'Figure7_signaling_gene_evidence')

# Full expression results with gene, program, manuscript and panel mapping.
claim_source=json.loads(args.claim_source.read_text())
body=[(p['index'],p['text']) for p in claim_source['paragraphs']]
source_manuscript_sha256=claim_source['source_manuscript_sha256']
colbase=['ensembl_gene_id','external_gene_name','baseMean','FoldChange','log2FoldChange','lfcSE','stat','pvalue','padj']
full=deg[colbase].merge(tpm[['ensembl_gene_id']+SAMPLES],on='ensembl_gene_id',validate='one_to_one')
full['mean_TPM_CON']=full[SAMPLES[:3]].mean(axis=1);full['mean_TPM_dcKO']=full[SAMPLES[3:]].mean(axis=1);full['mean_TPM_all']=full[SAMPLES].mean(axis=1)
full['strict_DEG']=(full.padj<.05)&(full.log2FoldChange.abs()>1)
full['padj_status']=np.where(full.padj.isna(),'not available',np.where(full.padj==0,'source numeric underflow','finite'))
full['program_membership']=full.external_gene_name.map(lambda g:'; '.join(k for k,v in PROGRAMS.items() if g in v))
full['figure_panels']=full.external_gene_name.map(lambda g:'; '.join(dict.fromkeys(PANEL[k] for k,v in all_groups.items() if g in v)))
full['main_text_paragraphs']=full.external_gene_name.map(lambda g:'; '.join(str(i) for i,txt in body if isinstance(g,str) and re.search(r'(?<![A-Za-z0-9])'+re.escape(g)+r'(?![A-Za-z0-9])',txt)))
full=full.rename(columns={'external_gene_name':'gene','FoldChange':'source_FoldChange'})
full.to_csv(S/'Supplementary_Data_S3_all_gene_statistics.csv',index=False,float_format='%.17g')
tables={}
displayrows=[]
for k,gs in all_groups.items():
    for g in gs:
        rows=full[full.gene==g]
        if len(rows): row=rows.iloc[0].to_dict()
        else: row={'gene':g,'padj_status':'not retained; no DESeq2 result'}
        row={'display_group':k,'display_panel':PANEL[k],**row};displayrows.append(row)
display=pd.DataFrame(displayrows)
display.to_csv(S/'Supplementary_Data_S3_displayed_genes.csv',index=False,float_format='%.17g')
pd.DataFrame(program_tests).to_csv(S/'Supplementary_Data_S3_program_tests.csv',index=False,float_format='%.17g')
pd.DataFrame(pvals).to_csv(S/'Supplementary_Data_S3_deconvolution_tests.csv',index=False,float_format='%.17g')

# Auxiliary count evidence is kept separate from Salmon/tximport model inputs.
counts=[]
for gene,gid in [('Sox2','ENSMUSG00000074637'),('Lgr5','ENSMUSG00000020140')]:
    for sample in ['CON1','CON2','CON3','dcKO1','dcKO2','dcKO3']:
        path=ROOT/'00.RawData/Bulk-RNAseq'/f'{sample}.STAR.counts'
        for line in path.open():
            if line.startswith(gid+'\t'):
                counts.append({'gene':gene,'ensembl_gene_id':gid,'sample':sample,'alignment_feature_count':int(line.rstrip().split('\t')[-1]),'source_file':path.name,'interpretation':'Auxiliary alignment counts; not DESeq2 input; no differential-expression test available'});break
assert len(counts)==12
pd.DataFrame(counts).to_csv(S/'Supplementary_Data_S3_nonretained_gene_counts.csv',index=False)

# Source claims table lists every original numeric-expression paragraph, plus all
# explicitly named genes, so unsupported/omitted citations can be reviewed directly.
claimrows=[]
for i,txt in body:
    if ('fold' in txt or 'TPM' in txt or 'transcripts per million' in txt or 'adjusted P' in txt) and i>=32:
        hits=full[full.main_text_paragraphs.map(lambda x:str(i) in x.split('; '))]
        for _,r in hits.iterrows():
            claimrows.append({'source_paragraph':i,'gene':r.gene,'log2FoldChange':r.log2FoldChange,'fold_dcKO_over_CON':2**r.log2FoldChange,'adjusted_P':r.padj,'padj_status':r.padj_status,'mean_TPM_CON':r.mean_TPM_CON,'mean_TPM_dcKO':r.mean_TPM_dcKO,'mean_TPM_all':r.mean_TPM_all,'figure_panels':r.figure_panels,'supporting_table':'Supplementary Data S3'})
claims=pd.DataFrame(claimrows);claims.to_csv(S/'Supplementary_Data_S3_text_claim_index.csv',index=False,float_format='%.17g')

def matrix(df):
    return [df.columns.tolist()]+df.astype(object).where(pd.notna(df),None).values.tolist()
compactcols=['display_group','gene','log2FoldChange','padj','padj_status','mean_TPM_CON','mean_TPM_dcKO']+SAMPLES+['display_panel','ensembl_gene_id']
readme=[['Field','Definition'],
 ['Comparison','dcKO versus CON; 3 pooled biological libraries per genotype.'],
 ['Gene-level inference','Archived DESeq2 likelihood-ratio test: full ~W1 + genotype; reduced ~W1. P values adjusted by Benjamini-Hochberg.'],
 ['Expression','TPM is descriptive, unadjusted expression. Heatmaps are row z-scores of log2(TPM + 1), sample SD (ddof=1).'],
 ['Effect size','log2FoldChange comes from the fitted model; not recomputed as a ratio of mean TPM. Linear fold = 2^log2FoldChange.'],
 ['No selection','All members of the four predefined programs are displayed, including effects below the strict DEG cutoff.'],
 ['Zero adjusted P','Source values of zero denote numerical underflow; figures show 0*. They are not literal zero probabilities.'],
 ['Missing adjusted P','A missing test result is retained as missing, not converted to zero.'],
 ['Program tests','Two-sided exact label permutation over all 20 allocations of 6 libraries; smallest attainable P is 0.10.'],
 ['Source data','DEG_table.csv and TPM_matrix.csv; original workbooks STAR_DEGs_ResOut_01_dcKOvsCON.xlsx and 00.DEGs_Out_TPM.xlsx.'],
 ['Program membership','programs.py is unchanged. ERK, stress and Discussion groups are display-only, not new scored programs.'],
 ['Nonretained genes','Sox2 and Lgr5 have low nonzero auxiliary alignment feature counts but no row in the retained DESeq2 table. These counts are not the Salmon/tximport model inputs.'],
 ['Animal sex','All experimental donors and all six RNA-seq libraries were male (author confirmation, 3 October 2026). Sex-marker expression is retained as measured; it does not establish mixed-sex pooling.'],
 ['Panel numbering','4C program overview; 4D sex markers; 5A identity; 5B human AI orthologs; 5C wound; 5D composition; 5E additional epithelial genes; 5F secretory/maturation TPM; 6A sterol; 6B IFN; 6C reference specificity; 7A ERK/MAPK; 7B p38/TGF-beta/BMP; 7C ER stress.'],
 ['Manuscript','Source paragraph numbers in Text claims are zero-based positions in the captured author DOCX, not printed paragraph or section numbers.'],
 ['Source manuscript SHA256',source_manuscript_sha256]]
payload={'sheets':[
 {'name':'Displayed genes','data':matrix(display[compactcols]),'numeric':['log2FoldChange','padj','mean_TPM_CON','mean_TPM_dcKO']+SAMPLES},
 {'name':'Text claims','data':matrix(claims),'numeric':['source_paragraph','log2FoldChange','fold_dcKO_over_CON','adjusted_P','mean_TPM_CON','mean_TPM_dcKO','mean_TPM_all']},
 {'name':'All genes','data':matrix(full),'numeric':[c for c in full if pd.api.types.is_numeric_dtype(full[c])]},
 {'name':'Program tests','data':matrix(pd.DataFrame(program_tests)),'numeric':[c for c in program_tests[0] if c!='program']},
 {'name':'Deconvolution tests','data':matrix(pd.DataFrame(pvals)),'numeric':[c for c in pvals[0] if c!='cell_type']},
 {'name':'Nonretained genes','data':matrix(pd.DataFrame(counts)),'numeric':['alignment_feature_count']},
 {'name':'Read me','data':readme,'numeric':[]}]}
(args.scratch/'workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False))
sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'DEG_table.csv',OUT/'TPM_matrix.csv',OUT/'module_scores.csv',OUT/'deconv_proportions.csv',OUT/'program_celltype_specificity.csv',ROOT/'scripts/JDR_reanalysis/programs.py']}
prov={'source_manuscript_sha256':source_manuscript_sha256,'sources':sources,'program_tests':program_tests,'figures':figure_meta,'gene_counts':{'retained':len(full),'finite_adjusted_P':int(deg.padj.notna().sum()),'source_zero_adjusted_P':int(sum(deg.padj==0)),'display_rows':len(display),'unique_display_genes':len(set(display.gene)),'claim_index_rows':len(claims)},'analysis':'Archived DESeq2 values retained; program scores and permutation tests recomputed and matched; no new count-model fit.'}
(DEST/'Figure_Refinement_Provenance.json').write_text(json.dumps(prov,indent=2))
print(json.dumps(prov['gene_counts']),flush=True)
