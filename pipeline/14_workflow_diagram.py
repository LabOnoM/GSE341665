#!/usr/bin/env python3
"""Current Figure S5 workflow, aligned with the retained reproduction scripts.

Panel A records archived primary analysis, without claiming a new raw-read fit.
Panel B maps the current executable analyses.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT=Path(__file__).resolve().parent/'figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})

def canvas(height):
    fig,ax=plt.subplots(figsize=(170/25.4,height/25.4))
    ax.set_xlim(0,170);ax.set_ylim(height,0);ax.axis('off')
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    return fig,ax

def box(ax,x,y,width,height,text,accent=False):
    ax.add_patch(FancyBboxPatch((x-width/2,y-height/2),width,height,
        boxstyle='round,pad=0.8,rounding_size=1.5',linewidth=.9,
        edgecolor='#0072B2' if accent else '#707070',
        facecolor='#EAF3F9' if accent else '#F7F7F7'))
    ax.text(x,y,text,ha='center',va='center',fontsize=10,linespacing=1.15)

def arrow(ax,start,end):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=10,
        linewidth=.9,color='#555555',connectionstyle='arc3,rad=0'))

def save(fig,name):
    for suffix in ['pdf','svg','png']:
        fig.savefig(OUT/f'{name}.{suffix}',dpi=300,facecolor='white')
    plt.close(fig)

def primary():
    fig,ax=canvas(204)
    nodes=[
      'Provider-filtered paired-end FASTQ, 150 bp\n3 CON + 3 dcKO; 21.0–22.5 million pairs per sample',
      'FastQC → Trim Galore 0.6.10 / Cutadapt 4.9\nPhred 20; adapter detection; post-trimming FastQC',
      'STAR 2.7.10b two-pass alignment\nEnsembl 116 GRCm39; sjdbOverhang 149\n93.66–93.84% uniquely mapped',
      'StringTie merge → Salmon 1.10.0 (ISR) → tximport\nGene-level counts and TPM',
      'Retain ≥10 reads in ≥3 samples: 15,803 genes\nRUVSeq: one factor of unwanted variation (W₁)',
      'DESeq2 likelihood-ratio test\nFull: ~W₁ + genotype; reduced: ~W₁',
      'BH-adjusted P < 0.05 and |log₂ fold change| > 1\n563 DEGs: 209 upregulated, 354 downregulated',
      'Primary GO / KEGG: uncorrected P < 0.05\nQualitative corroboration; not multiplicity-corrected']
    ys=[15+24*i for i in range(len(nodes))]
    for i,(y,text) in enumerate(zip(ys,nodes)):
        box(ax,85,y,154,20,text,accent=i==6)
        if i:arrow(ax,(85,ys[i-1]+10.8),(85,y-10.8))
    ax.text(85,198,'Archived primary processing; the current repository does not refit this model.',
        ha='center',va='center',fontsize=8.5)
    save(fig,'FigureS5A_pipeline')

def reanalysis():
    fig,ax=canvas(193)
    box(ax,49,16,83,25,'All retained gene statistics\nand full-precision TPM\nfrom archived primary results')
    box(ax,133,16,65,25,'GSE146855\nP7 incisor reference\n(public GEO)')
    box(ax,49,51,83,25,'01  Extract source workbooks\nPreserve model statistics\nand all 15,803 gene rows')
    box(ax,133,51,65,25,'04 / 05  Download; Seurat\nQC and clustering\nSpecificity-based annotation')
    box(ax,49,87,83,25,'03  Four program scores\n16 identity; 20 sterol\n14 wound/stress; 18 interferon')
    box(ax,133,87,65,25,'06  NNLS deconvolution\nReference cell-type fractions\nand program specificity')
    box(ax,49,122,83,25,'17 / 16  Pinned PanelApp 269 v4.37\nMouse orthologs; directional\nFisher enrichment (GREEN genes)')
    box(ax,133,122,65,25,'pipeline/programs.py\nSingle program definition\nImported by 03 / 06 / 07',accent=True)
    box(ax,85,157,154,22,'07  Current Figures 4–7 and complete Data S3 source CSVs\nExact program tests; descriptive fractions; all displayed genes',accent=True)
    for x in [49,133]:
        arrow(ax,(x,29.5),(x,37.5));arrow(ax,(x,64.5),(x,73.5))
    arrow(ax,(49,100.5),(49,108.5));arrow(ax,(49,135.5),(66,145))
    # Route reference output around the program-definition note.
    ax.plot([166.5,168.5,168.5],[87,87,138],color='#555555',linewidth=.9)
    arrow(ax,(168.5,138),(150,145));arrow(ax,(133,135.5),(103,145))
    ax.text(85,179,'export_additional_data.py → Additional files and Tables S4–S5',ha='center',va='center',fontsize=10)
    arrow(ax,(85,169),(85,174))
    ax.text(85,189,'Current reproduction uses archived gene-level results, without a new count-model fit.',ha='center',va='center',fontsize=8.5)
    save(fig,'FigureS5B_reanalysis')

if __name__=='__main__':primary();reanalysis()
