#!/usr/bin/env python3
"""Render Supplementary Figure S5 as one composite with panels A and B."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent / 'figures'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
                     'pdf.fonttype':42, 'svg.fonttype':'none'})


def render():
    fig, ax = plt.subplots(figsize=(180/25.4, 184/25.4))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(0,180), ylim=(184,0)); ax.axis('off')
    text_bounds = []

    def box(x,y,w,h,text,accent=False):
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,
            boxstyle='round,pad=0.5,rounding_size=1', linewidth=.75,
            edgecolor='#0072B2' if accent else '#707070',
            facecolor='#EAF3F9' if accent else '#F7F7F7'))
        obj=ax.text(x,y,text,ha='center',va='center',fontsize=9,linespacing=1.13)
        text_bounds.append((obj,(x-w/2,y-h/2,x+w/2,y+h/2)))

    def arrow(start,end):
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',
            mutation_scale=8,linewidth=.8,color='#555555'))

    def bus(xs,y,y0,y1):
        for x in xs:
            ax.plot([x,x],[y0,y],color='#555555',lw=.8)
        ax.plot([min(xs),max(xs)],[y,y],color='#555555',lw=.8)
        arrow((90,y),(90,y1))

    ax.text(4,5,'A',fontweight='bold',fontsize=12,va='center')
    ax.text(12,5,'RNA sequencing and gene-level analysis',fontweight='bold',fontsize=10,va='center')
    box(31,22,52,23,'Dental epithelium\n3 CON + 3 dcKO libraries\n150-bp paired-end reads')
    box(90,22,52,23,'FastQC and Trim Galore\nAdapter and quality\ntrimming')
    box(149,22,52,23,'STAR two-pass alignment\nGRCm39\nEnsembl release 116')
    arrow((58,22),(63,22)); arrow((117,22),(122,22))
    arrow((149,34),(149,39))
    box(149,52,52,23,'StringTie → Salmon\ntximport\nGene counts and TPM')
    box(90,52,52,23,'Low-count filtering\n15,803 genes retained\nRUVSeq: W₁ (k = 1)')
    box(31,52,52,23,'DESeq2 likelihood-ratio test\nFull: W₁ + genotype\nReduced: W₁')
    arrow((122,52),(117,52)); arrow((63,52),(58,52))
    arrow((31,64),(31,69))
    box(31,80,52,19,'563 differentially\nexpressed genes\n209 up; 354 down',True)
    box(100,80,65,19,'GO and KEGG\nenrichment analysis')
    arrow((58,80),(66.5,80))

    ax.text(4,99,'B',fontweight='bold',fontsize=12,va='center')
    ax.text(12,99,'Gene-program and reference-based analyses',fontweight='bold',fontsize=10,va='center')
    box(90,114,171,20,'Bulk expression data\nGene-level statistics, counts, TPM, effective lengths and W₁')
    box(31,146,52,29,'Four curated gene programs\nTMM → voom → ROAST\nW₁ + genotype\nDescriptive z scores',True)
    box(90,146,52,29,'AI-gene enrichment\nPanelApp 269, v4.37\nMGI mouse orthologs\nOne-sided Fisher test',True)
    box(149,146,52,29,'Reference-based analyses\nGSE146855 (P7 incisor)\nSeurat annotation\nNNLS; gene-set specificity',True)
    for x in (31,90,149):
        arrow((x,124.5),(x,131))
    bus((31,90,149),166,162,170)
    box(90,177,171,12,'Figures 4–7 and Additional-file tables',True)
    # Fail early if a label escapes its allotted box at the export size.
    fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    for obj,(x0,y0,x1,y1) in text_bounds:
        bb=obj.get_window_extent(renderer)
        corners=ax.transData.transform([(x0,y1),(x1,y0)])
        assert bb.x0 >= corners[0,0] and bb.x1 <= corners[1,0], obj.get_text()
        assert bb.y0 >= corners[0,1] and bb.y1 <= corners[1,1], obj.get_text()
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','svg','png'):
        fig.savefig(OUT/f'FigureS5_workflow.{ext}',dpi=300,facecolor='white')
    plt.close(fig)

if __name__ == '__main__':
    render()
