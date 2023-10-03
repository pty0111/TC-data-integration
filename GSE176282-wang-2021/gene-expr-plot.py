import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt
import palantir as pl
from matplotlib.backends.backend_pdf import PdfPages
import pyreadr as pr
import re
import os
matplotlib.rcParams['font.family'] = ['serif']

os.chdir('/Users/pty0111/CCR7_DC/GSE176282-wang-2021')

def plot_gene_expr(expr, vis, dim1, dim2, genes, file, n_cols, s=3, 
                   newGeneName=None):
    cmap = matplotlib.cm.Spectral_r
    fig = pl.plot.FigureGrid(len(genes), n_cols)
    for g, ax in zip(genes, fig):
        c = expr.loc[vis.index, g.upper()]
        ax.scatter(vis.loc[:, dim1], vis.loc[:, dim2], s=s, c=c, cmap=cmap)
        ax.set_axis_off()
        if newGeneName and g in newGeneName:
            figTitle = newGeneName[g]
        else:
            figTitle = g
        ax.set_title(figTitle)
        normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cb = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
        #cb.ax.set_title('Expression')
        cb.set_label('Expression')
    n_rows = math.ceil(len(genes) / n_cols)
    fig.figure.set_size_inches(7 * n_cols, 5 * n_rows)
    file.savefig(fig.figure, bbox_inches='tight')
    return 0


md = pd.read_csv('results/meta-data.csv', header=0, index_col=0)
# umap = pd.read_csv('results/UMAP.csv', header=0, index_col=0)
umap = pd.read_csv('results/UMAP-subset.csv', header=0, index_col=0)
imp_df = pd.read_csv('results/imputed-expr.csv', header=0, index_col=0)
# imp_df = pr.read_r('results/imputed-expr.rds')[None]
imp_df = imp_df.transpose()
# imp_df.index = imp_df.index.astype('int64')
unimp_df = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0)
# unimp_df = pr.read_r('results/unimputed-expr.rds')[None].transpose()
unimp_df = unimp_df.transpose()
unimp_df.index = unimp_df.index.str.replace(".", "-")

genes = imp_df.columns
os.mkdir("plots/gene-expr")
############################################################################
# main genes, all cells
# select cells
ci = md.index
geneList = [['Rorc', 'Rora', 'Cxcr6', 'Aire', "H2-Ab1"]]
ncols = 2
for gene_list in geneList:
    for gene in gene_list:
        if gene.upper() not in genes:
            print(gene)
# imputed
for i, title in zip(range(len(geneList)), ['main_markers']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-gene-expr-imputed-' + title + '.pdf')
    plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()

# unimputed
for i, title in zip(range(len(geneList)), ['main_markers']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-gene-expr-unimputed-' + title + '.pdf')
    plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
############################################################################
# TC genes, RORgt cells
# select cells
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in [1, 4, 11, 20]])
geneList = [['Mki67', 'Gal', 'Nrxn1', 'Aire', 'Kif21a', 'Pigr', 'Col17a1', 'Hk2', 'LTb', 'Dnase1l3', 'Ahcyl2', 'Nlrc5', 'Itgb8', 'Ccl22', 'Ccl5', 'Il2ra']]
ncols = 3
for gene_list in geneList:
    for gene in gene_list:
        if gene.upper() not in genes:
            print(gene)
for i, title in zip(range(len(geneList)), ['TC']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-imputed-' + title + '.pdf')
    plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
# unimputed
for i, title in zip(range(len(geneList)), ['TC']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-unimputed-' + title + '.pdf')
    plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
############################################################################
# JC1/2/3 genes, RORgt cells
# select cells
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in [1, 4, 11, 20]])
ncols = 3
geneList = [
    ['Slc7a10', 'Dcaf12l2', 'Olig1', 'Gal', 'Atp1b1', 'Dsg1b', 'Ttyh1', 'Tbx5', 'Cnr1', 'Ank', 
     'Fam81a', 'B3galt1', 'Ube2e2', 'Syt1', 'Zfand6'], # JC1
    ['Egfl6', 'Tnni1', '1110008L16Rik', 'Cep112', 'Asic1', 'Ly9', 'Fabp1', 'Col17a1', 'Pgam2',
     'Poc1a', 'Clic3', 'Prdm16', 'Ppp2r2c', 'Gstt2'], # JC2
    ['Gm26917', 'Ptbp2', 'Zc3h7a', 'Lcor', 'Nfat5', 'Smg1', 'Cep350', 'Mdm4', 'Chuk', 
     'Mapk8ip3', 'Prpf39', 'Eml5', 'Phip', 'Rnf111', 'Trpm7'], # JC3
    ['Cxcr6', 'Clnk', 'Fam184b', 'Klrb1b', 'Klrb1f', 'Chad', 'Apol7e', 'Ncr1', 
     'Il22', 'Arg1', 'Il2rb', 'Dgat1', 'Il18rap', 'Gzmb', 'Ccdc184'] # ILC3
    ]
for gene_list in geneList:
    for gene in gene_list:
        if gene.upper() not in genes:
            print(gene)
# imputed
for i, title in zip(range(len(geneList)), ['JC1', 'JC2', 'JC3', 'ILC3']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-imputed-' + title + '.pdf')
    plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
# unimputed
for i, title in zip(range(len(geneList)), ['JC1', 'JC2', 'JC3', 'ILC3']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-unimputed-' + title + '.pdf')
    plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()

################################################################################
# sum of gene expression

geneList = [['Rorc', 'Rora', 'Cxcr6', 'Aire', "H2-Ab1"],
            ['Mki67', 'Gal', 'Nrxn1', 'Aire', 'Kif21a', 'Pigr', 'Col17a1', 'Hk2', 'LTb', 'Dnase1l3', 'Ahcyl2', 'Nlrc5', 'Itgb8', 'Ccl22', 'Ccl5', 'Il2ra'], # TC
            ['Slc7a10', 'Dcaf12l2', 'Olig1', 'Gal', 'Atp1b1', 'Dsg1b', 'Ttyh1', 'Tbx5', 'Cnr1', 'Ank', 'Fam81a', 'B3galt1', 'Ube2e2', 'Syt1', 'Zfand6'], # JC1
            ['Egfl6', 'Tnni1', '1110008L16Rik', 'Cep112', 'Asic1', 'Ly9', 'Fabp1', 'Col17a1', 'Pgam2', 
             'Poc1a', 'Clic3', 'Prdm16', 'Ppp2r2c', 'Gstt2'], # JC2
            ['Gm26917', 'Ptbp2', 'Zc3h7a', 'Lcor', 'Nfat5', 'Smg1', 'Cep350', 'Mdm4', 'Chuk', 
             'Mapk8ip3', 'Prpf39', 'Eml5', 'Phip', 'Rnf111', 'Trpm7'], # JC3
            ['Cxcr6', 'Clnk', 'Fam184b', 'Klrb1b', 'Klrb1f', 'Chad', 'Apol7e', 'Ncr1', 
             'Il22', 'Arg1', 'Il2rb', 'Dgat1', 'Il18rap', 'Gzmb', 'Ccdc184'] # ILC3]
            ]
def plot_multigene_expr(expr, vis, dim1, dim2, genes, file, figTitle, s=3):
    cmap = matplotlib.cm.Spectral_r
    fig, ax = plt.subplots(figsize=(7, 5))
    c = expr[list(map(lambda x: x.upper(), genes))].sum(axis = 1)
    ax.scatter(vis.loc[:, dim1], vis.loc[:, dim2], s=s, c=c, cmap=cmap)
    ax.set_axis_off()
    ax.set_title(figTitle)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    file.savefig(fig.figure, bbox_inches='tight')
    return 0
# ci = md.index
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in [1, 4, 11, 20]])
for i, title in zip(range(1, len(geneList)), ['TC', 'JC1', 'JC2', 'JC3', 'ILC3']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+-combined-gene-expr-unimputed-' + title + '.pdf')
    plot_multigene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                        file=file, figTitle=title)
    file.close()
for i, title in zip(range(1, len(geneList)), ['TC', 'JC1', 'JC2', 'JC3', 'ILC3']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+-combined-gene-expr-imputed-' + title + '.pdf')
    plot_multigene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                        file=file, figTitle=title)
    file.close()
# ci = md.index
# load TC signature genes
TC_signatures = pd.read_csv('../MLN_RORgt_MHCII_multiome/Seurat/results/markers/C2-5_top130.csv', header=0)
TC_I = TC_signatures.loc[TC_signatures.celltype == 'TC I']
TC_II = TC_signatures.loc[TC_signatures.celltype == 'TC II']
TC_III = TC_signatures.loc[TC_signatures.celltype == 'TC III']
TC_IV = TC_signatures.loc[TC_signatures.celltype == 'TC IV']
TC_SIG = [list(TC_I.gene_name), list(TC_II.gene_name), list(TC_III.gene_name), list(TC_IV.gene_name)]
for gene_list in TC_SIG:
    for gene in gene_list:
        if gene.upper() not in imp_df.columns:
            print(gene)
            gene_list.remove(gene)

for gene_list in TC_SIG:
    print(len(gene_list))

for i, title in zip(range(len(TC_SIG)), ['TC I', 'TC II', 'TC III', 'TC IV']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+-combined-gene-expr-imputed-' + title + '.pdf')
    plot_multigene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=TC_SIG[i],
                        file=file, figTitle=title)
    file.close()
    
for i, title in zip(range(len(TC_SIG)), ['TC I', 'TC II', 'TC III', 'TC IV']):
    print(i)
    file = PdfPages('plots/gene-expr/UMAP-RORgt+-combined-gene-expr-unimputed-' + title + '.pdf')
    plot_multigene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=TC_SIG[i],
                        file=file, figTitle=title)
    file.close()