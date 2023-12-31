import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt
import palantir as pl
from matplotlib.backends.backend_pdf import PdfPages
import os
matplotlib.rcParams['font.family'] = ['serif']

os.chdir('/Users/parkt1/0-workspace/TC-matters-arising/GSE176282-wang-2021')

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
umap = pd.read_csv('results/UMAP.csv', header=0, index_col=0)
umap_s = pd.read_csv('results/UMAP-subset.csv', header=0, index_col=0)
imp_df = pd.read_csv('results/imputed-expr.csv', header=0, index_col=0)
imp_df = imp_df.transpose()
unimp_df = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0)
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

for i, title in zip(range(len(geneList)), ['main_markers']):
    # imputed
    file = PdfPages('plots/gene-expr/UMAP-gene-expr-imputed-' + title + '.pdf')
    plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
    # unimputed
    file = PdfPages('plots/gene-expr/UMAP-gene-expr-unimputed-' + title + '.pdf')
    plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()

############################################################################
# TC, JC1/2/3, ILC3 genes, RORgt cells
# select cells
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in [1, 4, 11, 20]])
ncols = 3
geneList = [
    ['Mki67', 'Gal', 'Nrxn1', 'Aire', 'Kif21a', 'Pigr', 'Col17a1', 'Hk2', 'LTb', 
     'Dnase1l3', 'Ahcyl2', 'Nlrc5', 'Itgb8', 'Ccl22', 'Ccl5', 'Il2ra'], # TC
    ['Slc7a10', 'Dcaf12l2', 'Olig1', 'Gal', 'Atp1b1', 'Dsg1b', 'Ttyh1', 'Tbx5', 'Cnr1', 'Ank', 
     'Fam81a', 'B3galt1', 'Ube2e2', 'Syt1', 'Zfand6'], # JC1
    ['Egfl6', 'Tnni1', '1110008L16Rik', 'Cep112', 'Asic1', 'Ly9', 'Fabp1', 'Col17a1', 'Pgam2',
     'Poc1a', 'Clic3', 'Prdm16', 'Ppp2r2c', 'Gstt2'], # JC2
    ['Gm26917', 'Ptbp2', 'Zc3h7a', 'Lcor', 'Nfat5', 'Smg1', 'Cep350', 'Mdm4', 'Chuk', 
     'Mapk8ip3', 'Prpf39', 'Eml5', 'Phip', 'Rnf111', 'Trpm7'], # JC3
    ['Cxcr6', 'Clnk', 'Fam184b', 'Klrb1b', 'Klrb1f', 'Chad', 'Apol7e', 'Ncr1', 
     'Il22', 'Arg1', 'Il2rb', 'Dgat1', 'Il18rap', 'Gzmb', 'Ccdc184'] # ILC3
    ]

for i, title in zip(range(len(geneList)), ['TC', 'JC1', 'JC2', 'JC3', 'ILC3']):
    # imputed
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-imputed-' + title + '.pdf')
    plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()
    # unimputed
    file = PdfPages('plots/gene-expr/UMAP-RORgt+cells-gene-expr-unimputed-' + title + '.pdf')
    plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList[i],
                   file=file, n_cols=ncols)
    file.close()

################################################################################
# sum of gene expression
# cells greyed out except for a set of cells
def plot_multigene_expr_grayout(expr, vis, vis2, dim1, dim2, colName, file, figTitle, s=3):
    cmap = matplotlib.cm.Spectral_r
    fig, ax = plt.subplots(figsize=(7, 5))
    c = expr[colName]
    ax.scatter(vis2.loc[:, dim1], vis2.loc[:, dim2], s=s, c='gray')
    ax.scatter(vis.loc[:, dim1], vis.loc[:, dim2], s=s, c=c, cmap=cmap)
    ax.set_axis_off()
    ax.set_title(figTitle)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    file.savefig(fig.figure, bbox_inches='tight')
    return 0

########## for JC and ILC module scores ############
clusters_to_keep = [1, 4, 11, 20]
clusters_to_grey = []
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_keep])
ci2 = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_grey])

cellgroup = 'RORgt+'
for colname, title in zip(['TC1', 'TC2', 'TC3', 'TC4', 'JC1', 'JC2', 'JC3', 'ILC31'], 
                          ['TC I', 'TC II', 'TC III', 'TC IV', 'JC1', 'JC2', 'JC3', 'ILC3']):
    print(colname)
    file = PdfPages('plots/gene-expr/UMAP-{}-combined-gene-expr-unimputed-{}.pdf'.format(cellgroup, title))
    plot_multigene_expr_grayout(expr=md.loc[ci, :], 
                                vis=umap_s.loc[ci, :], vis2=umap_s.loc[ci2, :], 
                                dim1='UMAP_1', dim2='UMAP_2', colName=colname,
                                file=file, figTitle=title, s=8)
    file.close()
