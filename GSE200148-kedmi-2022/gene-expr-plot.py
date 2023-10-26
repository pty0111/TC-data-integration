import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt
import palantir as pl
from matplotlib.backends.backend_pdf import PdfPages
import os
matplotlib.rcParams['font.family'] = ['serif']

os.chdir('/Users/parkt1/0-workspace/TC-matters-arising/GSE200148-kedmi-2022')

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
# main markers
geneList = ['Rorc', 'Rora', 'Cxcr6', "Aire", "H2-Ab1"]

title = 'main_markers'
ncols = 2
ci = md.index

file = PdfPages('plots/gene-expr/UMAP-gene-expr-imputed-{}.pdf'.format(title))
plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList,
                file=file, n_cols=ncols)
file.close()
file = PdfPages('plots/gene-expr/UMAP-gene-expr-unimputed-{}.pdf'.format(title))
plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList,
                file=file, n_cols=ncols)
file.close()

############################################################################################################
# cells greyed out except for a set of cells
def plot_gene_expr_grayout(expr, vis, vis2, dim1, dim2, genes, file, n_cols, figTitle, s=3, 
                   newGeneName=None):
    cmap = matplotlib.cm.Spectral_r
    fig = pl.plot.FigureGrid(len(genes), n_cols)
    for g, ax in zip(genes, fig):
        c = expr.loc[vis.index, g.upper()]
        ax.scatter(vis2.loc[:, dim1], vis2.loc[:, dim2], s=s, c='gray')
        ax.scatter(vis.loc[:, dim1], vis.loc[:, dim2], s=s, c=c, cmap=cmap)
        ax.set_axis_off()
        if newGeneName and g in newGeneName:
            figTitle = newGeneName[g]
        else:
            figTitle = g
        ax.set_title(figTitle)
        normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    n_rows = math.ceil(len(genes) / n_cols)
    fig.figure.set_size_inches(7 * n_cols, 5 * n_rows)
    file.savefig(fig.figure, bbox_inches='tight')
    return 0

geneList = ['Mki67', 'Cxcr6',
             'Nrxn1', 'Gal', 
             'Col17a1', 'Hk2',
             'Dnase1l3', 'Nlrc5',
             'Itgb8', 'Ccl22']

clusters_to_keep = [1, 12, 16, 18]
clusters_to_grey = [10]
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_keep])
ci2 = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_grey])
cellgroup = 'RORgt+'
ncols = 2
title = 'Fig2c'

file = PdfPages('plots/gene-expr/UMAP-{}-gene-expr-imputed-{}.pdf'.format(cellgroup, title))
plot_gene_expr_grayout(expr=imp_df.loc[ci, :], 
                            vis=umap_s.loc[ci, :], vis2=umap_s.loc[ci2, :], 
                            dim1='UMAP_1', dim2='UMAP_2', genes=geneList, n_cols = ncols,
                            file=file, figTitle=title)
file.close()
file = PdfPages('plots/gene-expr/UMAP-{}-gene-expr-unimputed-{}.pdf'.format(cellgroup, title))
plot_gene_expr_grayout(expr=unimp_df.loc[ci, :], 
                            vis=umap_s.loc[ci, :], vis2=umap_s.loc[ci2, :], 
                            dim1='UMAP_1', dim2='UMAP_2', genes=geneList, n_cols = ncols,
                            file=file, figTitle=title)
file.close()
############################################################################################################
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

########## for TC module scores ############
clusters_to_keep = [1, 12, 16, 18]
clusters_to_grey = [10]
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_keep])
ci2 = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_grey])
cellgroup = 'RORgt+'

for colname, title in zip(['TC1', 'TC2', 'TC3', 'TC4'], ['TC I', 'TC II', 'TC III', 'TC IV']):
    print(colname)
    file = PdfPages('plots/gene-expr/UMAP-{}-combined-gene-expr-unimputed-{}.pdf'.format(cellgroup, title))
    plot_multigene_expr_grayout(expr=md.loc[ci, :], 
                                vis=umap_s.loc[ci, :], vis2=umap_s.loc[ci2, :], 
                                dim1='UMAP_1', dim2='UMAP_2', colName=colname,
                                file=file, figTitle=title, s=8)
    file.close()

########## for JC and ILC module scores ############
clusters_to_keep = [1, 10, 12, 16, 18]
clusters_to_grey = []
ci = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_keep])
ci2 = pd.Index([idx for idx, cl in zip(md.index, md['Clusters']) if cl in clusters_to_grey])
cellgroup = 'RORgt+'

for colname, title in zip(['JC1', 'JC2', 'JC3', 'ILC31'], ['JC1', 'JC2', 'JC3', 'ILC3']):
    print(colname)
    file = PdfPages('plots/gene-expr/UMAP-{}-combined-gene-expr-unimputed-{}.pdf'.format(cellgroup, title))
    plot_multigene_expr_grayout(expr=md.loc[ci, :], 
                                vis=umap_s.loc[ci, :], vis2=umap_s.loc[ci2, :], 
                                dim1='UMAP_1', dim2='UMAP_2', colName=colname,
                                file=file, figTitle=title, s=8)
    file.close()
 