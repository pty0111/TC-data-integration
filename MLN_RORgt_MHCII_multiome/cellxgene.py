import pandas as pd
import anndata as ad
import pyreadr as pr
import re
import os

os.chdir('/Users/pty0111/CCR7_DC/MLN_RORgt_MHCII_multiome/')


expr = pr.read_r('Seurat/results/imputed-expr.rds')[None]
idx_rna = expr.columns
idx_atac = '3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19#' + expr.columns

genes = pd.read_csv('cellranger-arc/outs/filtered_feature_bc_matrix/features.tsv', sep = '\t', index_col = 0,
                    header = None, names = ['gene_name', 'type', 'chr', 'start', 'end']).loc[expr.index, :]
genes['gene_name'] = genes['gene_name'] + '_' + genes.index

meta_data = pd.read_csv('Seurat/results/meta-data.csv', header = 0, index_col = 0).loc[idx_rna, :]
meta_data['Clusters_RNA'] = meta_data['Clusters'].astype("category")
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data.drop(['orig.ident', 'Clusters'], axis = 1, inplace = True)
md_atac = pd.read_csv('ArchR/ArchR_scRNA-seq_subset/cellColData.csv', header = 0, index_col = 0).loc[idx_atac, ]
md_atac.index = idx_rna
meta_data['nFrags_ATAC'] = md_atac['nFrags']; meta_data['PromoterRatio_ATAC'] = md_atac['PromoterRatio']
meta_data['Clusters_ATAC'] = pd.Series([int(re.sub('C', '', cl)) for cl in md_atac['Clusters']], index = idx_rna).astype("category")

UMAP_rna = pd.read_csv('Seurat/results/UMAP.csv', header = 0, index_col = 0).loc[idx_rna, :]
FDL_rna = pd.read_csv('Seurat/results/FDL.csv', header = 0, index_col = 0).loc[idx_rna, :]
UMAP_atac = pd.read_csv('ArchR/ArchR_scRNA-seq_subset/Embeddings/UMAP.csv', header = 0, index_col = 0).loc[idx_atac, :]
FDL_atac = pd.read_csv('ArchR/ArchR_scRNA-seq_subset/Embeddings/FDL.csv', header = 0, index_col = 0).loc[idx_atac, :]

ad_obj = ad.AnnData(
	X = expr.transpose(), var = genes, obs = meta_data,
	obsm = {
	  'X_fdl_RNA': FDL_rna.values, 'X_umap_RNA' : UMAP_rna.values,
	  'X_fdl_ATAC': FDL_atac.values, 'X_umap_ATAC' : UMAP_atac.values
	}
)
ad_obj.uns['Clusters_RNA_colors'] = ['#324E72', '#A4E804', '#CB7E98', '#0089A3', '#404E55',
                                     '#FDE8DC', '#5B4534', '#922329', '#3A2465', '#99ADC0',
                                     '#BC23FF', '#72418F', '#201625', '#FFF69F', '#549E79',
                                     '#9B9700', '#D790FF', '#772600', '#6B002C', '#A05837',
                                     '#6367A9']
ad_obj.uns['Clusters_ATAC_colors'] = ['#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941',
                                      '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6',
                                      '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87',
                                      '#5A0007', '#809693', '#6A3A4C']
ad_obj.write_h5ad("Seurat/results/cellxgene.h5ad")


expr = pr.read_r('Seurat/results_before_filtering_by_scATAC-seq/imputed-expr.rds')[None]
idx_rna = expr.columns

genes = pd.read_csv('cellranger-arc/outs/filtered_feature_bc_matrix/features.tsv', sep = '\t', index_col = 0,
                    header = None, names = ['gene_name', 'type', 'chr', 'start', 'end']).loc[expr.index, :]
genes['gene_name'] = genes['gene_name'] + '_' + genes.index

meta_data = pd.read_csv('Seurat/results_before_filtering_by_scATAC-seq/meta-data.csv', header = 0, index_col = 0).loc[idx_rna, :]
meta_data['Clusters_RNA'] = meta_data['Clusters'].astype("category")
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data.drop(['orig.ident', 'Clusters'], axis = 1, inplace = True)

UMAP_rna = pd.read_csv('Seurat/results_before_filtering_by_scATAC-seq/UMAP.csv', header = 0, index_col = 0).loc[idx_rna, :]
FDL_rna = pd.read_csv('Seurat/results_before_filtering_by_scATAC-seq/FDL.csv', header = 0, index_col = 0).loc[idx_rna, :]

ad_obj = ad.AnnData(
	X = expr.transpose(), var = genes, obs = meta_data,
	obsm = {
	  'X_fdl_RNA': FDL_rna.values, 'X_umap_RNA' : UMAP_rna.values
	}
)
ad_obj.uns['Clusters_RNA_colors'] = ['#324E72', '#A4E804', '#CB7E98' , '#0089A3', '#404E55',
                                     '#FDE8DC', '#5B4534', '#922329' , '#3A2465', '#99ADC0',
                                     '#BC23FF', '#72418F', '#201625' , '#FFF69F', '#549E79',
                                     '#9B9700', '#D790FF', '#772600' , '#6B002C', '#A05837',
                                     '#6367A9', '#A77500', '#7900D7']
ad_obj.write_h5ad("Seurat/results_before_filtering_by_scATAC-seq/cellxgene_before_filtering_by_scATAC-seq.h5ad")


########################################################
# integrate with Lyu
meta_data = pd.read_csv('Seurat/results/integrated-with-Lyu/meta-data.csv', header=0, index_col=0)
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data['Clusters'] = pd.Categorical(meta_data['Clusters'])
header = ['nCount_RNA', 'nFeature_RNA', 'MtFrac_RNA', 'S.Score', 'G2M.Score', 
          'Phase', 'Clusters']
meta_data = meta_data.loc[:, header]
# meta_data = meta_data.iloc[:, [1, 2, 7, 9, 10, 11, 6, 8]]

UMAP = pd.read_csv('Seurat/results/integrated-with-Lyu/UMAP.csv', header=0, index_col=0)

# # unimputed
# expr = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0).transpose()
# expr.index = [ ind.replace(".", "-") for ind in expr.index ]
# ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
# ad_obj.write_h5ad('results/unimputed.h5ad')

# imputed
expr = pr.read_r('Seurat/results/integrated-with-Lyu/imputed-expr.rds')[None]
# expr = pd.read_csv('integrated-with-wang/results/imputed-expr.csv', header=0, index_col=0).transpose()
expr = expr.transpose()
ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
ad_obj.write_h5ad('Seurat/results/integrated-with-Lyu/imputed.h5ad')

########################################################
# integrate with Kedmi and Wang
meta_data = pd.read_csv('Seurat/results/integrated-with-Lyu/meta-data.csv', header=0, index_col=0)
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data['Clusters'] = pd.Categorical(meta_data['Clusters'])
header = ['nCount_RNA', 'nFeature_RNA', 'MtFrac_RNA', 'S.Score', 'G2M.Score', 
          'Phase', 'Clusters']
meta_data = meta_data.loc[:, header]
# meta_data = meta_data.iloc[:, [1, 2, 7, 9, 10, 11, 6, 8]]

UMAP = pd.read_csv('Seurat/results/integrated-with-Lyu/UMAP.csv', header=0, index_col=0)

# # unimputed
# expr = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0).transpose()
# expr.index = [ ind.replace(".", "-") for ind in expr.index ]
# ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
# ad_obj.write_h5ad('results/unimputed.h5ad')

# imputed
expr = pr.read_r('Seurat/results/integrated-with-Lyu/imputed-expr.rds')[None]
# expr = pd.read_csv('integrated-with-wang/results/imputed-expr.csv', header=0, index_col=0).transpose()
expr = expr.transpose()
ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
ad_obj.write_h5ad('Seurat/results/integrated-with-Lyu/imputed.h5ad')

########################################################
# integrate with Kedmi Wang Lyu
meta_data = pd.read_csv('Seurat/results/integrated-with-Kedmi-Wang-Lyu/meta-data.csv', header=0, index_col=0)
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data['Clusters'] = pd.Categorical(meta_data['Clusters'])
header = ['nCount_RNA', 'nFeature_RNA', 'MtFrac_RNA', 'S.Score', 'G2M.Score', 
          'Phase', 'Clusters']
meta_data = meta_data.loc[:, header]
# meta_data = meta_data.iloc[:, [1, 2, 7, 9, 10, 11, 6, 8]]

UMAP = pd.read_csv('Seurat/results/integrated-with-Kedmi-Wang-Lyu/UMAP.csv', header=0, index_col=0)

# # unimputed
# expr = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0).transpose()
# expr.index = [ ind.replace(".", "-") for ind in expr.index ]
# ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
# ad_obj.write_h5ad('results/unimputed.h5ad')

# imputed
expr = pr.read_r('Seurat/results/integrated-with-Kedmi-Wang-Lyu/imputed-expr.rds')[None]
# expr = pd.read_csv('integrated-with-wang/results/imputed-expr.csv', header=0, index_col=0).transpose()
expr = expr.transpose()

gene_info = pd.read_csv("gene.info.csv")
id_to_symbol = dict(zip(gene_info.iloc[:,0], gene_info.iloc[:,1].str.upper()))

# there are duplicate gene names, so add id
new_columns = [id_to_symbol[g].upper()+"_"+g.upper() if g in id_to_symbol else g.upper() for g in expr.columns]
import numpy as np
expr.columns = new_columns

ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
ad_obj.write_h5ad('Seurat/results/integrated-with-Kedmi-Wang-Lyu/imputed.h5ad')