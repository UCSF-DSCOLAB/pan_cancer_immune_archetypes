
# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from sklearn.metrics import davies_bouldin_score
from itertools import product

sc.settings.verbosity = 3
#sc.logging.print_versions()
sc.set_figure_params(dpi=200)
sc.settings.autoshow=False

def get_score(df,f):
  return df.values.tolist()


def get_anndata_formats(df):

  df = df.set_index(df['Populations'])
  df =df.drop(columns=['Populations'])
  df =df.T
  
  X =df.values
  X = X.astype(np.float)
  
  obs = df.index
  obs = obs.rename('Sample')
  obs = obs.to_frame()


  var = df.columns
  var = var.to_frame()
  var['ind'] = range(len(var))
  var = var.set_index(var['ind'])
  var = var.drop(columns=['ind'])
  #print (var)
  
  return (X,obs,var)

def get_DBI_score(title,X,adata,spread,resolutions,neighbors):
  '''

  Parameters
  ----------
  title : String for naming outfiles
  X : feature values
  adata : anndata object
  spread : UMAP param that scale embeddings
  resolutions : UMAP resolutions list
  neighbors : list of nearest neighbors

  Returns
  -------
  "*_UMAP_cluster_size_counts.tsv" the number of time a cluster size is found
  "*_UMAP_cluster_frequency_scores.tsv" DBI , cluster size and resolution of each test
  "*_UMAP_cluster_num_DBImin.tsv" Minimun DBI of each cluster size seen

  '''

  counts ={}
  for res, n in product(resolutions,neighbors):         
    print (n)
    sc.tl.pca(adata,svd_solver='arpack')
    sc.pp.neighbors(adata,n_neighbors=n)#,n_pcs=2)
    sc.tl.umap(adata, min_dist=0.25,maxiter=2500,spread=spread)
    sc.tl.louvain(adata,use_weights=True,resolution=res)
    print ('clusters = %d'%len(set(adata.obs['louvain'])))
    try:
      dbi= davies_bouldin_score(X,adata.obs['louvain'] )
      print ('DBI = %f'%dbi)

    except:
      dbi=100

    counts[(n,res)] = (len(set(adata.obs['louvain'])),dbi,res)

  df_counts=pd.DataFrame.from_dict(counts,orient='index',columns=['Clusters','DBI','Resolution'])
  cluster_counts =df_counts.groupby(['Clusters']).size()
  
  cluster_counts.to_csv(title+'_UMAP_cluster_size_counts.tsv',sep='\t',header=True)
  df_counts.to_csv(title+'_UMAP_cluster_frequency_scores.tsv',sep='\t',header=True)
  
  df_cluster_dbi = df_counts.loc[df_counts.groupby('Clusters')['DBI'].idxmin()]
  df_cluster_dbi.to_csv(title+'_UMAP_cluster_num_DBImin.tsv',sep='\t',header=True)


def main():
  
  # add a title for your output files
  title = "title"
 
  #add a list of resolutions to test
  resolutions = [0.3,0.5,1]
 
  #add the range of nearest neighbors to test
  neighbors = range(3,10)
  
  # load up data
  df = pd.read_csv('../files_used_for_plots/feature_matrix.txt', sep='\t',index_col=[0], header=0)
  print (df.shape)

  df = df.T
  df_GI = df[['Tcell','Myeloid','CD90+ CD44+ Stroma','Treg','CD4_calculated_score','CD8_calculated_score','Classical_mono_calculated_score','Macro_calculated_score','cDC1_calculated_score','cDC2_calculated_score']]

  df_GI= df_GI.dropna(axis ='rows',how='any')
  df_GI = df_GI.T

  df_GI['Populations'] = df_GI.index
  df_GI['C'] = np.arange(len(df_GI))
  df_GI = df_GI.set_index('C')
 
  X, obs, var = get_anndata_formats(df_GI)
  spread = (np.std(X)) # optimal value for spread in tl.umap https://github.com/lmcinnes/umap/issues/158
  adata = anndata.AnnData(X=X,obs=obs,var=var)
  adata.var_names = var['Populations']
  
  get_DBI_score(title,X,adata,spread,resolutions,neighbors)
  


if __name__ == '__main__':
     main()            

