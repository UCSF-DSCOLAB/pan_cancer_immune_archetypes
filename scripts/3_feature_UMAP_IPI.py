
# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
import pandas as pd
import anndata

sc.settings.verbosity = 3
sc.logging.print_versions()
sc.set_figure_params(dpi=200)
sc.settings.autoshow=False

results_file = './write/TcStMyTr_features.h5ad'


def get_indication(df):
  indication =[]
  for sample in df.index:
    name = sample[sample.rfind('IPI')+3:sample.find('.')-3]
    if name == 'MELB':
      name = 'MEL'
    #print '%s\t%s'%(name, sample)  
    indication.append(name.strip())
  return indication

def get_anndata_formats(df,df_features):

  features = []

  X = df.values
  indications = get_indication(df)
  
  obs = df.index
  obs = obs.rename('Sample')
  obs = obs.to_frame()
  obs['indication'] = indications
  for feature in df_features.columns:
    obs[feature] = pd.to_numeric(df_features[feature],errors='coerce')
    features.append(feature)

  var = df.columns
  var = var.to_frame(name="Populations")
  var['ind'] = range(len(var))
  var = var.set_index(var['ind'])
  print (var)
  
  return (X,obs,var,features)


  
def main():
  
  # load up data
  df = pd.read_csv('../files_used_for_plots/feature_matrix.txt', sep='\t',index_col=[0], header=0)
  
  #add title for output files
  title ='Title4'
  print (df)

  df = df.T
  
  #filter out Gross Immune populations (df_GI) and also filter other features to overlay

  df_GI = df[['Tcell','Myeloid','CD90+ CD44+ Stroma']]
  df_GI= df_GI.dropna(axis ='rows',how='any')
  df_GI = df_GI.apply(pd.to_numeric,errors='coerce')

  df_features = df.drop(columns=['Tcell','Myeloid','CD90+ CD44+ Stroma'])

  print(df_features)
  # make the features the same size as main clustering features
  df_features = df_features.filter(df_GI.index,axis=0)
  print (df_GI)
  X, obs, var ,features= get_anndata_formats(df_GI,df_features)

  spread = np.std(X) # optimal value for spread in tl.umap https://github.com/lmcinnes/umap/issues/158
  adata = anndata.AnnData(X=X,obs=obs,var=var)
  adata.var_names = var['Populations']

  print (adata.var.index)

  
  for n in [34]:
    print (n)
    sc.tl.pca(adata,svd_solver='arpack')
    sc.pp.neighbors(adata,n_neighbors=n)#,n_pcs=2)
    sc.tl.umap(adata, min_dist=0.25,maxiter=2500,spread=spread)
    sc.tl.louvain(adata,use_weights=True,resolution=0.4)

    # put clusters in required order and re label from 1-6
    adata.obs['louvain']  = adata.obs['louvain'].replace(['0','1','2','3','4','5'],['2','0','4','1','3','5'])# louvain param TcMySt
    adata.obs['louvain']  = adata.obs['louvain'].map({'0':'1','1':'2','2':'3','3':'4','4':'5' ,'5':'6'})
   
    #print outputs
    for file_format in ['pdf','svg']:
      #louvain clusters
      sc.pl.umap(adata,
                 color=['louvain'],save='_'+title+'_louvain_'+str(n)+'.'+file_format,
                 palette= ['Red','Green','MidnightBlue','Grey','Gold','DarkOrange','Black'])

      #clusters by indication 
      sc.pl.umap(adata,
                 color=['indication'],save='_'+title+'_IND_'+str(n)+'.'+file_format,
                 palette=['Gold','DarkTurquoise','Grey','Coral','DodgerBlue',
                          'OliveDrab','Red','MidnightBlue','Maroon','Yellow','DarkOrchid','DeepPink','Black'])
      
      for feature in features:
        sc.pl.umap(adata,
                   color=[feature],
                   save="_"+title+'_'+str(n)+'_'+feature+'.'+file_format)#,size=800)

      pops = ['Tcell','Myeloid','CD90+ CD44+ Stroma']
     
      sc.pl.violin(adata,
                        pops,
                        groupby='louvain',
                        save='_'+title+'_'+str(n)+'.'+file_format,
                        palette=['Red','Green','MidnightBlue','Grey','Gold','DarkOrange'])

      for pop in pops:
         sc.pl.umap(adata,
                    color=[pop],
                    save='_'+title+'_'+pop+'_'+str(n)+'.'+file_format)


    adata.write_csvs(str(n)+'_'+title+'.csv',skip_data=False)#(results_file)

if __name__ == '__main__':
     main()            

