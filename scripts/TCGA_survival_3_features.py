import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
                 

# http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival5.html Explanation of logrank test 


def get_TCGA_clinical_data(ind_list):

  all_ind =[]
  for indication in ind_list:
    if indication in ['GYN']:
      df = get_gyn_df()
    else:  
      df = pd.read_csv('TCGA/'+indication+'/'+indication+'_clinical_data.tsv',sep='\t',index_col=[0],header=0)
    
    all_ind.append(df)
  df = pd.concat(all_ind,axis=0)
  df['E'] = df['vital_status'].replace(regex={r'Alive':0,'Dead':1})
  
  return df

def get_gyn_df():
  gyn_all =[]
  for indication in ['OV','UCEC','UCS']:
    df = pd.read_csv('TCGA/'+indication+'/'+indication+'_clinical_data.tsv',sep='\t',index_col=[0],header=0)
    gyn_all.append(df)
  return pd.concat(gyn_all,axis=0)

def make_archetype_dict(df,df_cluster):
  arch_name ={
    1: 'Immune Rich',
    2: 'Immune Stromal Rich',
    3: 'Tcell Centric',
    4: 'Myeloid Centric',
    5: 'Immune Desert',
    6: 'Immune Stromal Desert'
  }
  
  cluster_dict = {}
  for clust in [1, 2, 3, 4, 5, 6]:
    arch_samples = df_cluster[df_cluster['louvain']==clust].index.values
    df_c = df.filter(arch_samples,axis=0)
    df_c['group'] = arch_name[clust]
    #df_c = df_c.join(df_cluster[['Sample.1','HNSC','KIRC','SKCM','BLCA','SARC','GYN','COAD']])
    cluster_dict[arch_name[clust]] = df_c
  return cluster_dict

def get_TCGA_archetypes():
  df = pd.read_csv('../files_used_for_plots/TCGA_archetypes.tsv',sep='\t',index_col=[0],header=0)
  return df

def get_survival_data(cluster_dict):
  s_times=[]
  for arch in [
               'Immune Rich',
               'Immune Stromal Rich',
               'Tcell Centric',
               'Myeloid Centric',
               'Immune Desert',
               'Immune Stromal Desert']:
  
  
    if not cluster_dict[arch].empty:
      c = cluster_dict[arch][['OS.time','E','group']]
      c = c.dropna(how='any')

      
    else:
      c = pd.DataFrame()
    s_times.append(c)
      
  survival_times = pd.concat(s_times,axis=0)

  return survival_times

def main():
  
  # Add indication/all_indication - name used in the file names
  indication_name ='BLCA'
  # indciation list options ['HNSC','KIRC','SKCM','BLCA','SARC','GYN','COAD']
  ind_list =['BLCA']
 
  df = get_TCGA_clinical_data(ind_list)
  df_cluster = get_TCGA_archetypes()
  cluster_dict = make_archetype_dict(df,df_cluster)
  survival_times = get_survival_data(cluster_dict)
 


  T = survival_times['OS.time']
  E = survival_times['E']
  groups = survival_times['group']
  ir = (groups == 'Immune Rich')
  isr = (groups == 'Immune Stromal Rich') 
  tc = (groups == 'Tcell Centric')
  my =(groups =='Myeloid Centric')
  imd = (groups == 'Immune Desert')
  isd = (groups == 'Immune Stromal Desert') 
  
  kmf = KaplanMeierFitter()


  kmf.fit(T[ir], E[ir], label='Immune Rich',alpha=1)
  ax= kmf.plot(color='r',show_censors=True,censor_styles={'ms':6,'marker':'|'})

  kmf.fit(T[isr], E[isr], label='Immune Stromal Rich',alpha=1)
  ax= kmf.plot(ax=ax ,color = 'g',show_censors=True,censor_styles={'ms':6,'marker':'|'})

  kmf.fit(T[tc], E[tc], label='Tcell Centric',alpha=1)
  ax= kmf.plot(ax=ax,color = 'b',show_censors=True,censor_styles={'ms':6,'marker':'|'})
 
    
  kmf.fit(T[my], E[my], label='Myeloid Centric',alpha=1)
  ax= kmf.plot(ax=ax,color='grey',show_censors=True,censor_styles={'ms':6,'marker':'|'})
 
    
  kmf.fit(T[imd], E[imd], label='Immune Desert',alpha=1)
  ax= kmf.plot(ax=ax,color ='yellow',show_censors=True,censor_styles={'ms':6,'marker':'|'})
  

  kmf.fit(T[isd], E[isd], label='Immune Stromal Desert',alpha=1)
  ax= kmf.plot(ax=ax,color='orange',show_censors=True,censor_styles={'ms':6,'marker':'|'})


  results = logrank_test(T[~ir],T[ir], E[~ir],E[ir], alpha=.99)
  results.print_summary()
  print (results.p_value)

  plt.title('5 Year Survival Curve for TCGA '+indication_name+' Immune Clusters',fontweight='bold')
  plt.text(1000,1, 'P-value: %.2E'%results.p_value, horizontalalignment='center',verticalalignment='center')#,fontweight='bold')
  plt.ylim(0,1.19)
  plt.xlim(0,2000)

  plt.savefig('TCGA_'+indication_name+'_immune_clusters_5yr_overall_survival.pdf',dpi=500,format='pdf')

if __name__ == '__main__':
     main() 