import pandas as pd
import matplotlib.pyplot as plt
from lifelines.statistics import multivariate_logrank_test
from lifelines import CoxPHFitter

# http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival5.html Explanation of logrank test 

def get_gyn_df():
  gyn_all =[]
  for indication in ['OV','UCEC','UCS']:
    df = pd.read_csv('/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/TCGA/'+indication+'/'+indication+'_clinical_data.tsv',sep='\t',index_col=[0],header=0)
    gyn_all.append(df)
  return pd.concat(gyn_all,axis=0)

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

def get_TCGA_archetypes():
  df = pd.read_csv('../files_used_for_plots/TCGA_archetypes_assigned.tsv',sep='\t',index_col=[0],header=0)
  return adjust_categorical_data(df,'Indication')


def adjust_categorical_data(df,category):
  ind_dummies = pd.get_dummies(df[category])
  return df.join(ind_dummies)

def make_archetype_dict(df,df_cluster,cutoff):
  cluster_dict = {}
  for clust in ['red','pink','maroon','green','light_green','blue','light_blue','grey','black','yellow','khaki','orange']:
    arch_samples = df_cluster[df_cluster['group']==clust].index.values
    df_c = df.filter(arch_samples,axis=0)
    df_c = df_c.join(df_cluster[['Indication','HNSC','KIRC','SKCM','BLCA','SARC','GYN','COAD','LIHC','LUAD','PAAD','GBM']])

    if len(df_c) < cutoff:
      df_c = pd.DataFrame()
    cluster_dict[clust] = df_c
  return cluster_dict

def get_survival_data(cluster_dict):
  s_times=[]
  for arch in [
               ('red','Immune Rich CD8 Macro Bias'),
               ('pink','Immune Rich CD8 Mono Bias'),
               ('maroon', 'Immune Rich CD4 Macro Bias'),
               ('green','Immune Stromal CD8 Bias'),
               ('light_green', 'Immune Stromal CD4 Macro Bias'),
               ('blue','Tcell Centric Macro Bias'),
               ('light_blue','Tcell Centric DC Rich'),
               ('grey','Myeloid Centric DC2 Bias'),
               ('black','Myeloid Centric DC1 Bias'),
               ('yellow', 'Immune Desert CD4 Macro Bias'),
               ('khaki','Immune Desert Mono Bias'),
               ('orange', 'Immune Desert CD8 Macro Bias')]:
  
  
    if not cluster_dict[arch[0]].empty:
      cluster_dict[arch[0]]['group'] = arch[1]
      c = cluster_dict[arch[0]][['OS.time','E','group','Indication','HNSC','KIRC','SKCM','BLCA','SARC','GYN','COAD','LIHC','LUAD','PAAD','GBM']]
      c = c.dropna(how='any')
      print (c)
      
    else:
      c = pd.DataFrame()
    s_times.append(c)
      
  survival_times = pd.concat(s_times,axis=0)
  return survival_times

def main():
  
  # Add indication/all_indication - name used in the file names
  indication_name ='All_indication'
  # indciation list 
  ind_list =['HNSC','KIRC','SKCM','BLCA','SARC','GYN','COAD','LIHC','LUAD','PAAD','GBM']
  df = get_TCGA_clinical_data(ind_list)
  df_cluster = get_TCGA_archetypes()
  cluster_dict = make_archetype_dict(df,df_cluster,15)
  survival_times = get_survival_data(cluster_dict)


  T = survival_times['OS.time']
  E = survival_times['E']
  groups = survival_times['group']
  irr = (groups == 'Immune Rich CD8 Macro Bias')
  irp = (groups == 'Immune Rich CD8 Mono Bias')
  irm = (groups == 'Immune Rich CD4 Macro Bias')  
  isrg = (groups == 'Immune Stromal CD8 Bias')
  isrlg = (groups == 'Immune Stromal CD4 Macro Bias')  
  tcb = (groups == 'Tcell Centric Macro Bias')
  tclb = (groups == 'Tcell Centric DC Rich')
  myg =(groups =='Myeloid Centric DC2 Bias')
  myb =(groups =='Myeloid Centric DC1 Bias')
  imdy = (groups == 'Immune Desert CD4 Macro Bias')
  imdk = (groups == 'Immune Desert Mono Bias')
  imdo = (groups == 'Immune Desert CD8 Macro Bias') 
  

  cph = CoxPHFitter()
  arch = []
  first =0
  results_df_list=[]
  for arch_s in [
                 (irr,'Immune Rich CD8 Macro Bias','#ff0000'),
                 (irp,'Immune Rich CD8 Mono Bias','#ff69B4'),
                 (irm, 'Immune Rich CD4 Macro Bias','#800000'),
                 (isrg,'Immune Stromal CD8 Bias','#008000'),
                 (isrlg, 'Immune Stromal CD4 Macro Bias','#9acd32'),
                 (tcb,'Tcell Centric Macro Bias','#191970'),
                 (tclb,'Tcell Centric DC Rich','#1e90ff'),
                 (myg,'Myeloid Centric DC2 Bias','#808080'),
                 (myb,'Myeloid Centric DC1 Bias','#000000'),
                 (imdy, 'Immune Desert CD4 Macro Bias','#ffff00'),
                 (imdk,'Immune Desert Mono Bias','#c0ad8c'),
                 (imdo, 'Immune Desert CD8 Macro Bias','#ffa500')]:  


      arch.append((arch_s[0],arch_s[1]))
      #print (arch_s)
      data =survival_times[arch_s[0]]
      results_df_list.append(data)
      data = data.drop(['group','Indication'],axis=1)
      data = data.loc[:, (data != 0).any(axis=0)] #remove columns with all zeroes
      drop_list = [col for col, val in data.sum().iteritems() if val < 15 and col !='E']
      data = data.drop(drop_list, axis=1)
      
      if arch_s[0].any():
        if first == 0:
          cph.fit(data,duration_col='OS.time', event_col='E',step_size=(0.1),show_progress=True)
          cph.print_summary()
          ax =cph.baseline_survival_.plot(color=arch_s[2],label=arch_s[1])
          first = 1
        else:
          cph.fit(data,duration_col='OS.time', event_col='E',step_size=(0.1), show_progress=True)
          cph.print_summary()      
          ax =cph.baseline_survival_.plot(ax=ax,color=arch_s[2],label=arch_s[1])

  results_df = pd.concat(results_df_list,axis=0)
  T = results_df['OS.time']
  E = results_df['E']
  groups = results_df['group']  
  results = multivariate_logrank_test(T,groups,E, alpha=.99)
  results.print_summary()
  print (results.p_value)
  

  L= plt.legend()
  for i,n in enumerate(arch):
    L.get_texts()[i].set_text(arch[i][1])
  plt.title('5 Year Survival Curve for TCGA '+indication_name+' \n Immune Archetype Signatures',fontweight='bold')
  plt.text(300,1.1, 'P-value: %.2E'%results.p_value, horizontalalignment='center',verticalalignment='center')#,fontweight='bold')
  plt.ylim(0,1.19)
  plt.xlim(0,2000)
 
  #plt.legend(loc='upper right',fontsize=5)
  plt.savefig('TCGA_'+indication_name+'_immune_assigned_archetype_sig_5yr_overall_survival.pdf',dpi=500,format='pdf')
  plt.close()
  
             
      #plt.show()

if __name__ == '__main__':
     main() 