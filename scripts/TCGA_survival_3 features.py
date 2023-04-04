import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test,logrank_test
import glob
import os
from lifelines.utils import median_survival_times                  

# http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival5.html Explanation of logrank test 
def get_high_and_low_samples(df):
  num =int(df.shape[0])/4
  high_df = df.nlargest(num,'Score')
  high = high_df.index.values.tolist()
  low_df = df.nsmallest(num,'Score')
  low = low_df.index.values.tolist()
  print len(low)
  print len(high)
  return (high,low)

def get_map_names():
  
  tot =[]
  inpath ="/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/PanCan_new_data/TCGA/TCGA_name_maps"
  for infile in glob.glob(os.path.join(inpath,  '*_name_map.tsv')):
    df_map = pd.read_csv(infile, sep='\t',index_col=[0],header=0)
    df_map = df_map.set_index('new_name')
    tot.append(df_map)
  return pd.concat(tot,axis=0)

def main():
  all_ind =[]
  for indication in ['HNSC','KIRC','SKCM','BLCA','SARC','OV','UCS','UCEC','COAD','LIHC','LUAD','PAAD','GBM']:
    df = pd.read_csv('/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/TCGA/'+indication+'/'+indication+'_clinical_data_20171025_xena.tsv',sep='\t',index_col=[0],header=0)
    #df = df[df['OS.time']<= 2000]
    #print indication
    #print df.shape
    all_ind.append(df)

  

    df['E'] = df['vital_status'].replace(regex={r'Alive':0,'Dead':1})
    print df.shape
    indication = 'All_Indications'
    tot_medians ={}
    #for indication in ['BLCA','COAD','GYN','HNSC','KIRC', 'SARC','SKCM']:
    df_cluster = pd.read_csv('/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/PanCan_new_data/TSNE/925_TCGA_TcMySt_logCPM_April29_2020_res_0.5/obs.csv',sep=',',index_col=[0],header=0)
    #df_cluster = df_cluster.filter(regex=indication,axis=0)
    print df_cluster
    df_cluster = df_cluster['louvain']
    df_cluster = df_cluster.to_frame()
    df_map = get_map_names()
    cluster_dict = {}
    for clust in [1, 2, 3, 4, 5, 6]:
      #print df_cluster['louvain']
      cluster = df_cluster.loc[df_cluster['louvain'] == clust]
      print cluster.shape
      names = df_map.filter(cluster.index.values, axis=0)
      names = names.set_index('tcga')
      cluster_names = names.index.values
      df_c = df.filter(cluster_names,axis=0)
      print df_c
      cluster_dict[clust] = df_c
  
     
    cluster_dict[1]['group'] = 'Immune Rich'
    imm_rich= cluster_dict[1][['OS.time','E','group']]
    imm_rich = imm_rich.dropna(how='any')
    print imm_rich
    
    cluster_dict[3]['group'] = 'Tcell Centric'
    tcell= cluster_dict[3][['OS.time','E','group']]
    tcell = tcell.dropna(how='any')
    print tcell
    
    cluster_dict[5]['group'] = 'Immune Desert'
    imm_des= cluster_dict[5][['OS.time','E','group']]
    imm_des = imm_des.dropna(how='any')
    print imm_des
    
    cluster_dict[2]['group'] = 'Immune Stromal Rich'
    imm_st_rich= cluster_dict[2][['OS.time','E','group']]
    imm_st_rich = imm_st_rich.dropna(how='any')
    print imm_st_rich
    
    cluster_dict[4]['group'] = 'Myeloid Centric'
    myeloid = cluster_dict[4][['OS.time','E','group']]
    myeloid = myeloid.dropna(how='any')
    print myeloid
    
    cluster_dict[6]['group'] = 'Immune Stromal Desert'
    imm_st_des= cluster_dict[6][['OS.time','E','group']]
    imm_st_des = imm_st_des.dropna(how='any')
    print imm_st_des
    
    survival_times = pd.concat([imm_rich,imm_st_rich,tcell,myeloid,imm_des,imm_st_des],axis=0)
    survival_times.to_csv('TCGA_'+indication+'_clusters_survival_PFI_logCPM_May4_2020.tsv',sep='\t',header=True)
  #survival_times = pd.concat([high,low,low_hi,hi_low],axis=0)
  
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
  
    # kmf.fit(T[ih], E[ih], label='High '+signature+' Signature',alpha=1)
    # ax= kmf.plot()
    # kmf.fit(T[ix], E[ix], label='Low '+signature+' Signature',alpha=1)
    # ax= kmf.plot(ax=ax)
  
    kmf.fit(T[ir], E[ir], label='Immune Rich',alpha=1)
    ax= kmf.plot(color='r',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      ir_med = int(kmf.median_)
    except:
      ir_med = np.nan
    kmf.fit(T[isr], E[isr], label='Immune Stromal Rich',alpha=1)
    ax= kmf.plot(ax=ax ,color = 'g',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      isr_med = int(kmf.median_)
    except:
      isr_med = np.nan  
    kmf.fit(T[tc], E[tc], label='Tcell Centric',alpha=1)
    ax= kmf.plot(ax=ax,color = 'b',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      tc_med = int(kmf.median_)
    except:
      tc_med = np.nan
      
    kmf.fit(T[my], E[my], label='Myeloid Centric',alpha=1)
    ax= kmf.plot(ax=ax,color='grey',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      mc_med = int(kmf.median_)
    except:
      mc_med = np.nan
      
    kmf.fit(T[imd], E[imd], label='Immune Desert',alpha=1)
    ax= kmf.plot(ax=ax,color ='yellow',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      id_med = int(kmf.median_)
    except:
      id_med = np.nan
  
    kmf.fit(T[isd], E[isd], label='Immune Stromal Desert',alpha=1)
    ax= kmf.plot(ax=ax,color='orange',show_censors=True,censor_styles={'ms':6,'marker':'|'})
    try:
      isd_med = int(kmf.median_)
    except:
      isd_med = np.nan
    #results = logrank_test(high_['OS.time'], low_samples['OS.time'], event_observed_A=high_samples['E'],event_observed_B=low_samples['E'], alpha=.99)
    #results = multivariate_logrank_test(T,groups,E, alpha=.99)
    results = logrank_test(T[~ir],T[ir], E[~ir],E[ir], alpha=.99)
    results.print_summary()
    print results.p_value
    #ax= kmf.plot(ax=ax)
    plt.title('5 Year Survival Curve for TCGA '+indication+' Immune Clusters',fontweight='bold')
    plt.text(1000,1, 'P-value: %.2E'%results.p_value, horizontalalignment='center',verticalalignment='center')#,fontweight='bold')
    plt.ylim(0,1.19)
    plt.xlim(0,2000)
    #
    #plt.legend(loc='upper right',fontsize=10)
    plt.savefig('TCGA_'+indication+'_immune_clusters_5yr_overall_survival_April2_23.png',dpi=500)
    plt.savefig('TCGA_'+indication+'_immune_clusters_5yr_overall_survival_April2_23.svg',dpi=500,format='svg')
    plt.savefig('TCGA_'+indication+'_immune_clusters_5yr_overall_survival_April2_23.pdf',dpi=500,format='pdf')
    # out = open('TCGA'+indication+'_immune_clusters_10yr_overall_survival_Dec10_medians.tsv','w')
    medians= {'Immune Rich':ir_med,'Tcell Centric':tc_med,'Immune Desert':id_med,'Immune Stromal Rich':isr_med,'Myeloid Centric':mc_med,'Immune Stromal Desert':isd_med}
    tot_medians[indication] = medians
    print tot_medians  
    median_df = pd.DataFrame.from_dict(tot_medians,orient='columns')
    print median_df
    median_df.to_csv('TCGA_'+indication+'_immune_clusters_5yr_overall_survival_April3_medians.tsv',sep='\t',header=True)
  
  
             
      #plt.show()

if __name__ == '__main__':
     main() 