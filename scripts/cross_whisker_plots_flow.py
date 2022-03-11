import random
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

random.seed(21212)

def make_cross_whisker_plot(medians,iqr_tcga,iqr_ipi):
  '''
  Makes a cross-whisker plot that shows the median value, the interquartile range (IQR)
  in both the x and y directions and the identity line  
  '''
  
  colors = { 
        'LUNG':'MidnightBlue',
        'KID':'Red',
        'HNSC':'OliveDrab',
        'CRC':'DarkTurquoise',
        'GYN':'Coral',
        'BLAD':'Gold',
        'MEL':'Maroon',
        'HEP':'DodgerBlue',
        'SRC':'DeepPink',
        'PNET':'DarkOrchid',
        'PDAC' : 'Yellow',
        'GALL':'Green',
        'ADR':'Purple',
        'GSTR':'Yellow',
        'SI' : 'Black',
        'GBM': 'Grey'
        }
  

  indications = medians.keys()
  fig, ax = plt.subplots(1,1, figsize=(10,10))
  legend = []

  for i in indications:
      ax.scatter(*medians[i], color=colors[i], s=250, marker="s")  # The point itself
      # x iqrs
      ax.plot([iqr_tcga[i][0],iqr_tcga[i][1]], [medians[i][1], medians[i][1]], color=colors[i],linestyle="--")
      ax.plot([iqr_tcga[i][0],iqr_tcga[i][0]], [medians[i][1]-1, medians[i][1]+1], color=colors[i],linestyle="--")
      ax.plot([iqr_tcga[i][1],iqr_tcga[i][1]], [medians[i][1]-1, medians[i][1]+1], color=colors[i])
      
      # # y iqrs
      ax.plot([medians[i][0], medians[i][0]],[iqr_ipi[i][0],iqr_ipi[i][1]], color=colors[i],linestyle="--")
      ax.plot([medians[i][0]-1, medians[i][0]+1],[iqr_ipi[i][0],iqr_ipi[i][0]], color=colors[i],linestyle="--")
      ax.plot([medians[i][0]-1, medians[i][0]+1],[iqr_ipi[i][1],iqr_ipi[i][1]], color=colors[i])
              
      legend.append(matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor=colors[i], marker='s', label=i, markersize=10, alpha=1))
  x = np.linspace(0,100,100)
  y = x
  plt.plot(x, y, '-k', label='y=x')
  #ax.legend(handles=legend, loc='lower right')
  ax.set_xlim(0, 100)
  ax.set_ylim(0, 100)
  ax.set_xlabel('Flow',fontweight='bold',fontsize=20)
  ax.set_ylabel('IPI',fontweight='bold',fontsize=20)
  fig.suptitle('Comparison between IPI new alignment and Flow Tcell median scores and IQR',fontweight='bold',fontsize=15)
  plt.savefig('ipi_flow_tcell_new_alignment_comparison_Mar5_2021.pdf',format='pdf', dpi=500)
  plt.savefig('ipi_flow_tcell_new_alignment_comparison_Mar5_2021.svg',format='svg', dpi=500)
  #plt.show()


def get_indication(df):
  indication =[]
  for sample in df.index:
    name = sample[sample.rfind('IPI')+3:sample.find('.')-3]
    if name == 'MELB':
      name = 'MEL'
    indication.append(name.strip())
  return indication

def get_medians(df_t,df_i):
    
  tcga_med = df_t.median(axis=0,numeric_only=True)
  ipi_med = df_i.median(axis=0,numeric_only=True)
  return (tcga_med.values[0],ipi_med.values[0])

def get_IQR(df):
  q1 = df.quantile(0.25,axis=0, numeric_only=True)
  q3 = df.quantile(0.75,axis=0, numeric_only=True)

  return (q1.values[0],q3.values[0])

def get_data(df_ipi,df_flow):
  
  medians ={}
  iqr_flow ={}
  iqr_ipi= {}
  
  ind_ipi = get_indication(df_ipi)
  df_ipi['indication'] =ind_ipi
  #print df_ipi
  ind_flow = get_indication(df_flow)
  df_flow['indication'] =ind_flow
  #print df_tcga
  for ind in ['HEP','HNSC','KID','MEL','BLAD','LUNG','SRC','GYN','CRC','PDAC','GBM','PNET']:
    df_f = df_flow[df_flow.indication==ind]
    df_i = df_ipi[df_ipi.indication==ind]
    medians[ind] = get_medians(df_f,df_i)
    iqr_flow[ind] = get_IQR(df_f)
    iqr_ipi[ind] = get_IQR(df_i)
  return(medians,iqr_flow,iqr_ipi)
    
  
def main():
  df_ipi = pd.read_csv('../files_used_for_plots/tcell_score_percentile_of_percentiles.tsv',sep='\t',index_col=[0],header=0)
  df_flow = pd.read_csv('../files_used_for_plots/Flow_score_Sept10_2020_Tcell_percentile_of_percent.tsv',sep='\t',index_col=[0],header=0)
  medians,iqr_flow,iqr_ipi = get_data(df_ipi,df_flow)
  make_cross_whisker_plot(medians,iqr_flow,iqr_ipi)


if __name__ == '__main__':
     main()  