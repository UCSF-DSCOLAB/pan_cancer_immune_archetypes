import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import sqrt, log
from scipy.stats import zscore
from scipy.stats import percentileofscore
from collections import Counter


def get_indication(df):
  indication =[]
  for sample in df.index:
    name = sample[sample.rfind('IPI')+3:sample.find('.')-3]
    if name == 'MELB':
      name = 'MEL'
    indication.append(name.strip())
  return indication

def make_plot(s,title,MIN,MAX,comp):
    '''
    input:  series of flow scores
            Min Max for y axis
    output: boxplot of flow scores by indication sorted for highest median to lowest
    '''
    
    df = s.to_frame(name='Population_Percent')
    sns.set_style("whitegrid")
    f, ax = plt.subplots(figsize=(20,10))
    indication = get_indication(df) 
    print (len(indication))
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
    df['Indication'] = indication

    my_order = df.groupby(by=["Indication"]).median().iloc[::-1]
    my_order = my_order.sort_values(by=["Population_Percent"],ascending=False).index
    print (my_order.tolist())
   
    ax = sns.boxplot(x="Indication", y="Population_Percent", data=df, whis=np.inf,palette=colors,order=my_order,saturation=0.45)
    ax = sns.swarmplot(x='Indication',y='Population_Percent',data=df, size =10, palette=colors,order=my_order,edgecolor='gray')

    plt.ylabel('Score',fontsize= 30, fontweight='bold')
    plt.xlabel('')
    plt.ylim((MIN,MAX))
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=15, fontweight='bold',rotation=0)
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=25, fontweight='bold',rotation=0)
    plt.title(comp+' score from Flow Cytometry' ,fontsize=20,fontweight='bold')
    #plt.show()
    plt.savefig(title+'.pdf',format='pdf',dpi=500)
    plt.savefig(title+'.svg',format='svg',dpi=500)
    plt.close()
    return my_order 

def main():
  
  title="Flow_score_Sept10_2020"
 
  for comp in ['Tcell','Myeloid','Stroma']:
    df = pd.read_csv('../files_used_for_plots/feature_flow_population_percents/'+comp+'_flow_population_percents.tsv',sep='\t',index_col=[0],header=0,skipinitialspace=True)
    df = df.dropna(how='any')
    df =pd.to_numeric(df['flow_'+comp], errors='coerce')
   
    #convert population percent to percentile
    perc_df = df.apply(lambda x: percentileofscore(df.values,x))
   
    perc_df.to_csv(title+'_'+comp+'_percentile_of_percent.tsv',sep='\t',header=True)
    make_plot(perc_df,comp+'_'+title+'_percentile_of_percent',-10,110,comp)
  print ("Done")
  

if __name__ == '__main__':
     main() 