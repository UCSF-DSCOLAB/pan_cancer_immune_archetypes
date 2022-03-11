import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns; sns.set()
import sys

def get_indication(df):
  indication =[]
  for sample in df.index:
    name = sample[sample.rfind('IPI')+3:sample.find('.')-3]
    if name == 'MELB':
      name = 'MEL'
    indication.append(name.strip())
  return indication

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

def rho(x, y):
    return (stats.spearmanr(x, y)[0], stats.spearmanr(x, y)[1]) 
  
def plot_lin_reg(x,y,df,flow):
  sns.set_style("whitegrid")
  #print df.shape
  x = np.concatenate(x.values,axis=0)
  y = np.concatenate(y.values,axis=0)
  slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
  print("slope: %f    intercept: %f" % (slope, intercept))
  print("r-squared: %f" % r_value**2)
  print("std err: %f" % std_err)
  print("P Value: %f" % p_value)
  rho_corr , rho_p_value = rho(x,y)
  print("Spearman Rho: %f" % rho_corr)
  print("Spearman P value: %f" % rho_p_value)
  
  indications = get_indication(df)
  df['Indication'] = indications
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
        'GBM': 'Grey'
        }
  c1 = pd.Series(indications).map(colors)

  df['size'] = 20
  ax =sns.regplot(x='Score', y=flow, data = df ,color='black')#,stat_func=r2),kind="reg",
  ax = sns.scatterplot(x='Score', y=flow, data = df,hue='Indication',palette=colors, legend=False,s=300, linewidths= 25,edgecolor='black')

  plt.xlim((-5,105))
  plt.ylim((-5,100))
  
  plt.text(20,80, 'Spearman Rho: %.2f'%rho_corr, horizontalalignment='center',verticalalignment='center',fontsize= 13,fontweight='bold')
  plt.xlabel('Score',fontsize= 12, fontweight='bold')
  plt.ylabel('flow',fontsize= 8, fontweight='bold')

  plt.savefig('flow_vs_score.svg',format='svg',dpi=500)
  plt.savefig('flow_vs_score.pdf',format='pdf',dpi=500)

  plt.clf()
  return (slope,intercept)


  
def main():
  #use flow populations and matching score to make the model

  corr_file = sys.argv[1]
  flow = 'flow_Tcell'
  df= pd.read_csv(corr_file,sep='\t',index_col=[0],header=0,skipinitialspace=True)
  
  df_known = df[[flow,'Score']]
  df_known = df_known.dropna(how='any')
  X = df_known.filter(['Score'])
  Y = df_known.filter([flow])

  slope ,intercept = plot_lin_reg(X,Y,df_known,flow)


if __name__ == "__main__":
    main()