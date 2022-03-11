import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import sys
from math import sqrt, log
from scipy.stats import zscore
from scipy.stats import percentileofscore
import optparse

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def option_parser():
  '''
  Required flags for getting a gene signature score
  '''
  description = 'Get gene signature score'
  usage = "usage: %prog -o <output file name> -g <gene signature file name> -f <gene expresson tsv>  [options] ALL OPTIONS ARE REQUIRED"
  parser = optparse.OptionParser(usage=usage, description=description)
  parser.add_option('-t', help='title used for all output files', dest='t',
                    action='store', metavar="<title>")
  parser.add_option('-g', help='path to a file with a list of genes (HUGO Names) in gene signature.', dest='g',
                    action='store', metavar="<gene sig>")
  parser.add_option('-f', help='path to normalized gene expression file (TPM or logCPM) with samples in the columns and genes in the rows', dest='f',
                    action='store', metavar="<gene file>")
  
  return parser


def initialize():
  # Get input arguments
  parser = option_parser()
  (opt, args) = parser.parse_args()
  title = opt.t
  gene_sig = opt.g
  gene_file =opt.f
  

  if not (title and gene_sig and gene_file):
    print ('\n')
    parser.print_help()
    print ('\n')
    sys.exit('Error: Not all required inputs have been provided')
  
  return(title,gene_sig,gene_file)


def get_indication(df):
  indication =[]
  for sample in df.index:
    name = sample[sample.rfind('IPI')+3:sample.find('.')-3]
    if name == 'MELB':
      name = 'MEL'
    indication.append(name.strip())
  return indication

def get_sample_type(df):
  sample_type=[]
  for sample in df.index:
    name = sample[sample.find('.')+1:sample.find('rna')-2]
    sample_type.append(name.strip())
  return sample_type

    
def make_plot(s,title,MIN,MAX):
  '''
  Makes a box plot of scores by indication
  '''
  df = s.to_frame(name='Score')
  df.to_csv(title+'.tsv',sep='\t', header=True)
  sns.set_style("whitegrid")
  f, ax = plt.subplots(figsize=(20,10))
  sample_type = get_sample_type(df)
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
      'GBM': 'Thistle',
  }  
  df['Indication'] = indication
  df['sample_type'] = sample_type
  
  #print (df)
  df_T = df[df.sample_type == 'T']
  df_N = df[df.sample_type == 'N']
  #print (df_T.shape)
  #print (df_N.shape)
  my_order = df.groupby(by=["Indication"]).median().iloc[::-1]
  my_order = my_order.sort_values(by=["Score"],ascending=False).index

  ax = sns.boxplot(x="Indication", y="Score", data=df, whis=np.inf,palette=colors,order=my_order,saturation=0.45)
  ax = sns.swarmplot(x='Indication',y='Score',data=df_T, size =10, palette=colors,order=my_order,edgecolor='gray')
 
  plt.ylabel('Score',fontsize= 15, fontweight='bold')
  plt.xlabel('')
  plt.ylim((MIN,MAX))
  plt.setp(ax.yaxis.get_majorticklabels(), fontsize=15, fontweight='bold',rotation=0)
  plt.setp(ax.xaxis.get_majorticklabels(), fontsize=15, fontweight='bold',rotation=0)
  plt.title(title ,fontsize=20,fontweight='bold')
  #plt.show()
  plt.savefig(title+'.pdf',format='pdf',dpi=500)
  plt.savefig(title+'.svg',format='svg',dpi=500)
  plt.savefig(title+'.png',dpi=500)


def convert_array_to_percentiles(df):
  percentile_array = []
  for column in df:
    arr = df[column]
    percentiles = arr.apply(lambda x: percentileofscore(arr,x))
    percentile_array.append(percentiles)
  return pd.concat(percentile_array, axis=1)#join_axes=[percentile_array[0].index])


def get_sig_genes(f_name):
  '''
  Read the gene file and return gene_list 
  '''
  return [line.rstrip() for line in open(f_name)]



def get_score(df,genes):
  '''
  Returns a series of the mean of percentiles of the gene signature
  per sample
  '''
  up = df.filter(genes,axis=0)
  return up.T.mean(axis=1)

def drop_all_zero_samples(df):
  '''
  Remove samples that are all zero counts
  '''

  drop_list =[]
  for index, row in df.iterrows():
    if len(set(row.values)) ==1:
      drop_list.append(index)
  print ( drop_list ) 
  return df.drop(drop_list,axis="index")



def main():
  
  title,gene_sig, gene_file = initialize()
  #get a list of gene signature genes
  sig_genes = get_sig_genes(gene_sig)
  
  #filter dataframe
  df_all = pd.read_csv(gene_file, sep='\t',index_col=[0], header=0)
  df_all = df_all[~df_all.index.duplicated()]
  df = df_all.filter(sig_genes,axis=0)
  print (df)
  df = drop_all_zero_samples(df)
  
  
  perc_df = convert_array_to_percentiles(df.T)
  print (perc_df.shape)
  
  mean_s = get_score(perc_df.T,sig_genes)
  
  mean_perc_df = mean_s.apply(lambda x: percentileofscore(mean_s,x))
  
  make_plot(mean_perc_df,title+'_fmcs_EHK8-10_percentile_of_percentiles',-10,110)

  print ("Done")
if __name__ == '__main__':
     main() 