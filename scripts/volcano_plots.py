#!/usr/bin/env python

import glob
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import optparse
from adjustText import adjust_text
import networkx as nx
'''
Septemeber 2017
makes a volcano plot from a standard Limma DGE output
e.g.
        logFC	    |   AveExpr	  |       t	    |  P.Value	|  adj.P.Val	|      B
PLK1	7.374675473	| 1.415522735	| 6.886987802	| 3.92E-06	| 0.050802876	| 3.971754014
ZWINT	5.772193613	| 1.340166868	| 5.872941797	| 2.48E-05	| 0.065626979	| 2.450364751

Expects Hugo names. There are function to convert from ensembl names (commented out)
Inputs are pvalue cutoff input file name output file name and title of plot


******************************************************************
Things to do
Add

adjust_text(texts, force_points=0.2, force_text=0.2,
            expand_points=(1, 1), expand_text=(1, 1),
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
******************************************************************
'''
def option_parser():
  '''
  Required flags for making a volcano plot
  '''
  description = 'Make Volcano Plots'
  usage = "usage: %prog -i <input file name> [options] "
  parser = optparse.OptionParser(usage=usage, description=description)
  parser.add_option('-i', help='Path and name of input file', dest='i',
                    action='store', metavar="<input file>")
  parser.add_option('-o', help='Path and names of output files e.g DGE_treg_vs_tcell (a .png,pdf &svg file extension will automatically be added)', dest='o',
                    action='store', metavar="<output file>")
  parser.add_option('-p', help='P-value cutoff DEFAULT 0.005', dest='p',
                    default= 0.005,action='store', metavar="<p-value>")
  parser.add_option('-t', help='Title in quotes e.g. "Tregs vs Tcells" ', dest='t',
                    action='store', metavar="<Title>")
    
  return parser

    
def make_volcano_plot(df,pvalue,title,outfile):


  plt.figure(figsize=(10,8))
  plt.scatter(df['logFC'],df['log10Pvalue'],color='b')
  
  # x and y values that should be highlighted
  x =  df['logFC'][abs(df['P.Value']) < pvalue]
  y =  df['log10Pvalue'][abs(df['P.Value']) < pvalue]
  labels = df['log10Pvalue'][abs(df['P.Value']) < pvalue].index
  
  plt.scatter(x,y,color='r',alpha=0.7,s=75,label='P-Value . 0.005')
  for label, x1, y1 in zip(labels, x, y):
    plt.annotate(
        label,
        xy=(x1, y1), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.001', fc='white', alpha=0.3),fontsize=8, fontweight='bold',
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    
  plt.xlim((-7,5))
  plt.ylim((-1,7))
  plt.xlabel('logFC')
  plt.ylabel('-log10 P-Value')
  plt.title(title) 
  #plt.show()
  plt.savefig(outfile+'.png',dpi=500)
  plt.savefig(outfile+'.svg',dpi=500)
  plt.savefig(outfile+'.pdf',dpi=500)
  return        

def initialize():
  # Get input arguments
  parser = option_parser()
  (opt, args) = parser.parse_args()
  infile = opt.i
  outfile = opt.o
  pvalue = float(opt.p)
  title = opt.t

  if not (infile and outfile and title):
    print ('\n')
    parser.print_help()
    print ('\n')
    sys.exit('Error: Not all required inputs have been provided')
    # print '\n'        
    # return
 
  return(infile, outfile,pvalue,title)

def main():
  infile,outfile,pvalue,title = initialize() 
  df = pd.DataFrame.from_csv(infile,sep='\t',header=0)
  df['log10Pvalue'] = -(np.log10(df['P value of: E vs. A'].values))
  df = df.filter(items=['log E vs. A','log10Pvalue','P value of: E vs. A'])
  df = df.apply(pd.to_numeric)
  make_volcano_plot(df,pvalue,title,outfile)

  

if __name__ == "__main__":
    main()             