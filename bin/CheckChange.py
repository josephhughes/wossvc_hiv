#!/usr/bin/env python3.6

import sys
import os
import argparse
from collections import defaultdict 
import re
#import pandas as pd 
#import seaborn as sns
#import numpy as np
#import datetime


# examle:
# python3 CheckChange.py  -r registry.txt -i run_10-11-19.txt -o registry_10-11-2019.txt
# python3 CheckChange.py  -r registry_10-11-2019.txt -i run_12-11-19.txt -o registry_12-11-2019.txt
# python3 ../CheckChange.py  -r samplesUpTo2009.fa.results.txt -i samplesUpTo2010.fa.results.txt -o registry_2010.txt
# TO DO: summary of cluster growth overtime


parser = argparse.ArgumentParser()
parser.add_argument('--registry',"-r", help='file containing the registry')
parser.add_argument('--input',"-i", help='the most recent results of a clustering run csv file with column labels.')
parser.add_argument('--output',"-o", help='name for the new registry, ideally add the date dd-mm-yyyy.')
parser.add_argument('--log',"-l", help='log output with the key changes')
parser.add_argument('--date',"-d", help='date of the new samples added to the clustering dd-mm-yyyy')

args = parser.parse_args()


if(args.registry=="" or args.input==""):
	print("Error, registry and input must be specified with -r and -i options.")
	exit()

newreg = open(args.output,"w")    
## Read input files
# the registry file that can contain multiple columns
# get what the last cluster number was
regData = {}
regDict = defaultdict(list) 
colnames = []
data = []
lastclusterid=0
with open(args.registry) as f:
    header = f.readline()
    header=header.rstrip("\n\r")
    colnames=header.split("\t")
#    print(header)
    for line in f:
      if line!="" and line!="\n":
        line=line.rstrip("\n\r")
        linelist=line.split("\t")
        records=len(linelist)-1
        regDict[linelist[-1]].append(linelist[0]) 
        cid=re.split(r'(\d+)',linelist[-1])
        if int(cid[1])>lastclusterid:
          lastclusterid=int(cid[1])
        regData[linelist[0]]=line
        #print(">",linelist[0],"<")
f.close()


# new file of cluster ids
inDict=defaultdict(list)
with open(args.input) as f:
    header = f.readline()
    header=header.rstrip("\n\r")
    names=header.split("\t")
#    now = datetime.datetime.now()
#    time_now = (str('{:02d}'.format(now.day)) + "-" +str('{:02d}'.format(now.month)) + "-" + str(now.year))

    #print(header)
    colnames.append(args.date)
    for line in f:
      if line!="" and line!="\n":
        line=line.rstrip("\n\r")
        linelist=line.split("\t")
        inDict[linelist[1]].append(linelist[0]) 
f.close()

#print(regDict) 
#print(inDict)

# new cluster numbers and for remapping cluster names
cluster=defaultdict(list)
# dictionary of growing clusters
addition=defaultdict(list)
# dictionary of new clusters
newcluster=defaultdict(list)
# dict of the seqid that have been seen
seen=defaultdict(list)
newreg.write('\t'.join(map(str, colnames)) + "\n")
for regkey in regDict.keys():
  for inkey in inDict.keys():
    #print(regDict[regkey],regkey,"<=>",inkey,inDict[inkey])
    #print(list(set(inDict[inkey])-set(regDict[regkey])))
    difs=list(set(inDict[inkey])-set(regDict[regkey]))
    dif2=list(set(regDict[regkey])-set(inDict[inkey]))
    if len(difs)==0 and len(dif2)==0:
#      print("Cluster ",inkey," from input is ",regkey," from registry")
      for i in regDict[regkey]: 
        # add a statement for cases where the id might already have been seen
#        print(i,"\t",regkey) 
        newline=regData[i]+"\t"+regkey
        vals = re.split(r'\t+', newline)
        if vals[0] not in seen:
          data.append(vals)
          seen[vals[0]]=regkey
          newreg.write(regData[i]+"\t"+regkey+"\n")
          #print("SEEN >"+vals[0]+ "< " +regkey)
        
        
    # if all inDict seqIDs are in regDict, then it is an addition of seqIDs to a cluster
    # it also needs to be added to the largest pre-existing cluster
    #elif len(dif2)==0 and len(regDict[regkey])>=len(difs):
    elif len(dif2)==0:
#      print("Looks like a new addition of ",difs,"to cluster ",regkey," in registry")
      for i in regDict[regkey]: 
        # add a statement for cases where the id might already have been seen
        
        # if the sequences being added to a cluster already had a cluster id, then 
        # need to make sure that they also get detected
        # need to check that subsumed sequences into a cluster don't get called again
        
    #    print(regData[i],"\t",regkey) 
        newline=regData[i]+"\t"+regkey
   #     print(newline)
        vals = re.split(r'\t+', newline)
        if vals[0] not in seen:
          newreg.write(regData[i]+"\t"+regkey+"\n") 
          data.append(vals)
          seen[vals[0]]=regkey
          #print("SEEN >"+vals[0]+ "< " +regkey)
          
        
      for n in difs:
  #      print(n)
        # add a statement for cases where the id might already have been seen
        if n in regData.keys():
   #       print(n)
          newline=regData[n]+"\t"+regkey
          vals = re.split(r'\t+', newline)
   #       print(regData[n]+ " " +regkey)
          if vals[0] not in seen:
            newreg.write(regData[n]+"\t"+regkey+"\n") 
            data.append(vals)
            #print("SEEN >"+vals[0]+ "< " +regkey)
            seen[vals[0]]=regkey

        else:# if sequence is new to this run
 #         print(n,"n-value\n") # number of time to print depends on number of records that needs to be worked out
          newline=n+"\t"+ "NaN\t" * records+regkey
          #seen[n]=regkey
          
  #        print(newline)
          vals = re.split(r'\t+', newline)
          if vals[0] not in seen:
            
            newreg.write(n+"\t"+"NaN\t" * records+regkey+"\n") 
            data.append(vals)   
            #print("SEEN >"+vals[0]+ "< " +regkey)
            seen[vals[0]]=regkey
            addition[regkey].append(vals[0])  
 

# now there are a certain number of new seqIDs in new clusters and we need to start
# the clusterid from lastclusterid
#print("+++++++++++++++++")
#print("Incrementing from the last cluster id in the registry which is :",lastclusterid)

for inkey in inDict.keys():
  #print(inDict[inkey])
  for seqid in inDict[inkey]:
    if seqid not in seen.keys():
      #print(seqid)
  #    print("Not seen "+seqid+ " " + inkey)
      if inkey in cluster.keys():
        newline=seqid+"\t"+ "NaN\t" * records+cluster[inkey]
        vals = re.split(r'\t+', newline)
        if vals[0] not in seen:
          newreg.write(seqid+"\t"+ "NaN\t" * records+cluster[inkey]+"\n")
   #       print(newline)
          data.append(vals) 
          seen[vals[0]]=cluster[inkey]
          newcluster[cluster[inkey]].append(vals[0])
      else:
        lastclusterid+=1
        clustername="c"+str(lastclusterid)
        cluster[inkey]=clustername
        newline=seqid+"\t"+ "NaN\t" * records+cluster[inkey]
        vals = re.split(r'\t+', newline)
        if vals[0] not in seen:
          newreg.write(seqid+"\t"+ "NaN\t" * records+cluster[inkey]+"\n")
   #       print(newline)
          data.append(vals)
          seen[vals[0]]=cluster[inkey] 
          newcluster[cluster[inkey]].append(vals[0])
      #if 
      #newline=seqid+"\t"+ "NaN\t" * records+regkey
newreg.close()             

# little summary
logsummary = open(args.log,"w")    
for clstr in addition.keys():
    #print("The cluster ",clstr," has grown by ",(len(addition[clstr]))," with the addition of ",addition[clstr])
    logsummary.write("The cluster "+str(clstr)+" has grown by "+str(len(addition[clstr]))+" with the addition of "+str(addition[clstr])+"\n")
  
    #print("There are ",len(newcluster.keys())," new clusters with the following sizes")
    logsummary.write("The cluster "+str(clstr)+" has grown by "+str(len(addition[clstr]))+" with the addition of "+str(addition[clstr])+"\n")
for clstr in newcluster.keys():
  #print(clstr," with ",len(newcluster[clstr])," members = ",newcluster[clstr])
  logsummary.write("New clusters"+str(clstr)+" with "+str(len(newcluster[clstr]))+" members = "+str(newcluster[clstr])+"\n")
#print(data)

#now = datetime.datetime.now()
#time_now = (str('{:02d}'.format(now.day)) + "-" +str('{:02d}'.format(now.month)) + "-" + str(now.year))

# df = pd.DataFrame(data, columns = colnames) 
# df = df.replace('NaN',np.NaN)
# df.to_csv(time_now +"_registryUpdate.txt", sep='\t',index=False)
# df=df.melt(id_vars=['SequenceID'])
# 
# df=df.dropna(how='any')
# 
# df=df.groupby(['variable','value']).count().reset_index()
# #print(df.columns)
# df.rename(columns={'variable':'date',
#                           'value':'cluster',
#                           'SequenceID':'count'}, 
#                  inplace=True)
# df['date']=pd.to_datetime(df.date, errors='coerce', format='%d%m%Y')
# count_out=args.output+"_count.txt"

#df.to_csv(count_out, sep='\t',index=False)
#print(df)
  
# sns_plot = sns.lineplot(x="date", y="count", hue="cluster", data=df)
# sns_plot.figure.savefig("growth.png")
