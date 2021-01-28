#!/usr/bin/env python3.6

import argparse
import json
import os
import datetime
import re

parser = argparse.ArgumentParser()
parser.add_argument('--json', required=True, help='input json file containing HIV-TRACE cluster output')
parser.add_argument('--out', required=True, help='output text-tab for seq ids and cluster')

args = parser.parse_args()

json_in = ""


# get user input - json returned from HIV-TRACE
if args.json:
    json_in = args.json
else:
    print ("json file cannot be read by parser")

# generate name for .txt output file
#base = os.path.basename(json_in)
#output = os.path.splitext(base)[0] + ".txt"
output = args.out

# parse the json for cluster information, store in a dict as cluster to list of samples
hivtrace_clusters = {}
with open(json_in) as json_file:
    data = json.load(json_file)


traceResults = data['trace_results']

for counter, cluster in enumerate(traceResults['Nodes']['cluster']['values']):
    print(counter,"\t",cluster,"\t",traceResults['Nodes']['id'][counter])
    seq_id=traceResults['Nodes']['id'][counter]
    
    result = re.sub(r"\|DUPLICATE \d+","", seq_id)
    #print(result)
    seq_id=result
    if cluster in hivtrace_clusters:
        hivtrace_clusters[cluster].append(seq_id)
    else:
        hivtrace_clusters[cluster] = [seq_id]


now = datetime.datetime.now()
time_now = (str('{:02d}'.format(now.day)) + "-" +str('{:02d}'.format(now.month)) + "-" + str(now.year))
# write to file. for each cluster, loop through the list of sequence IDs and write to file
with open(output, 'w+') as out:
    out.write("SeqID\t" + time_now + "\n")
    for key, value in hivtrace_clusters.items():
        for i in value:
            out.write(i + "\tc" + str(key) + "\n")


