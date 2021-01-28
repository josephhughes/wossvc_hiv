#!/usr/bin/env python3.6
'''
Short Python script to parse json returned from sierrapy query and store sample metadata in tab delimited text file.
Fields include sampleID, month, actualMonth and year plus a column for each drug and corresponding resistance mutation (if any).
'''

import json
import argparse
import datetime


parser = argparse.ArgumentParser()
parser.add_argument('--json', required=True, help='input json file containing query results')
parser.add_argument('--output', required=False, help='name of tab-delimited text file containing sample metadata')
args = parser.parse_args()

# get user input and set file names
json_in = ""
output_file=""

if args.json:
    json_in = args.json
else:
    print ("json file cannot be read by parser")

now = datetime.datetime.now()

time_now = (str('{:02d}'.format(now.day)) + "-" +str('{:02d}'.format(now.month)) + "-" + str(now.year))

if args.output:
    output_file = args.output
else:
    output_file = time_now + "_DRM-overview.txt"


# parse sierrapy json output for relevant info
header_list = ['SampleID', 'ATV', 'DRV', 'FPV', 'IDV' , 'LPV', 'NFV', 'SQV', 'TPV', 'ABC', 'AZT', 'D4T', 'DDI', 'FTC', 'LMV', 'TDF', 'DOR', 'EFV', 'ETR', 'NVP', 'RPV']

append_dict={}

with open(json_in) as json_file:
    data = json.load(json_file)
    for i in data:
        ID = i["inputSequence"]["header"]
        append_dict[ID] = {}
        append_dict[ID]['SampleID']=ID
        for j in i["drugResistance"]:
            for k in j["drugScores"]:
                drug_name = k["drug"]["name"]
                drug_score = k["score"]
                if drug_score != 0.0:
                    mutation_list=[]
                    for m in k["partialScores"]:
                        for key,value in m.items():
                            if key =='mutations':
                                for x, y in value[0].items():
                                    if x == 'text':
                                        mutation_list.append(y)
                                        
                    append_dict[ID][drug_name] = mutation_list
                else:
                    append_dict[ID][drug_name] = 0

f = open(output_file, "w")
f.write('\t'.join(header_list) + '\n')
for ids in append_dict:
  f.write(ids)
  for cols in header_list[1:]:
    f.write("\t"+str(append_dict[ids][cols]))
  f.write("\n")
f.close()