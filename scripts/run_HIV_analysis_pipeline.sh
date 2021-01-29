#!/bin/bash

script_dir=$(pwd)
#script_dir='/home2/db/HIV-cluster'

# Provide as input:
# -f the set of reference sequences
# -c option for if you want cluster analysis performed, if used, then the following is needed
# -d the set of previously processed sequences in hivtrace, i.e, curated set of subtype C
# -r the associated registry, text-tab file with cluster assignment for past runs


FASTA=""
REF="${script_dir}/data/HIV_aligned_references.fasta"
CLUSTER="false"
DB=""
REGISTRY=""

usage="run_HIV_analysis_pipeline.sh -- program to assign subtypes, detect DRMs and conduct HIVtrace clustering:\n
        -h      show this help text
        -f      input sequences in fasta format
        -t      date for the output folder dd-mm-yyyy
        -d      [OPTIONAL] if wanting hivtrace clustering then path to database of fasta sequences previously processed by HIVtrace
        -r      [OPTIONAL] only necessary if -d is provided, then path to the registry for the database of sequences already provided"

while getopts ":h:f:t:d:r:" opt;
do
        case $opt in
                h)
                        echo -e "$usage"
                        exit 1
                        ;;
                f)
                        echo -e "\nInput fasta: $OPTARG \n" 
                        FASTA=$OPTARG
                        ;;
#                 c)
#                         echo "Clustering analysis of sequences will be performed"
#                         CLUSTER="true"
#                         #DB="/home1/hugh01j/Documents/HIV_Samantha/HIV_detection_pipeline/test_data/Inregistry.fa"
#                         #REGISTRY=$OPTARG
#                         ;;
                t)
                        echo -e "\nDate to use for the output resuld folder: $OPTARG \n" 
                        datelabel=$OPTARG
                        ;;
                d)
                        echo -e "\nInput database of sequences: $OPTARG \n" 
                        DB=$OPTARG
                        ;;
                r)
                        echo -e "\nInput registry associated with database sequences: $OPTARG \n" 
                        REGISTRY=$OPTARG
                        ;;
               \?)
                        echo "invalid option -$OPTARG"
                        echo $usage
                        ;;
                :)
                        echo "Arguments required"
                        exit 1
                        ;;
        esac
done

if [ $OPTIND -eq 1 ]; then echo "$usage"; fi

# CHECK USER INPUT PARAMETERS!!
if [ -z "$FASTA" ];
then
        echo "Missing parameter - fasta file!"
        exit 0
fi

mac2unix $FASTA
dos2unix $FASTA
wait

# use code from getSubtypeInformation.sh to perform subtype query

echo "Number of sequences: "
grep -c '>' $FASTA
echo ""

# create a results folder and tidy output files
if [ -z "$datelabel" ]
then
  now="$(date +'%d-%m-%Y')"
  echo "No datelabel provided, using $now"
else
  now=$datelabel
fi


RESULTS="Results_${now}"
mkdir ${RESULTS}
cp $FASTA ${RESULTS}/.
cd ${RESULTS}


BASE=$(basename $FASTA)
FILENAME=$(basename "${BASE%.*}")


JSON="$now.$FILENAME.json"
TEXTFILE="$now.$FILENAME.txt"

${script_dir}/bin/perform_query.py --fasta $BASE --json $JSON

wait
echo "Subtype query complete. Parsing json..."


${script_dir}/bin/parse_json_write_docx.py --json $JSON --output ${now}_subtypes.txt --reports
wait

${script_dir}/bin/parse_json_store_metadata.py --json $JSON
wait

echo "Finished generating reports."

OUTPUT="$now.$FILENAME.aln.fasta"
# Run alignment
echo "Running mafft alignment. Will save to :  $OUTPUT"

mafft --add $BASE --reorder $REF > $OUTPUT
wait

# Generate phylogeny
echo "Running phylogeny using Raxml and Jukes Cantor model with bootstrap analysis..."
raxmlHPC -f a -m GTRGAMMA --JC69 -T 8 -p 12345 -x 12345  -# 100 -s $OUTPUT -n $now
TREE="RAxML_bipartitionsBranchLabels.$now"
wait
# Visualise phylogeny as pdf
${script_dir}/bin/visualise_phylogeny.py --tree $TREE --reroot


RAX="RAxML"
mkdir $RAX
mv RAxML* $RAX

TREE_PNG="${RAX}/RAxML_tree-rerooted.png"
cd ..

# CLUSTERING
if [ -z "$DB" ] && [ ]-z "$REGISTRY" ];
then
  # CHECK USER INPUT PARAMETERS!!
  echo "Missing parameter - the database of past sequences and associated registry!"
  echo "Not performing hivtrace"
  # write run report
  echo ${script_dir}/bin/write_run_overview.py --fasta $BASE --subtypes ${now}_subtypes.txt --tree ${TREE_PNG}
  ${script_dir}/bin/write_run_overview.py --fasta $BASE --subtypes ${now}_subtypes.txt --tree ${TREE_PNG}

  exit 0
else
  cp $DB ${RESULTS}/.
  cp $REGISTRY ${RESULTS}/.
  cd ${RESULTS}/
  dbname=$(basename $DB)
  registry=$(basename $REGISTRY)
  echo "Running cluster analysis..."
  echo "Extracting Subtype C sequences..."

  # split input data based on subtype & add the subtype C samples to the main database of sequences
  newdb="sequence_db.fa"
  echo ${script_dir}/bin/filterSubtypeC.py --fasta $BASE --subtypes ${now}_subtypes.txt --db ${dbname} --out ${now}_${newdb}
  ${script_dir}/bin/filterSubtypeC.py --fasta $BASE --subtypes ${now}_subtypes.txt --db ${dbname} --out ${now}_${newdb}
  HIVTRACE_OUTPUT="hivtrace_${now}"
  mkdir ${HIVTRACE_OUTPUT}

  # run HIV-TRACE + parse output
  echo "Calling HIV-TRACE..."
  echo hivtrace -i ${now}_${newdb} -a resolve -r HXB2_pol -t 0.015 -m 500 -g .05 -s wheeler -o ${now}.json
  hivtrace -i ${now}_${newdb} -a resolve -r HXB2_pol -t 0.015 -m 500 -g .05 -s wheeler -o ${now}.json
  mv hivtrace.log ${HIVTRACE_OUTPUT}
  mv ${now}.json ${HIVTRACE_OUTPUT}
  mv ${now}_${newdb}_output.fasta ${HIVTRACE_OUTPUT}
  
  # call script to parse json output
  #echo  ${script_dir}/bin/parse_hivtrace_output.py --json ${HIVTRACE_OUTPUT}/${now}.json --out ${HIVTRACE_OUTPUT}/${now}.txt
  #${script_dir}/bin/parse_hivtrace_output.py --json ${HIVTRACE_OUTPUT}/${now}.json --out ${HIVTRACE_OUTPUT}/${now}.txt
  
  # get cluster assignments for singletons as well
  hivnetworkcsv -i ${now}_${newdb}_user.tn93output.csv  -t 0.015 -f plain -J -q --singletons include -c ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.csv > ${now}_withsingletons.json
  mv ${now}_${newdb}_user.tn93output.csv ${HIVTRACE_OUTPUT}
  #sed 's/,/\tc/' ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.csv > ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.txt  
  sed $'s/,/\tc/' ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.csv > ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.txt  
  
  
  printf "\nReached the end of clustering pipeline.\n\n.Updating registry...\n"
  
  # this provides input to Joseph's script - run this here
  echo ${script_dir}/bin/CheckChange.py --registry ${registry} --input ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.txt --output ${now}.registry.txt --log ${now}.log --date $now
  ${script_dir}/bin/CheckChange.py --registry ${registry} --input ${HIVTRACE_OUTPUT}/${now}_${newdb}_cluster.txt --output ${now}.registry.txt --log ${now}.log --date $now

  # collate output and generate overall report for this run

  
  # add line to write run report - WITH CLUSTER INFO!
  echo ${script_dir}/bin/write_run_overview.py --fasta $BASE --subtypes ${now}_subtypes.txt --tree $TREE_PNG --clusters 01-01-2020.log --reg ${now}.registry.txt
  ${script_dir}/bin/write_run_overview.py --fasta $BASE --subtypes ${now}_subtypes.txt --tree $TREE_PNG --clusters 01-01-2020.log --reg ${now}.registry.txt

  # Generate an ML phylogeny for all the sequences
  echo "Running phylogeny using Raxml and GTR+GAMMA model on all the sequences WITHOUT bootstrap analysis..."
#  raxmlHPC -f a -m GTRGAMMA -T 8 -p 12345 -x 12345 -s ${now}_${newdb} -n ${now}_${newdb}
  wait


  cd ..

fi



