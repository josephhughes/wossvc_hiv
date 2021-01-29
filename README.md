# wossvc_hiv

Pipeline for detecting clusters of cases in HIV data from WoSSVC

## Install
```
git clone https://github.com/josephhughes/wossvc_hiv.git
cd wossvc_hiv
conda env create -f environment.yml
conda activate wossvc_hiv
```


## The pipeline 

### Processing the sequences from a run

1. Query sierrapy to obtain the subtype assignment and drug resistance information.

2. Add the sequences to an alignment with a reference set of sequences using mafft.

3. Reconstruct the phylogeny using raxml

4. Plot the phylogeny using ete3.

5. Provide a run report word document.

### Processing the subtype C sequences through hivtrace [optional]

1. Filter subtype C sequences

2. Run hivtrace

3. Obtain cluster assignment from hivnetworkcsv

4. Compare the changes in cluster assignment to the previous cluster assignment

5. Add the clustering information to the report.



## Running the pipeline
```
run_HIV_analysis_pipeline.sh -- program to assign subtypes, detect DRMs and conduct HIVtrace clustering:\n
        -h      show this help text
        -f      input sequences in fasta format
        -t      date for the output folder dd-mm-yyyy
        -d      [OPTIONAL] if wanting hivtrace clustering then path to database of fasta sequences previously processed by HIVtrace
        -r      [OPTIONAL] only necessary if -d is provided, then path to the registry for the database of sequences already provided"
```

### Example:

```
./scripts/run_HIV_analysis_pipeline.sh -f test_data/2018_seqs.fa -r test_data/Pre2018_registry.txt -d test_data/Pre2018_seqs.fa -t 01-01-2019
```

2018_seqs.fa is a set of sequences that need to be processed. Here these are a set of sequences from 2018.

Pre2018_registry.txt is a text-tab delimited file of previously clustered sequences. 
This file represents samples that are from before 2018 and have been assigned as subtype C.

Pre2018_seqs.fa is the set of fasta files that are in the registry.

01-01-2019 is the date for the output files and date added to the registry.

### Follow on example:

```
./scripts/run_HIV_analysis_pipeline.sh -f test_data/2019_seqs.fa -r Results_01-01-2019/01-01-2019.registry.txt -d Results_01-01-2019/01-01-2019_sequence_db.fa -t 01-01-2020
```

Here, we are using as an input the registry and compiled sequences from the previous run 
to increment the sequence file and registry.

## Outputs

### Results_dd-mm-yyyy

`dd-mm-yyyy_RunOverview.docx`
> This is a report for the run which contains the number of sequences assigned to each subtype, the phylogeny for the sequences in the run and the results of the clustering

`dd-mm-yyyy_DRM-overview.txt`  
> Text-tab delimited file with the resistance associated mutations detected by sierrapy.

`dd-mm-yyyy.log`
> Summary results of the clustering

`dd-mm-yyyy.registry.txt`
> The previous and new registry depending on date (text-tab delimited) - IMPORTANT

`dd-mm-yyyy_sequence_db.fa`
> The previous and new dataset of subtype C sequences - IMPORTANT

### dd-mm-yyyy_reports

> Output of sierrapy formatted for word

### hivtrace_dd-mm-yyy

> Outputs from hivtrace

### RAxML

> Output from the phylogenetic analysis

## Dependencies

wossvc_hiv uses the following tools:

```
mafft
raxml
conda
biopython
python3.6
python-docx
tn93
pyqt
ete3 (https://github.com/etetoolkit/ete)
sierrapy (https://github.com/hivdb/sierra-client/tree/master/python)
hivtrace (https://github.com/veg/hivtrace)
```