SMC-RNA Challenge Entry: uhrigs-8396553
========================================================

This is a CWL workflow which was submitted to the [SMC-RNA Challenge][smcrna] in the fusion detection category.

Example and reference data files are stored on [Synapse][data] and are
defined in synapse_inputs.json.  
A Synapse user account is required to access data from Synapse.

Usage
--------------------------------------------------------

```bash
# Clone this repo
git clone https://github.com/smc-rna-challenge/uhrigs-8396553
cd uhrigs-8396553

# Install the requirements
pip install -r requirements.txt

# Download the test data to ./data
# You'll need a synapse.org account.
python download.py

# Run the workflow
cwltool main.cwl
```

Resource Usage
--------------------------------------------------------

Most SMC-RNA entries were tested on machines with 16 vCPUs and 60 GB of RAM.
Many entries did not use the full amount of resources.
Some required more memory and were allowed 104 GB of RAM.


Docker images
--------------------------------------------------------


- [quay.io/smc-rna-challenge/uhrigs-8396553-smc-rna-challenge-3:1.0.0](https://quay.io/smc-rna-challenge/uhrigs-8396553-smc-rna-challenge-3:1.0.0)





Details
========================================================

Challenge participants were asked to complete a quiz to describe some of the
details about their entry:



#### What methods does the tool use?
Split-read mapping



#### What alignment algorithms does the tool use?
STAR



#### Are there non-default parameters used for the alignment algorithms?
Yes



#### Are there filters used to generate the final results?
Blacklists



#### Are there other inputs other than sequencing reads and gene annotation?
Yes



#### Does the tool report the following types of fusions?
Inter-chromosomal, Intra-chromosomal, read-through, Self-fusion



#### Does the tool report fusions with their junctions in the following regions?
Intergenic



#### Does the tool report only the consensus isoform of a gene fusion, or all the isoforms?
All the isoforms



#### Does the tool report any scores for gene fusions?
Other



#### Does the tool have a cutoff of number of reads to call fusions?
Yes



#### Does your method require isoform models to be supplied or is it reference free, able to construct transcripts de novo?
Isoform Models must be supplied



#### If your method detects fusions, does it incorporate information about the expression levels of the 5’ or 3’ partner genes?
Yes



#### Other than transcript models, does your method make use of any prior information not available in the data? If so, does it use:
highly recurrent fusions from CancerGeneCensus, blacklist



#### Is your approach a stand-alone algorithm or would you consider it to be an ensemble of multiple methods that are later combined in some way?
Stand-alone method



#### Does your method require any specific parameters for tuning (e.g. mean and standard deviation of expected fragment sizes)?
No



#### Do you anticipate your algorithm would need to be parameterized for different read-lengths?
No



#### If your method detects cancer-associated fusions, does it account for potential artifacts? If so, which ones (check any that apply):
immunoglobulin (IG) rearrangements in tumor samples infiltrated by immune cells, nuclear mitochondrial DNA, read-through chimeras produced by co-transcription of adjacent genes, recurrent alignment artifacts



#### If your method detects cancer-associated fusions, does it filter out those found in normal (or other control) samples?
Yes



#### If your method predicts gene fusions, do you use canonical exon boundaries from a database?
Yes



#### If your method predicts gene fusions, do you output isoforms that could contain non-canonical exon boundaries?
Yes



#### Has your method been designed to address any of the following explicitly (check any that apply):
Over-represented k-mers, Duplicate reads, Large amounts of ribosomal RNA, Sequencing depth of the library



#### If your method infers expression levels, does it simultaneously identify transcripts and their expression levels or does it perform these steps sequentially?
Question not applicable to my method.



#### Is your method alignment independent?
No



#### Does your method use the quality of the read mapping at all (e.g. filter out poorly mapped reads)?
Yes



#### Does your method use information about strand-specific sequences?
Yes



#### Does your method output transcript-length- and library size- adjusted estimates of abundance (such as RPKM, FPKM, or TPM)?
No



#### We are developing a new DREAM Challenge on visualization of transcript isoforms and gene fusions. This visualization challenge is a companion to the currently running 2016 DREAM SMC-RNA (somatic mutation calling-RNA) Challenge that you are participating in. Please find more information here - http://biovis.net/2016/dream/. Do note that the content of this site will continue to evolve. Would you like to share your results from this challenge with the SMC-RNA BioVis Data Visualization DREAM Challenge?
No





[smcrna]: https://www.synapse.org/#!Synapse:syn2813589/wiki/401435
[data]: https://www.synapse.org/#!Synapse:syn9878868
