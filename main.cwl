arguments: []
baseCommand: [bash, run_arriba.sh, '8', '50000', /opt/arriba]
class: CommandLineTool
cwlVersion: sbg:draft-2
dct:creator: {'@id': 'http://orcid.org/0000-0002-7681-6415', 'foaf:mbox': uhrigs@synapse.org,
  'foaf:name': uhrigs}
description: ''
doc: 'SMC-RNA challenge fusion detection submission

  '
hints:
- {class: 'sbg:CPURequirement', value: 8}
- {class: 'sbg:MemRequirement', value: 50000}
- {class: DockerRequirement, dockerImageId: '', dockerPull: 'quay.io/smc-rna-challenge/uhrigs-8396553-smc-rna-challenge-3:1.0.0'}
id: https://cgc-api.sbgenomics.com/v2/apps/uhrigs/smc-rna-challenge-3/arriba-0-8/14/raw/
inputs:
- description: path to tar archive containing STAR index
  id: '#STAR_INDEX_TAR'
  inputBinding: {position: 1, 'sbg:cmdInclude': true, separate: true}
  label: STAR_INDEX_TAR
  sbg:fileTypes: TAR,TAR.GZ
  type: [File]
- description: path to reference annotation in GTF format (may be gzip-compressed)
  id: '#REFERENCE_GTF'
  inputBinding: {position: 2, 'sbg:cmdInclude': true, separate: true}
  label: REFERENCE_GTF
  sbg:fileTypes: GTF.GZ,GTF
  type: [File]
- description: path to reference genome in FastA format (may be gzip-compressed)
  id: '#REFERENCE_GENOME'
  inputBinding: {position: 3, 'sbg:cmdInclude': true, separate: true}
  label: REFERENCE_GENOME
  sbg:fileTypes: FA.GZ,FA
  type: [File]
- description: path to gzip-compressed FastQ file containing first mates
  id: '#TUMOR_FASTQ_1'
  inputBinding: {position: 4, 'sbg:cmdInclude': true, separate: true}
  label: TUMOR_FASTQ_1
  sbg:fileTypes: FQ.GZ,FASTQ.GZ
  type: [File]
- description: path to gzip-compressed FastQ file containing second mates
  id: '#TUMOR_FASTQ_2'
  inputBinding: {position: 5, 'sbg:cmdInclude': true, separate: true}
  label: TUMOR_FASTQ_2
  sbg:fileTypes: FQ.GZ,FASTQ.GZ
  type: [File]
label: arriba-0.8
outputs:
- description: predicted fusions in BEDPE format
  id: '#OUTPUT'
  label: OUTPUT
  outputBinding: {glob: fusions.bedpe}
  sbg:fileTypes: BEDPE
  type: [File]
requirements:
- class: CreateFileRequirement
  fileDef:
  - {fileContent: "#!/bin/bash\n\nset -o pipefail\n\nif [ $# -ne 8 ]; then\n     \
      \   echo \"Usage: $(basename $0) threads memory tools_dir STAR_index.tar.gz\
      \ annotation.gtf.gz assembly.fa.gz read1.fastq.gz read2.fastq.gz\" 1>&2\n  \
      \      exit 1\nfi\n\nset -x\n\n# fetch arguments\nTHREADS=\"$1\"\nMEMORY=\"\
      $2\"\nTOOLS_DIR=\"$3\"\nSTAR_INDEX=\"$4\"\nREFERENCE_GTF=\"$5\"\nREFERENCE_GENOME=\"\
      $6\"\nREAD1=\"$7\"\nREAD2=\"$8\"\nFUSIONS_OUT=\"fusions\"\n\nfunction autounzip()\
      \ {\n        if [[ \"$1\" =~ \\.gz$ ]]; then\n                \"$TOOLS_DIR/pigz\"\
      \ -d -c \"$1\"\n        else\n                cat \"$1\"\n        fi\n}\n\n\
      # extract STAR index\nmkdir STAR_index || exit 1\nautounzip \"$STAR_INDEX\"\
      \ | tar -x -C STAR_index -f - || exit 1\n\n# run alignment\n\"$TOOLS_DIR/STAR\"\
      \ \\\n        --runThreadN \"$THREADS\" \\\n        --genomeDir \"$(find -name\
      \ SAindex -printf %h)\" --genomeLoad NoSharedMemory \\\n        --readFilesIn\
      \ \"$READ1\" \"$READ2\" --readFilesCommand zcat \\\n        --outStd BAM_Unsorted\
      \ --outSAMtype BAM Unsorted SortedByCoordinate \\\n        --outFilterMultimapNmax\
      \ 1 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 \\\n        --alignIntronMax\
      \ 500000 --alignMatesGapMax 500000 \\\n        --alignSJstitchMismatchNmax 5\
      \ -1 5 5 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreSeparation 1 --chimScoreJunctionNonGTAG\
      \ 0 --chimJunctionOverhangMin 15 --chimSegmentReadGapMax 3 \\\n        --limitBAMsortRAM\
      \ ${MEMORY}000000 |\n\"$TOOLS_DIR/extract_read-through_fusions\" -g \"$REFERENCE_GTF\"\
      \ -G \"gene_name=gene_name gene_id=gene_id transcript_id=transcript_id gene_status=gene_biotype\
      \ status_KNOWN=protein_coding gene_type=gene_biotype type_protein_coding=protein_coding\
      \ feature_exon=exon feature_UTR=UTR feature_gene=gene\" > read_through.bam ||\
      \ exit 1\n\n# index normal alignments\nmv Aligned.sortedByCoord.out.bam rna.bam\
      \ || exit 1\n\"$TOOLS_DIR/samtools\" index rna.bam & PID_SAMTOOLS_INDEX=$!\n\
      \n# convert chimeric alignments from SAM to BAM\n\"$TOOLS_DIR/samtools\" view\
      \ -Sb Chimeric.out.sam > chimeric.bam & PID_SAMTOOLS_VIEW=$!\n\n# index reference\
      \ genome\nautounzip \"$REFERENCE_GENOME\" > assembly.fa || exit 1\n\"$TOOLS_DIR/samtools\"\
      \ faidx assembly.fa || exit 1\n\nwait $PID_SAMTOOLS_INDEX || exit 1\nwait $PID_SAMTOOLS_VIEW\
      \ || exit 1\n\n# call arriba\n\"$TOOLS_DIR/arriba\" \\\n        -c chimeric.bam\
      \ \\\n        -r read_through.bam \\\n        -x rna.bam \\\n        -o \"$FUSIONS_OUT.tsv\"\
      \ \\\n        -O \"$FUSIONS_OUT.discarded.tsv\" \\\n        -a assembly.fa \\\
      \n        -g \"$REFERENCE_GTF\" \\\n        -G \"gene_name=gene_name gene_id=gene_id\
      \ transcript_id=transcript_id gene_status=gene_biotype status_KNOWN=protein_coding\
      \ gene_type=gene_biotype type_protein_coding=protein_coding feature_exon=exon\
      \ feature_UTR=UTR feature_gene=gene\" \\\n        -b \"$TOOLS_DIR/database/blacklist_hs37d5_gencode19_2017-01-09.tsv.gz\"\
      \ \\\n        -k \"$TOOLS_DIR/database/known_fusions_CancerGeneCensus_gencode19_2017-01-16.tsv.gz\"\
      \ \\\n        -f \"spliced\" || exit 1\n\n# SMC-RNA Challenge-specific filters\n\
      \"$TOOLS_DIR/awk\" '\n        ($1 ~ /^RP(L|S)[0-9]/ || $2 ~ /^RP(L|S)[0-9]/)\
      \ &&\n        $12 > 0 && $13 > 0 && $14 > 0 &&\n        $7 == \"splice-site\"\
      \ && $8 == \"splice-site\" &&\n        $9 != \"deletion/read-through\"\n' \"\
      $FUSIONS_OUT.discarded.tsv\" >> \"$FUSIONS_OUT.tsv\"\n\"$TOOLS_DIR/awk\" '\n\
      \        {\n                homolog1 = $1; homolog2 = $2;\n                sub(/-.*/,\
      \ \"\", homolog1); sub(/-.*/, \"\", homolog2);\n                sub(/([LP]?[0-9]+[A-Z]?(|[0-9]+[A-Z]?))$/,\
      \ \"\", homolog1);\n                sub(/([LP]?[0-9]+[A-Z]?(|[0-9]+[A-Z]?))$/,\
      \ \"\", homolog2);\n        }\n        ($7 == \"splice-site\" && $8 == \"splice-site\"\
      \ || $7 == \"exon\" && $8 == \"exon\" || $12+$13 == 0) &&\n        $3 ~ /\\\
      +\\/\\+|-\\/-|.\\/\\./ && $4 ~ /\\+\\/\\+|-\\/-|.\\/\\./ &&\n        !($9 ~\
      \ /5.-5.|3.-3./) &&\n        !duplicate[$5,$6]++ && !duplicate[$6,$5]++ &&\n\
      \        !duplicate[$1,$2]++ && !duplicate[$2,$1]++ &&\n        homolog1 !=\
      \ homolog2 &&\n        !($1 ~ /^HIST[0-9]/ && $2 ~ /^HIST[0-9]/) &&\n      \
      \  !($1 ~ /^COL[0-9]/ && $2 ~ /^COL[0-9]/) &&\n        !($1 ~ /^NBPF[0-9]/ &&\
      \ $2 ~ /^NBPF[0-9]/) &&\n        !($1 ~ /^AC[0-9]+\\.[0-9]+$/ && $2 ~ /^AC[0-9]+\\\
      .[0-9]+$/) &&\n        $1 != $2\n' \"$FUSIONS_OUT.tsv\" > \"$FUSIONS_OUT.filtered.tsv\"\
      \ || exit 1\n\n# convert to BEDPE\n\"$TOOLS_DIR/awk\" '\nBEGIN{ OFS=\"\\t\"\
      \ }\n{\n        chromosome1 = $5; position1 = $5;\n        sub(/:.*/, \"\",\
      \ chromosome1); sub(/.*:/, \"\", position1);\n        chromosome2 = $6; position2\
      \ = $6;\n        sub(/:.*/, \"\", chromosome2); sub(/.*:/, \"\", position2);\n\
      \        if ($10 == \"downstream\") position1--;\n        if ($11 == \"downstream\"\
      ) position2--;\n        strand1 = $3; if (strand1 ~ /\\.$/) { sub(/\\/.*/, \"\
      \", strand1) } else { sub(/.*\\//, \"\", strand1) }\n        strand2 = $4; if\
      \ (strand2 ~ /\\.$/) { sub(/\\/.*/, \"\", strand2) } else { sub(/.*\\//, \"\"\
      , strand2) }\n        filters = $18;\n        gsub(/a-z_()/, \"\", filters);\n\
      \        split(filters, filtered_reads, \",\");\n        reads = $12+$13+$14;\n\
      \        for (i in filtered_reads) reads += int(filtered_reads[i]);\n      \
      \  print chromosome1, position1-1, position1, chromosome2, position2-1, position2,\
      \ $1\">>\"$2, $15, strand1, strand2, reads;\n}' \"$FUSIONS_OUT.filtered.tsv\"\
      \ > \"$FUSIONS_OUT.bedpe\" || exit 1", filename: run_arriba.sh}
sbg:cmdPreview: bash run_arriba.sh 8 50000 /opt/arriba  /path/to/STAR_index_dir.ext  /path/to/annotation.ext  /path/to/assembly.ext  /path/to/fastq1.ext  /path/to/fastq2.ext
sbg:contributors: [uhrigs]
sbg:createdBy: uhrigs
sbg:createdOn: 1488647203
sbg:id: uhrigs/smc-rna-challenge-3/arriba-0-8/14
sbg:image_url: null
sbg:job:
  allocatedResources: {cpu: 8, mem: 50000}
  inputs:
    REFERENCE_GENOME:
      class: File
      path: /path/to/assembly.ext
      secondaryFiles:
      - {path: .fai}
      size: 0
    REFERENCE_GTF:
      class: File
      path: /path/to/annotation.ext
      secondaryFiles: []
      size: 0
    STAR_INDEX_TAR:
      class: File
      path: /path/to/STAR_index_dir.ext
      secondaryFiles: []
      size: 0
    TUMOR_FASTQ_1:
      class: File
      path: /path/to/fastq1.ext
      secondaryFiles: []
      size: 0
    TUMOR_FASTQ_2:
      class: File
      path: /path/to/fastq2.ext
      secondaryFiles: []
      size: 0
sbg:latestRevision: 14
sbg:license: MIT License
sbg:modifiedBy: uhrigs
sbg:modifiedOn: 1488975714
sbg:project: uhrigs/smc-rna-challenge-3
sbg:projectName: smc-rna-challenge-3
sbg:revision: 14
sbg:revisionsInfo:
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488647203, 'sbg:revision': 0, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488749749, 'sbg:revision': 1, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488750044, 'sbg:revision': 2, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488751914, 'sbg:revision': 3, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488752178, 'sbg:revision': 4, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488754124, 'sbg:revision': 5, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488756175, 'sbg:revision': 6, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488790987, 'sbg:revision': 7, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488793496, 'sbg:revision': 8, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488798798, 'sbg:revision': 9, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488806646, 'sbg:revision': 10, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488811876, 'sbg:revision': 11, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488841163, 'sbg:revision': 12, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488902467, 'sbg:revision': 13, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': uhrigs, 'sbg:modifiedOn': 1488975714, 'sbg:revision': 14, 'sbg:revisionNotes': null}
sbg:sbgMaintained: false
sbg:toolAuthor: Sebastian Uhrig, DKFZ Heidelberg
sbg:toolkitVersion: '0.8'
sbg:validationErrors: []
stdin: ''
stdout: ''
successCodes: []
temporaryFailCodes: []
