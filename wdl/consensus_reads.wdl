## Jorge de la Barrera, 2019

workflow consensus_reads{

  String sample_id
  String rootdir
  String queue
  String accounting
  Object wf
  Object genome
  String raw_dir
  String references_dir
  String bin_dir
  String enrichment_intervals
  Array[String] reads_in


  scatter (sample_reads_in in reads_in) { 
    call uBAM_2_mappedBAM {
      input :
        rootdir = rootdir,
        execution = wf.uBAM_2_mappedBAM,
        queue=queue,
        accounting=accounting,
        genome = genome,
        raw_dir = raw_dir,
        references_dir = references_dir,
        reads_in = sample_reads_in
    }
  }

  call MarkDuplicates {
    input :
      sample_id = sample_id,
      rootdir = rootdir,
      execution = wf.MarkDuplicates,
      queue=queue,
      accounting=accounting,
      alignment_in =  uBAM_2_mappedBAM.alignment,
      alignment_index_in =  uBAM_2_mappedBAM.alignment_index
  }

  call GroupReadsByUmi {
    input :
      rootdir = rootdir,
      execution = wf.GroupReadsByUmi,
      queue=queue,
      accounting=accounting,
      alignment_in =  MarkDuplicates.alignment,
      bin_dir = bin_dir
  }

  call TemplateSortBam {
    input :
      rootdir = rootdir,
      execution = wf.TemplateSortBam,
      queue=queue,
      accounting=accounting,
      reads_in=GroupReadsByUmi.alignment,
      bin_dir = bin_dir
  }

  call CallMolecularConsensusReads {
    input :
      rootdir = rootdir,
      execution = wf.CallMolecularConsensusReads,
      queue=queue,
      accounting=accounting,
      reads_in=TemplateSortBam.templated_sorted,
      bin_dir = bin_dir
  }

  call mapping {
    input :
      rootdir = rootdir,
      execution = wf.mapping,
      queue=queue,
      accounting=accounting,
      genome = genome,
      reads_in = CallMolecularConsensusReads.consensus_reads,
      references_dir = references_dir
  }

  call FilterConsensusReads {
    input :
      rootdir = rootdir,
      execution = wf.FilterConsensusReads,
      queue=queue,
      accounting=accounting,
      genome = genome,
      alignment_in = mapping.alignment,
      alignment_index_in = mapping.alignment_index,
      bin_dir = bin_dir,
      references_dir = references_dir
  }

  call ClipBam {
    input :
      rootdir = rootdir,
      execution = wf.ClipBam,
      queue=queue,
      accounting=accounting,
      genome = genome,
      sample_id = sample_id,
      alignment_in = FilterConsensusReads.filtered_alignment,
      bin_dir = bin_dir,
      references_dir = references_dir
  }

  # Stats

  # uBAM mapped
  call getBarcodes {
    input :
      rootdir = rootdir,
      execution = wf.getBarcodes,
      queue=queue,
      accounting=accounting,
      alignment_in=GroupReadsByUmi.alignment
  }

  call QualityScoreDistribution {
    input :
      rootdir = rootdir,
      execution = wf.QualityScoreDistribution,
      queue=queue,
      accounting=accounting,
      alignment_in=GroupReadsByUmi.alignment
  }

  call CollectHsMetrics as raw_HsMetrics {
    input :
      rootdir = rootdir,
      execution = wf.CollectHsMetrics,
      queue=queue,
      accounting=accounting,
      genome=genome,
      alignment_in=GroupReadsByUmi.alignment,
      min_base_qual = 0,
      min_mapping_qual = 0,
      references_dir = references_dir,
      enrichment_intervals = enrichment_intervals
  }

    call CollectHsMetrics as consensus_HsMetrics {
      input :
        rootdir = rootdir,
        execution = wf.CollectHsMetrics,
        queue=queue,
        accounting=accounting,
        genome=genome,
        alignment_in=out_alignment,
        min_base_qual = 0,
        min_mapping_qual = 0,
        references_dir = references_dir,
        enrichment_intervals = enrichment_intervals
  }

  call BAM_flagstat {
    input:
      rootdir = rootdir,
      execution =wf.BAM_flagstat,
      queue = queue,
      accounting = accounting,
      bamfile = out_alignment,
      bamfile_index = out_alignment_index
  }
  output {
    #File alignment = uBAM_2_mappedBAM.alignment
    #File alignment_index = uBAM_2_mappedBAM.alignment_index
    Array[File] raw_aligment_log = uBAM_2_mappedBAM.aligment_log
    File markduplicates_metrics = MarkDuplicates.markduplicates_metrics 
    File umi_metrics = MarkDuplicates.umi_metrics
    File hist = GroupReadsByUmi.hist
    File consensus_alignment_log = mapping.alignment_log
    File flagstat = BAM_flagstat.flagstat
    File RX_barcode = getBarcodes.RX_barcode
    File qual_score_dist = QualityScoreDistribution.qual_score_dist
    File raw_hs_metrics = raw_HsMetrics.hs_metrics
    File raw_target_coverage = raw_HsMetrics.target_coverage
    File raw_base_coverage = raw_HsMetrics.base_coverage
    File consensus_hs_metrics = consensus_HsMetrics.hs_metrics
    File consensus_target_coverage = consensus_HsMetrics.target_coverage
    File consensus_base_coverage = consensus_HsMetrics.base_coverage

    File out_alignment = ClipBam.clipped_alignment
    File out_alignment_index = ClipBam.clipped_alignment_index
  }
}

task uBAM_2_mappedBAM {

  String rootdir
  Object execution
  String queue
  String accounting
  Object genome
  String raw_dir
  String references_dir
  String reads_in

  String file_prefix  = basename (reads_in, ".bam") + ".mapped"

  String ref_fasta = "${references_dir}/${genome.path}/${genome.ref_fasta}"

  String bwa_commandline="-K 100000000 -p -Y -v 3 -t ${execution.cpu} ${ref_fasta}"
  
  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    # job
    picard ${execution.java_args_SamToFastq} SamToFastq \
      TMP_DIR=$TMPDIR \
      I=${reads_in} \
      F=/dev/stdout \
      INTERLEAVE=true \
      INCLUDE_NON_PF_READS=false \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=X| \
    bwa mem ${bwa_commandline} /dev/stdin - 2>${file_prefix}.bwa.stderr.log | \
    picard ${execution.java_args_MergeBamAlignment} MergeBamAlignment \
        TMP_DIR=$TMPDIR \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${reads_in} \
        R=${ref_fasta} \
        O=${file_prefix}.bam \
        EXPECTED_ORIENTATIONS=FR \
        PAIRED_RUN=true \
        SORT_ORDER="coordinate" \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        ALIGNER_PROPER_PAIR_FLAGS=false \
        CREATE_INDEX=true && \
    samtools quickcheck -q ${file_prefix}.bam && \
    grep -m1 "read .* ALT contigs" ${file_prefix}.bwa.stderr.log | grep -v "read 0 ALT contigs"

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${file_prefix}.bam"
    File alignment_index = "${file_prefix}.bai"
    File aligment_log = "${file_prefix}.bwa.stderr.log"
  }
}

task MarkDuplicates {

  String sample_id
  String rootdir
  Object execution
  String queue
  String accounting
  Array[File] alignment_in
  Array[File] alignment_index_in

  String file_prefix = "${sample_id}.markdup"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    picard ${execution.java_args_MarkDuplicates} UmiAwareMarkDuplicatesWithMateCigar \
      TMP_DIR=$TMPDIR \
      I=${sep=' I=' alignment_in} \
      O=${file_prefix}.bam \
      M=${file_prefix}.markduplicates_metrics \
      UMI_METRICS=${file_prefix}.umi_metrics \
      MAX_EDIT_DISTANCE_TO_JOIN=1 \
      TAGGING_POLICY=All \
      TAG_DUPLICATE_SET_MEMBERS=true \
      REMOVE_SEQUENCING_DUPLICATES=true \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ADD_PG_TAG_TO_READS=false && \
    samtools quickcheck -q ${file_prefix}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File markduplicates_metrics = "${file_prefix}.markduplicates_metrics"
    File umi_metrics = "${file_prefix}.umi_metrics"
    File alignment = "${file_prefix}.bam"
  }
}

task GroupReadsByUmi {
  String rootdir
  Object execution
  String queue
  String accounting
  File alignment_in
  String bin_dir

  String file_prefix  = basename (alignment_in, ".bam") + ".grouped"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    java ${execution.java_args_GroupReadsByUmi} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    GroupReadsByUmi \
      -i ${alignment_in} \
      -f ${file_prefix}.hist \
      -o ${file_prefix}.bam \
      --include-non-pf-reads=false \
      --strategy=${execution.grouping_strategy} \
      --edits=${execution.edit_distance} && \
    samtools quickcheck -q ${file_prefix}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"
      
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File hist = "${file_prefix}.hist"
    File alignment = "${file_prefix}.bam"
  }
}

task TemplateSortBam {
  String rootdir
  Object execution
  String queue
  String accounting
  File reads_in
  String bin_dir

  String file_prefix  = basename (reads_in, ".bam") + ".sorted"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    java ${execution.java_args_SortBam} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    SortBam \
    --input=${reads_in} \
    --sort-order=TemplateCoordinate \
    --output ${file_prefix}.bam
      
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File templated_sorted = "${file_prefix}.bam"
  }
}

task CallMolecularConsensusReads {
  String rootdir
  Object execution
  String queue
  String accounting
  File reads_in
  String bin_dir

  String file_prefix  = basename (reads_in, ".bam") + ".consensus"

  command {
    # Load enviroment
    source /programs/GATK/env.sh


    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    java ${execution.java_args_CallMolecularConsensusReads} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    CallMolecularConsensusReads \
      -i ${reads_in} \
      -o ${file_prefix}.bam \
      --min-reads=1
      
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File consensus_reads = "${file_prefix}.bam"
  }
}

task mapping {
  String rootdir
  Object execution
  String queue
  String accounting
  Object genome
  File reads_in
  String references_dir

  String ref_fasta = "${references_dir}/${genome.path}/${genome.ref_fasta}"

  String bwa_commandline="-K 100000000 -p -Y -v 3 -t ${execution.cpu} ${ref_fasta}"

  String file_prefix  = basename (reads_in, ".bam") + ".mapped"


  command {

    # Load enviroment
    source /programs/GATK/env.sh


    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    # job
    picard ${execution.java_args_SamToFastq} SamToFastq \
      TMP_DIR=$TMPDIR \
      I=${reads_in} \
      F=/dev/stdout \
      INTERLEAVE=true \
      INCLUDE_NON_PF_READS=false| \
    bwa mem ${bwa_commandline} /dev/stdin - 2>${file_prefix}.bwa.stderr.log | \
    picard ${execution.java_args_MergeBamAlignment} MergeBamAlignment \
        TMP_DIR=$TMPDIR \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${reads_in} \
        R=${ref_fasta} \
        O=${file_prefix}.bam \
        EXPECTED_ORIENTATIONS=FR \
        PAIRED_RUN=true \
        SORT_ORDER="coordinate" \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        ALIGNER_PROPER_PAIR_FLAGS=false \
        CREATE_INDEX=true && \
    samtools quickcheck -q ${file_prefix}.bam && \
    grep -m1 "read .* ALT contigs" ${file_prefix}.bwa.stderr.log | grep -v "read 0 ALT contigs"

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${file_prefix}.bam"
    File alignment_index = "${file_prefix}.bai"
    File alignment_log = "${file_prefix}.bwa.stderr.log"
  }
}

task FilterConsensusReads {
  String rootdir
  Object execution
  String queue
  String accounting
  Object genome
  File alignment_in
  File alignment_index_in
  String bin_dir
  String references_dir

  String ref_fasta = "${references_dir}/${genome.path}/${genome.ref_fasta}"

  String file_prefix  = basename (alignment_in, ".bam") + ".filtered"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    java ${execution.java_args_FilterConsensusReads} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    FilterConsensusReads \
      -i ${alignment_in} \
      -o ${file_prefix}.bam \
      --ref=${ref_fasta} \
      --min-reads=${execution.supporting_min_reads} \
      --min-base-quality=${execution.supporting_min_base_quality} \
      --reverse-per-base-tags=${execution.reverse_per_base_tags} \
      --require-single-strand-agreement=${execution.require_single_strand_agreement} && \
    samtools quickcheck -q ${file_prefix}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File filtered_alignment = "${file_prefix}.bam"
  }
}

task ClipBam {
  String rootdir
  Object execution
  String queue
  String accounting
  Object genome
  String sample_id
  File alignment_in
  String bin_dir
  String references_dir

  String ref_fasta = "${references_dir}/${genome.path}/${genome.ref_fasta}"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    java ${execution.java_args_ClipBam} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    ClipBam \
      -i ${alignment_in} \
      -o ${sample_id}.bam \
      --ref=${ref_fasta} \
      --clip-overlapping-reads=${execution.clip_overlapping_reads} \
      --clipping-mode=${execution.clipping_mode} && \
    samtools quickcheck -q ${sample_id}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output { 
    File clipped_alignment = "${sample_id}.bam"
    File clipped_alignment_index = "${sample_id}.bai"
  }
}

task BAM_flagstat {

  String rootdir
  Object execution
  String queue
  String accounting
  File bamfile
  File bamfile_index

  String file_prefix = basename( bamfile , ".bam")

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    samtools flagstat ${bamfile} > ${file_prefix}.flagstat

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
  File flagstat = "${file_prefix}.flagstat"
  }
}

task getBarcodes {
  String rootdir
  Object execution
  String queue
  String accounting
  File alignment_in

  String file_prefix = basename( alignment_in , ".bam")

  command <<<

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh


    # 0x0040: first in pair
    # 0x0100: not primary alignment 


    samtools view -f 0x0040 -F 0x0100 ${alignment_in}|awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^RX:Z:"){ s = substr($i,6,length($i)-5); }; print s }}'|sort|uniq > ${file_prefix}.RX.barcode
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  >>>
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File RX_barcode = "${file_prefix}.RX.barcode"
  }
}

task QualityScoreDistribution {
  String rootdir
  Object execution
  String queue
  String accounting
  File alignment_in

  String file_prefix = basename( alignment_in , ".bam")

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # Debug
    #sh ${rootdir}/0000_pipeline/_debug.sh

    picard ${execution.java_args_QualityScoreDistribution} QualityScoreDistribution \
      I=${alignment_in} \
      O=${file_prefix}.qual_score_dist \
      CHART=${file_prefix}.qual_score_dist.pdf \
      PF_READS_ONLY=true

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File qual_score_dist = "${file_prefix}.qual_score_dist"
  }
}

task CollectHsMetrics {
  String rootdir
  Object execution
  String queue
  String accounting
  Object genome
  File alignment_in
  Int min_base_qual
  Int? min_mapping_qual=20
  String references_dir
  String enrichment_intervals

  String file_prefix = basename(alignment_in, ".bam")

  String ref_fasta = "${references_dir}/${genome.path}/${genome.ref_fasta}"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    picard ${execution.java_args_CollectHsMetrics} CollectHsMetrics \
      I=${alignment_in} \
      O=${file_prefix}.hs_metrics \
      R=${ref_fasta} \
      BAIT_INTERVALS=${enrichment_intervals} \
      TARGET_INTERVALS=${enrichment_intervals} \
      PER_TARGET_COVERAGE=${file_prefix}.target_coverage \
      PER_BASE_COVERAGE=${file_prefix}.base_coverage  \
      COVERAGE_CAP=${execution.covarage_cap} \
      MINIMUM_BASE_QUALITY=${min_base_qual} \
      MINIMUM_MAPPING_QUALITY=${min_mapping_qual}  

    # Saca el RC
    echo "ExitCode:$?"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File hs_metrics ="${file_prefix}.hs_metrics"
    File target_coverage = "${file_prefix}.target_coverage"
    File base_coverage = "${file_prefix}.base_coverage"
  }
}



