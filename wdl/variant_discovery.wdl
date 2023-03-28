## Jorge de la Barrera, 2020
##
## Identify variants for a single sample

workflow variant_discovery{

  String sample_id
  String queue
  String accounting
  Object wf
  Object genome
  String references_dir
  String intervals_fullpath
  String gnomad
  Boolean make_bamout
  Boolean run_ob_filter
  String alignment_in


  call SomaticVC {
    input :
      sample_id = sample_id,
      execution = wf.SomaticVC,
      queue=queue,
      accounting=accounting,
      genome = genome,
      references_dir = references_dir,
      intervals_fullpath = intervals_fullpath,
      gnomad = gnomad,
      make_bamout = make_bamout,
      run_ob_filter = run_ob_filter,
      alignment_in = alignment_in
  }

  call create_gVCF {
    input :
      sample_id = sample_id,
      execution = wf.create_gVCF,
      queue = queue,
      accounting = accounting,
      genome = genome,
      references_dir = references_dir,
      intervals_fullpath = intervals_fullpath,
      alignment_in = alignment_in
  }

  call Bam2Cram {
    input :
      sname = sample_id,
      execution = wf.create_gVCF,
      queue = queue,
      accounting = accounting,
      genome = genome,
      references_dir = references_dir,
      alignment_in = alignment_in
  }

  output {
    File vcf_out = SomaticVC.vcf_out
    File vcf_out_idx = SomaticVC.vcf_out_idx
    File vcf_out_stats = SomaticVC.vcf_out_stats
    #File bamout = SomaticVC.bamout
    #File bamout_idx = SomaticVC.bamout_idx
    #File f1r2 = SomaticVC.f1r2
    #File vcf_germinal_out = GerminalVC.vcf_out
    #File vcf_germinal_out_idx = GerminalVC.vcf_out_idx
    #File vcf_germinal_hq_out = GerminalVC.vcf_hq_out
    File gvcf_out = create_gVCF.gvcf_out
    File cram_file = Bam2Cram.cram_file
  }
}
  
task SomaticVC {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String references_dir
  String intervals_fullpath
  String gnomad
  Boolean make_bamout
  Boolean run_ob_filter
  String alignment_in
  
  String genome_fasta = references_dir + "/" + genome.path + "/" + genome.ref_fasta
  String gnomad_fullpath = references_dir + "/" + gnomad
  
  command {
      # Load enviroment
      source /programs/GATK/env_vep.sh

      set -e

      gatk --java-options "${execution.java_args_SomaticVC}" Mutect2 \
          ${"-L " + intervals_fullpath} \
          ${"-R " + genome_fasta} \
          ${"-I " + alignment_in} \
          ${"--germline-resource " + gnomad_fullpath} \
          -O "${sample_id}.tmp.vcf.gz" \
          ${true='--bam-output bamout.bam' false='' make_bamout} \
          ${true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter}

      mv ${sample_id}.tmp.vcf.gz.stats ${sample_id}.vcf.gz.stats

      gatk RenameSampleInVcf \
      -I ${sample_id}.tmp.vcf.gz \
      -O ${sample_id}.vcf.gz \
      --NEW_SAMPLE_NAME ${sample_id}

    bcftools index -t ${sample_id}.vcf.gz

  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File vcf_out = "${sample_id}.vcf.gz"
    File vcf_out_idx = "${sample_id}.vcf.gz.tbi"
    File vcf_out_stats = "${sample_id}.vcf.gz.stats"
    #File bamout = "bamout.bam"
    #File bamout_idx = "bamout.bai"
    #File f1r2 = "f1r2.tar.gz"
  }
}

task GerminalVC {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String references_dir
  String intervals_fullpath
  String alignment_in
  
  String genome_fasta = references_dir + "/" + genome.path + "/" + genome.ref_fasta
  String dbsnp = references_dir + "/" + genome.path + "/" + genome.dbSNP_vcf
  
  command <<<
      # Load enviroment
      source /programs/GATK/env_vep.sh

      set -e

      gatk --java-options "${execution.java_args_HaplotypeCaller}" HaplotypeCaller \
        ${"-I " + alignment_in} \
        ${"-R " + genome_fasta} \
        -O "${sample_id}.tmp.vcf.gz" \
        ${"-L " + intervals_fullpath} \
        --max-reads-per-alignment-start 0 \
        ${"--dbsnp " + dbsnp}

      gatk RenameSampleInVcf \
        -I ${sample_id}.tmp.vcf.gz \
        -O ${sample_id}.germline.vcf.gz \
        --NEW_SAMPLE_NAME ${sample_id}

      bcftools index -t ${sample_id}.germline.vcf.gz
      bcftools isec ${sample_id}.germline.vcf.gz /data3/genomes/homo_sapiens/GATK_bundle/BROAD_hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -n =2 -w 1 > ${sample_id}.germline.HQ.vcf.gz

    >>>
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File vcf_out = "${sample_id}.germline.vcf.gz"
    File vcf_out_idx = "${sample_id}.germline.vcf.gz.tbi"
    File vcf_hq_out = "${sample_id}.germline.HQ.vcf.gz"
  }
}

task create_gVCF {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String references_dir
  String intervals_fullpath
  String alignment_in
  
  String genome_fasta = references_dir + "/" + genome.path + "/" + genome.ref_fasta
  String dbsnp = references_dir + "/" + genome.path + "/" + genome.dbSNP_vcf
  
  command {
      # Load enviroment
      source /programs/GATK/env.sh

      set -e

      gatk --java-options "${execution.java_args_HaplotypeCaller}" HaplotypeCaller \
        ${"-I " + alignment_in} \
        ${"-R " + genome_fasta} \
        -O "${sample_id}.tmp.g.vcf.gz" \
        ${"-L " + intervals_fullpath} \
        --max-reads-per-alignment-start 0 \
        -ERC GVCF \
        --pair-hmm-implementation ${execution.gatk_gkl_pairhmm_implementation} \
        --native-pair-hmm-threads ${execution.gatk_gkl_pairhmm_threads} \
        --smith-waterman ${execution.smith_waterman_implementation} && \
        gatk RenameSampleInVcf \
          -I ${sample_id}.tmp.g.vcf.gz \
          -O ${sample_id}.g.vcf.gz \
          --NEW_SAMPLE_NAME ${sample_id} && \
        rm ${sample_id}.tmp.g.vcf.gz       

  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File gvcf_out = "${sample_id}.g.vcf.gz"
  }
}

task Bam2Cram {

    Object execution
    String queue
    String accounting
    Object genome
    String sname
    String references_dir
    File alignment_in
    #File alignment_index_in

    File genome_fasta = references_dir + "/" + genome.path + "/" + genome.ref_fasta
    command {

    source /programs/GATK/env.sh

    #samtools
    samtools view -C -T ${genome_fasta} -o ${sname}.cram ${alignment_in}
    samtools index ${sname}.cram

    exit $rc

    }
    runtime {
        backend : 'SGE_nope'
        cpu : execution.cpu
        memory : execution.memory
        accounting : accounting
        sge_project : queue
    }
    output {
    File cram_file = "${sname}.cram"
    File cram_index = "${sname}.cram.crai"
    }
}