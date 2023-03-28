## Jorge de la Barrera, 2021
##
## Identify germline variants for a single sample

workflow germline_calling{

  String sample_id
  String rootdir
  String queue
  String accounting
  Object wf
  Object genome
  String references_dir
  String intervals_fullpath
  String alignment_in


  call GerminalVC {
    input :
      sample_id = sample_id,
      rootdir = rootdir,
      execution = wf.GerminalVC,
      queue=queue,
      accounting=accounting,
      genome = genome,
      references_dir = references_dir,
      intervals_fullpath = intervals_fullpath,
      alignment_in = alignment_in
  }
  output {
    File vcf_out = GerminalVC.vcf_out
    File vcf_out_idx = GerminalVC.vcf_out_idx
    File vcf_hq_out = GerminalVC.vcf_hq_out 
  }
}
  
task GerminalVC {
  String sample_id
  String rootdir
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

    gatk --java-options "${execution.java_args_GerminalVC}" HaplotypeCaller \
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
