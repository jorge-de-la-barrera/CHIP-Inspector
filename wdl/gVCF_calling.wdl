## Jorge de la Barrera, 2021
##
## Identify germline variants for a single sample

workflow gVCF_calling{

  String sample_id
  String rootdir
  String queue
  String accounting
  Object wf
  Object genome
  String references_dir
  String intervals_fullpath
  String alignment_in


  call HaplotypeCalling {
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
    File gvcf_out = HaplotypeCalling.gvcf_out
    File gvcf_index = HaplotypeCalling.gvcf_index
  }
}
  
task HaplotypeCalling {
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
  
command <<<
    # Load enviroment
    source /programs/GATK/env_vep.sh

    set -e

    gatk --java-options "${execution.java_args_HaplotypeCalling}" HaplotypeCaller \
      ${"-I " + alignment_in} \
      ${"-R " + genome_fasta} \
      -O "${sample_id}.tmp.vcf.gz" \
      ${"-L " + intervals_fullpath} \
      --max-reads-per-alignment-start 0 \
      -ERC GVCF \
      --enable-all-annotations true 
  >>>
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File gvcf_out = "${sample_id}.germline.g.vcf.gz"
    File gvcf_index = "${sample_id}.germline.g.vcf.gz.tbi"
  }
}
