#configfile: "snakemake_config.yaml"

if not os.path.exists(config["output_dir"]):
    os.makedirs(config["output_dir"])

result_dir = os.path.join(config["output_dir"], "Results/")
filtered_vcf = config["basename"] + "_filtered_" + config["variant_caller"] + ".vcf"
raw_vcf_input = config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz"

if "vcf_filter" in config:
  vcf_filter_params = ' '.join(['%s %s' % (key, value) for (key, value) in config["vcf_filter"].items()])
else:
  vcf_filter_params = ""

if "flanking_sequences" in config:
  flanking_sequences = ' '.join(['%s %s' % (key, value) for (key, value) in config["flanking_sequences"].items()])
else:
  flanking_sequences = ""

if "pair_file" in config:
  pair_file = "-p " + config["pair_file"]
else:
  pair_file = ""

if config["method"] == "frequency":
  filter_script = "method_frequency/vcf_freq_filter.py"
  select = "method_frequency/select_frequency_markers.py"
  configure = "method_frequency/configure_dist_set.R"

if config["method"] == "snps":
  filter_script = "method_snps/vcf_filter.py"
  select = "method_snps/new_select_markers.py"
  configure = "method_snps/configure_snp_set.R"

rule all:
  input:
    result_dir + "snp_set_" + config["basename"] + ".vcf"

rule RG_tag_to_header:
  input:
    config["input_dir"]
  output:
    config["output_dir"] + config["basename"] + ".fasta.gz",
    config["output_dir"] + "RG_tag_list_" + config["basename"] + ".txt"
  shell:
    "python3 normal_RG_tag_to_header.py -d {input} -o {output[0]}"

rule Bwa_mem:
  input:
    config["reference_genome"],
    config["output_dir"] + config["basename"] + ".fasta.gz"
  threads:
    8
  output:
    config["output_dir"] + config["basename"] + ".sam"
  shell:
    "bwa mem -t {threads} -C {input[0]} {input[1]} > {output}"

rule RG_to_sam_file_header:
  input:
    config["output_dir"] + config["basename"] + ".sam",
    config["output_dir"] + "RG_tag_list_" + config["basename"] + ".txt"
  output:
    config["output_dir"] + config["basename"] + ".bam",
    config["output_dir"] + "sorted_" + config["basename"] + ".bam"
  shell:
    "python3 RG_to_sam.py -s {input[0]} -r {input[1]} -o {output[0]}"

if config["variant_caller"] == "samtools":

  rule samtools:
    input:
      config["reference_genome"],
      config["output_dir"] + "sorted_" + config["basename"] + ".bam"
    threads:
      8
    output:
      config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz"
    shell:
        "samtools mpileup -uf {input[0]} {input[1]} -d 1000000 --output-tags DP,AD,ADF,ADR,SP | bcftools call -mv --threads {threads} | bgzip -c > {output}"

elif config["variant_caller"] == "freebayes":

  rule freebayes:
    input:
      config["reference_genome"][:-3],
      config["output_dir"] + "sorted_" + config["basename"] + ".bam"
    threads:
      8
    output:
      config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz"
    shell:
        "freebayes -f {input[0]} {input[1]} | bgzip -c > {output}"

rule filter_vcf:
  input:
    config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz",
  params:
    filter_script,
    vcf_filter_params,
    pair_file
  output:
      config["output_dir"] + filtered_vcf
  shell:
    "python3 {params[0]} -v {input[0]} -o {output} {params[1]} {params[2]}"

rule bgzip:
  input:
    config["output_dir"] + filtered_vcf
  output:
    config["output_dir"] + filtered_vcf + ".gz"
  shell:
    "bgzip -c {input} > {output}"


rule tabix:
  input:
    config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz"
  output:
    config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz" + ".tbi"
  shell:
    "tabix -p vcf {input}"

rule get_markers:
  input:
    config["output_dir"] + filtered_vcf + ".gz",
  params:
    select,
    pair_file
  output:
    config["output_dir"] + filtered_vcf[:-4] + ".csv"
  shell:
    "python3 {params[0]} -v {input[0]} -o {output} {params[1]}"

rule make_snp_set:
  input:
    config["output_dir"] + filtered_vcf[:-4] + ".csv",
  params:
    configure,
    result_dir,
    config["configure_snp_set"]["genetic_distance"],
  output:
    result_dir + "snps_" + str(config["configure_snp_set"]["genetic_distance"]) + "_" + filtered_vcf[:-4] + ".csv"
  shell:
    "Rscript {params[0]} {input[0]} {params[1]} {params[2]}"


rule get_flanking_sequences:
  input:
    config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz",
    config["reference_genome"][:-3],
    result_dir + "snps_" + str(config["configure_snp_set"]["genetic_distance"]) + "_" + filtered_vcf[:-4] + ".csv",
    config["output_dir"] + config["basename"] + "_" + config["variant_caller"] + ".vcf.gz" + ".tbi"
  params:
    flanking_sequences
  output:
    result_dir + "flanking_sequences_" + config["basename"] + ".tsv",
    result_dir + "snp_set_" + config["basename"] + ".vcf"
  shell:
    "python3 flanking_sequences.py -v {input[0]} -r {input[1]} -s {input[2]} -f {output[0]} -c {output[1]} {params}"

