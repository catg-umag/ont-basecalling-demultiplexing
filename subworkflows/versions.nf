workflow CollectVersions {
  take:
    basecalled_ubam

  main:
    (dorado & samtools & fastQC & nanoq & nanoPlot & pycoQC & toulligQC)
      | mix
      | set { software_versions }

    software_versions
      | collectFile( name: 'software_versions.yaml', newLine: true, sort: true) {
          "${it[0]}: \"${it[1]}\""
        }
      | set { software_versions_combined }

    doradoModel(basecalled_ubam)
  
  emit:
    software_versions = software_versions_combined
    model_versions = doradoModel.out
}


process dorado {
  label 'dorado'

  output:
  tuple val('Dorado'), stdout

  script:
  """
  dorado --version 2>&1 | tr -d '\n'
  """
}


process samtools {
  label 'samtools'

  output:
  tuple val('Samtools'), stdout

  script:
  """
  samtools --version | head -n 1 | grep -Eo '[0-9.]+' | tr -d '\n'
  """
}


process fastQC {
  label 'fastqc'

  when:
  'fastqc' in params.qc_tools

  output:
  tuple val('FastQC'), stdout

  script:
  """
  fastqc --version | grep -Eo '[0-9.]+' | tr -d '\n'
  """
}


process nanoPlot {
  label 'nanoplot'

  when:
  'nanoplot' in params.qc_tools

  output:
  tuple val('NanoPlot'), stdout

  script:
  """
  NanoPlot --version | grep -Eo '[0-9.]+' | tr -d '\n'
  """
}


process nanoq {
  label 'nanoq'

  when:
  'nanoq' in params.qc_tools

  output:
  tuple val('nanoq'), stdout

  script:
  """
  nanoq --version | grep -Eo '[0-9.]+' | tr -d '\n'
  """
}


process pycoQC {
  label 'pycoqc'

  when:
  'pycoqc' in params.qc_tools

  output:
  tuple val('pycoQC'), stdout

  script:
  """
  pycoQC --version | grep -Eo '[0-9.]+' | tr -d '\n'
  """
}


process toulligQC {
  label 'toulligqc'

  when:
  'toulligqc' in params.qc_tools
  
  output:
  tuple val('toulligQC'), stdout
  
  script:
  """
  toulligqc --version | tr -d '\n'
  """
}


process doradoModel {
  label 'samtools'

  input:
  path(bam)
  
  output:
  path('dorado_model.tsv')
  
  script:
  """
  model_version=\$(samtools view -H ${bam} | grep -Po '(?<=basecall_model=)([^ ]+)' | uniq)
  echo "Software\tModel\tVersion" > dorado_model.tsv
  echo "Dorado\tBasecalling\t\${model_version}" >> dorado_model.tsv 
  """
}
