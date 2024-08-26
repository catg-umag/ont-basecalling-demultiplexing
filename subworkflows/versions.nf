workflow CollectVersions {
  take:
    basecalled_ubam

  main:
    (dorado & samtools & fastQC & nanoq & nanoPlot & pycoQC & toulligQC)
      | mix
      | set { software_versions }

    software_versions
      | collectFile( name: 'software_versions.yaml', newLine: false, sort: true) {
          "${it[0]}: ${it[1]}"
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
  dorado --version 2>&1
  """
}


process samtools {
  label 'samtools'

  output:
  tuple val('Samtools'), stdout

  script:
  """
  samtools --version | head -n 1 | grep -Eo '[0-9.]+'
  """
}


process fastQC {
  label 'fastqc'

  output:
  tuple val('FastQC'), stdout

  script:
  """
  fastqc --version | grep -Eo '[0-9.]+'
  """
}


process nanoPlot {
  label 'nanoplot'

  output:
  tuple val('NanoPlot'), stdout

  script:
  """
  NanoPlot --version | grep -Eo '[0-9.]+'
  """
}


process nanoq {
  label 'nanoq'

  output:
  tuple val('nanoq'), stdout

  script:
  """
  nanoq --version | grep -Eo '[0-9.]+'
  """
}


process pycoQC {
  label 'pycoqc'

  output:
  tuple val('PycoQC'), stdout

  script:
  """
  pycoQC --version | grep -Eo '[0-9.]+'
  """
}


process toulligQC {
  label 'toulligqc'
  
  output:
  tuple val('ToulligQC'), stdout
  
  script:
  """
  toulligqc --version
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
  model_version=\$(samtools view -H ${bam} | grep -Po '(?<=basecall_model=)([^ ]+)')
  echo "Software\tModel\tVersion" > dorado_model.tsv
  echo "Dorado\tBasecalling\t\${model_version}" >> dorado_model.tsv 
  """
}
