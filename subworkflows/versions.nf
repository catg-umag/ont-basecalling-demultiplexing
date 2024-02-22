workflow CollectVersions {
  main:
      dorado()
        | mix(samtools())
        | mix(fastQC())
        | mix(nanoPlot())
        | mix(pycoQC())
        | set { software_versions }

      software_versions
        | map { it.text }
        | collectFile(
            name: 'software_versions.yaml',
            newLine: false,
            sort: true
          )

      doradoModel()
  
  emit:
    software_versions = software_versions.collect()
    model_versions = doradoModel.out
}


process dorado {
  label 'dorado'

  output:
  path('dorado_version.yaml')

  script:
  """
  version=\$(dorado --version 2>&1)
  echo "Dorado: '\$version'" > dorado_version.yaml
  """
}


process samtools {
  label 'samtools'

  output:
  path('samtools_version.yaml')

  script:
  """
  version=\$(samtools --version | head -n 1 | grep -Eo '[0-9.]+')
  echo "Samtools: '\$version'" > samtools_version.yaml
  """
}


process fastQC {
  label 'fastqc'

  output:
  path('fastqc_version.yaml')

  script:
  """
  version=\$(fastqc --version | grep -Eo '[0-9.]+')
  echo "FastQC: '\$version'" > fastqc_version.yaml
  """
}


process nanoPlot {
  label 'nanoplot'

  output:
  path('nanoplot_version.yaml')

  script:
  """
  version=\$(NanoPlot --version | grep -Eo '[0-9.]+') 
  echo "NanoPlot: '\$version'" > nanoplot_version.yaml
  """
}


process pycoQC {
  label 'pycoqc'

  output:
  path('pycoqc_version.yaml')

  script:
  """
  version=\$(pycoQC --version | grep -Eo '[0-9.]+') 
  echo "PycoQC: '\$version'" > pycoqc_version.yaml
  """
}


process doradoModel {
  output:
  path('dorado_model.tsv')
  
  script:
  model_version = params.dorado_basecalling_model
  """
  echo "Software\tModel\tVersion" > dorado_model.tsv
  echo "Dorado\tBasecalling\t${model_version}" >> dorado_model.tsv 
  """
}