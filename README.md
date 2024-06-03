# ONT Basecalling / Demux Pipeline

Small pipeline to perform basecalling and demultiplexing (optional) of ONT data, collect QC metrics and generate a MultiQC report.
It uses Dorado for basecalling and demultiplexing.

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.04)
- [Apptainer](https://apptainer.org/) / Singularity
- Dorado (0.5.3 tested). It can be used via container, or installed locally from https://github.com/nanoporetech/dorado.

## Usage

- Clone this repository
- **If you want to demultiplex:** create a `samples.csv` file with at least the `barcode` and `sample` columns. The `barcode` column should contain the barcode used for demultiplexing (with the leading zero, e.g. `barcode01`), and the `sample` column should contain the sample name (this name with be used on the report and as name for FASTQ file).
- Copy `params.example.yml` (for example to `./my_params.yml`) and modify it according to your needs. Remember to point `sample_data` parameter to the file created at the previous step.
- Run the pipeline passing your params file to `-params-file` option:

```
nextflow run ont-basecalling-demultiplexing/ -profile apptainer -params-file my_params.yml
```

## Parameters

| Parameter                                | Required | Default                              | Description                                                                                                             |
| ---------------------------------------- | -------- | ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------- |
| `experiment_name`                        | False    | -                                    | Name of the experiment, used for final reports (title and filename).                                                    |
| `data_dir`                               | True     | -                                    | Path to the folder containing the POD5 files.                                                                           |
| `sample_data`                            | True     | `input/samples.csv`                  | Path to the CSV file containing the sample data (required if demultiplexing).                                           |
| `output_dir`                             | False    | `demultiplex_results`                | Path to the folder where the results will be saved.                                                                     |
| `fastq_output`                           | False    | `true`                               | If `true`, the pipeline will generate FASTQ files (if not, it would be UBAM files).                                     |
| `qscore_filter`                          | False    | `10`                                 | Minimum QScore for the "pass" data, used for demultiplexing.                                                            |
| `dorado_basecalling_model`               | False    | `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` | Model used for basecalling.                                                                                             |
| `dorado_basecalling_extra_config`        | False    | -                                    | Extra configuration for Dorado basecalling.                                                                             |
| `dorado_basecalling_gpus`                | False    | `1`                                  | Number of GPUs to use for basecalling.                                                                                  |
| `skip_demultiplexingskip_demultiplexing` | False    | `false`                              | If `true`, the pipeline will not perform demultiplexing                                                                 |
| `dorado_demux_kit`                       | False    | `EXP-NBD196`                         | Kit used for demultiplexing.                                                                                            |
| `dorado_demux_both_ends`                 | False    | `false`                              | If `true`, the pipeline will demultiplex using barcodes from both sides (5' and 3').                                    |
| `dorado_demux_extra_config`              | False    | -                                    | Extra configuration for Dorado demultiplexing.                                                                          |
| `dorado_demux_cpus`                      | False    | `16`                                 | Number of CPUs to use for demultiplexing.                                                                               |
| `use_dorado_container`                   | False    | `true`                               | If `true`, the pipeline will use Dorado via container (~3.5GB download). If `false`, it will expect to find it locally. |

## Considerations

- It is possible to run the pipeline either in SLURM clusters using `-profile slurm`.
- Basecalling and demultiplexing are performed on separated steps to allow for a better control of the resources used by each process, and to prevent a whole basecalling redo in case of a failure during demultiplexing, wrong kit specified, etc.
- The basecalling process uses GPU, so make sure to have one available. If using SLURM, the job will be submitted with `--gres=gpu:X` option.
- Demultiplexing step won't use GPU, only CPU.
