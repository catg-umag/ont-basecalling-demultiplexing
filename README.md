# ONT Basecalling / Demux Pipeline

Small pipeline to perform basecalling and (optional) demultiplexing of ONT data, collect QC metrics and generate a MultiQC report.
Uses Dorado for basecalling and demultiplexing.

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.04)
- [Apptainer](https://apptainer.org/) / Singularity
- Dorado (0.7.3 tested). It can be used via container, or installed locally from https://github.com/nanoporetech/dorado.

## Quick Start
1. Clone this repository:
	```bash
	git clone https://github.com/your-repo/ont-basecalling-demultiplexing.git
	```
2. Demultiplexing setup (optional):
	- If demultiplexing is needed, create a samples.csv file containing at least the barcode and sample columns.
	- Ensure the barcode column includes the barcode identifier (e.g., barcode01), and the sample column lists the sample name, which will be used in reports and as the FASTQ filename.
3. Configure parameters:
	- Copy the example parameters file:
		```bash
		cp params.example.yml my_params.yml
		```
	- Modify my_params.yml according to your needs. Ensure that the sample_data parameter points to your samples.csv file if you are demultiplexing.
4. Run the pipeline:
	```bash
	nextflow run ont-basecalling-demultiplexing/ -profile apptainer -params-file my_params.yml
	```

## Pipeline Parameters

| Parameter                         | Required | Default                            | Description                                                                                     |
| --------------------------------- | -------- | ---------------------------------- | ----------------------------------------------------------------------------------------------- |
| `experiment_name`                 | No       | -                                  | Name of the experiment, used for final reports (title and filename).                            |
| `data_dir`                        | Yes      | -                                  | Path to the directory containing POD5 files.                                                    |
| `sample_data`                     | Yes      | `input/samples.csv`                | Path to the CSV file containing the sample data (required if demultiplexing).                   |
| `output_dir`                      | No       | `demultiplex_results`              | Directory for saving results.                                                                   |
| `fastq_output`                    | No       | `true`                             | Generates FASTQ files if `true`; otherwise, generates UBAM files.                               |
| `qscore_filter`                   | No       | `10`                               | Minimum QScore threshold for "pass" data, used in demultiplexing.                               |
| `dorado_basecalling_model`        | No       | `sup`                              | Model used for basecalling.                                                                     |
| `dorado_basecalling_extra_config` | No       | -                                  | Additional configuration options for Dorado basecalling.                                        |
| `dorado_basecalling_gpus`         | No       | `1`                                | Number of GPUs to allocate for basecalling.                                                     |
| `skip_demultiplexing`             | No       | `false`                            | Skips demultiplexing if `true`.                                                                 |
| `dorado_demux_kit`                | No       | `EXP-NBD196`                       | Kit identifier used for demultiplexing.                                                         |
| `dorado_demux_both_ends`          | No       | `false`                            | Demultiplexes using barcodes on both ends (5' and 3') if `true`.                                |
| `dorado_demux_extra_config`       | No       | -                                  | Additional configuration options for Dorado demultiplexing.                                     |
| `use_dorado_container`            | No       | `true`                             | Uses Dorado via container if `true`; expects a local installation if `false`.                   |
| `qc_tools`                        | No       | `['fastqc', 'nanoq', 'toulligqc']` | Specifies which QC tools to run. Options: 'nanoq', 'nanoplot', 'fastqc', 'toulligqc', 'pycoqc'. |

## Considerations

- The pipeline is compatible with SLURM clusters; use `-profile slurm`.
- GPU resources are required for basecalling. On SLURM, ensure jobs request GPUs with the `--gres=gpu:X` option.
