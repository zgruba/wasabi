# WASABI – Wasserstein-based Annotation of Spectra Against Base Index

Tool dedicated for spectrum annotation leveraging q-Wasserstein algorithm, accompanied by a database generation module based on IsoSpec with *wasabi\_db.hpp* flexible library capable of utilising the q-Wasserstein algorithm independently of IsoSpec, making it applicable to a broader range of datasets beyond mass spectrometry.

## Requirements
The project depends on several external libraries and frameworks. OpenMS is required for LC--MS data handling and feature processing, OpenMP is used to easily enable parallelism, IsoSpec provides isotopic distribution calculations, and Boost is utilised for general-purpose C++ utilities and data structures.

## Usage

This section describes the intended usage of the software package developed in this work, with particular emphasis on its command-line utilities and their integration within the LC–MS/MS feature annotation workflow.

### 1. CSV preparation for DeepLC

The repository provides a dedicated utility for preparing input files compatible with DeepLC. This tool parses FASTA file and generate all possible peptides with availble PTMs under the preselected constraints and formats it into a comma-separated values file suitable for retention time prediction.

The executable *prepare\_csv\_for\_deeplc* operates produces a standardised CSV representation required by DeepLC.

**Command-line usage:**

```
prepare_csv_for_deeplc <fasta> <output.csv> [options]


Options:
  -h [ --help ]              Show help message
  -t [ --threads ] arg (=1)  Number of threads (default: 1)
  -d [ --with-decoys ]       Include decoys (toggle)
```

**Output:**
The resulting CSV file contains peptide sequences and associated metadata formatted according to DeepLC specifications, and can be used directly as input for retention time prediction.

### 2. WASABI feature annotation

The core functionality of the package is provided by the *wasabi\_annotate* executable, which performs peptide annotation of LC-MS features. The tool integrates several processing stages, including in silico database generation from a FASTA file, computation of theoretical isotopic distributions using IsoSpec, and candidate searching based on the q-Wasserstein distance.

**Command-line usage:**
```
wasabi_annotate <features> <mzML> <fasta> <output> --multipliers <list> [options]


Options:
  -h [ --help ]                         Show help message
  -t [ --threads ] arg (=1)             Number of threads (default: 1)
  -m [ --multipliers ] arg              Multipliers (space-separated list, 
                                        required)
  --queue-threshold arg (=0.14999999999999999)
                                        Queue threshold (default: 0.15)
  --search-threshold arg (=0.14999999999999999)
                                        Search threshold (default: 0.15)
  --rt-threshold arg (=1750)            RT threshold (default: 1750)
  --reduced-decoys                      Use reduced decoys
  -d [ --with-decoys ]                  Include decoys
  --use-rt-filtering                    Enable RT filtering
  --rt-mapping arg                      RT mapping file (required if 
                                        --use-rt-filtering is set)

Positional arguments:
  --feature_path arg                    Feature file (required)
  --mzml_path arg                       mzML file (required)
  --fasta_path arg                      FASTA file (required)
  --output_path arg                     Output path (required)
```

**Output:**
For each processed LC--MS feature, the output directory contains CSV files listing candidate peptide annotations together with their associated scores, ranks, and relevant metadata.

## Citing

If you use tools from this package, please cite:

Kozaryna, Z. WASABI – Wasserstein-based Annotation of Spectra Against Base Index (2025)