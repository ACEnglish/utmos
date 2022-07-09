# UTMOS

A reimplementation of [SVCollector](https://github.com/fritzsedlazeck/SVCollector)

SVCollector is a tool for solving the maximum-coverage problem for sample selection. Utmos is a python port of that code
that leverages scikit-allel, numpy, and joblib. Utmos is designed for extremely large cohorts by allowing subsets of
variants to be parsed and stored as small(ish) intermediate files before concatenating during the selection step.

## Install

After cloning the repository, 
```bash
cd utmos/
python3 -m pip install . 
```

## Quick Start

Simplest use case:

```bash
utmos select input1.vcf
```

Conversion of a vcf to joblib

```bash
utmos convert input1.vcf input1.jl
```

See `--help` of commands for more details.

## utmos convert

Pulls information from `vcf[.gz]` into numpy arrays and saves using joblib.

usage: convert [-h] [-c COMPRESS] in_file out_file

positional arguments:
  in_file               Input VCF
  out_file              Output joblib

optional arguments:
  -h, --help            show this help message and exit
  -c COMPRESS, --compress COMPRESS
                        joblib compress level 1-9 (5)
```

Future features
* `--af` atempts to pull AF, calculates if AF not available

## utmos select

Select samples for validation and resequencing

```
usage: select [-h] [-o OUT] [-k KEEP] [--safe] in_files [in_files ...]

positional arguments:
  in_files              Input VCF or jl files

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     Output file (stdout)
  -k KEEP, --keep KEEP  Number of samples to select as a percent if <1 or
                        count if >=1 (0.02)
  --safe                Ensure input files have same sample names
```

Future features:
* `--mode` : greedy (default) random, topN 
* `--af` : take allele frequency into account
* `--include` : file (or comma-separated list) of samples to force inclusion
* `--exclude` : file (or comma-separated list) of samples to exclude from inclusion
* `--weights` : file of samples and their weights


## Performace metrics
Utmos is upto 4 times faster than SVCollector but uses more memory

--mode [greedy, topN, random]

## utmos plot

Future feature. Will make plots for the `select` output. Will have to figure out `--meta`
