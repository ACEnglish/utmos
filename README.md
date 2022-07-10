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

```bash
utmos select input1.vcf
```

See `--help` of commands for more details.

## utmos convert

Pulls information from `vcf[.gz]` into numpy arrays and saves using joblib.

```
usage: convert [-h] [--lowmem] [-c COMPRESS] in_file out_file

positional arguments:
  in_file               Input VCF
  out_file              Output joblib

optional arguments:
  -h, --help            show this help message and exit
  --lowmem              Lower memory usage with hdf5 temporary files (False)
  -c COMPRESS, --compress COMPRESS
                        joblib compress level 1-9 (5)
```

Future features
* `--af` atempts to pull AF, calculates if AF not available

## utmos select

Select samples for validation and resequencing

```
usage: select [-h] [--lowmem] [-o OUT] [-c COUNT] [--safe] [--include INCLUDE]
              [--exclude EXCLUDE]
              in_files [in_files ...]

positional arguments:
  in_files              Input VCF or jl files

optional arguments:
  -h, --help            show this help message and exit
  --lowmem              Use temporary files to lower memory usage during vcf
                        conversion (False)
  -o OUT, --out OUT     Output file (stdout)
  -c COUNT, --count COUNT
                        Number of samples to select as a percent if <1 or
                        count if >=1 (0.02)
  --safe                Ensure input files have same sample names
  --include INCLUDE     Filename with or Comma-separated list of samples to
                        force selection
  --exclude EXCLUDE     Filename with or Comma-separated list of samples to
                        exclude selection
```

Future features:
* `--mode` : greedy (default), random, topN 
* `--af` : take allele frequency into account
* `--weights` : file of samples and their weights

## utmos plot

Future feature. Will make plots for the `select` output. Will have to figure out `--meta`

## Performace metrics
Using chr22 from 1kgp genotype
[link](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz)

2,504 samples x 1,103,547 variants

Utmos runtime:
```
real	31m32.103s
user	14m39.722s
sys	7m20.701s
```

SVCollector runtime: (including 30 seconds to uncompress the VCF)
```
real	61m34.008s
user	49m28.622s
sys	1m24.693s
```

## Dockerfile

A Dockerfile exists to build an image of utmos. To make a Docker image, clone the repository and run
```bash
docker build -t utmos .
```

You can then run utmos through docker using
```bash
docker run -v `pwd`:/data -it utmos
```
Where `pwd` can be whatever directory you'd like to mount in the docker to the path `/data/`, which is the working
directory for the utmos run. You can provide parameters directly to the entry point.
```bash
docker run -v `pwd`:/data -it utmos convert example.vcf.gz example.jl
```

