# UTMOS
![coverage](imgs/coverage.svg)
![pylint](imgs/pylint.svg)

Maximum-coverage algorithm to select samples for validation and resequencing.
This is a reimplementation of [SVCollector](https://github.com/fritzsedlazeck/SVCollector)

Utmos is a tool for solving the maximum-coverage problem for sample selection. Utmos leverages scikit-allele, joblib,
numpy, and hdf5 to allow inputs from extremely large cohorts. Subsets of variants can be converted into smaller
processing-ready files, which allows the data intake step to become massively parallel. Then, during selection, Utmos
concatenates the subsets and allows control over memory consumption for easier analysis.

## Install

Download a [release](https://github.com/ACEnglish/utmos/releases) and run
```
python3 -m pip install Utmos-<version>.tar.gz
```

Alternatively, build from the repository, 
```bash
git clone https://github.com/ACEnglish/utmos.git
cd utmos/
python3 -m pip install . 
```

## Quick Start

```bash
utmos select input1.vcf > sample_report.txt
```

See `--help` of commands for more details.

## Output format

Tab-delimited file of:

<table><tr><th>Column</th><th>Description</th>
<tr><td>sample</td><td>Sample name</td></tr>
<tr><td>var_count</td><td>Number of variants in the sample</td></tr>
<tr><td>new_count</td><td>Number of new variants contributed by the sample</td></tr>
<tr><td>tot_captured</td><td>Running total of number of variants captured by all samples upto this point</td></tr>
<tr><td>pct_captured</td><td>Percent of all variants captured by samples upto this point</td></tr>
</table>

## utmos convert

Pulls genotype presence/absence information from `vcf[.gz]` into numpy arrays, calcluates allele frequencies
and saves using joblib. This step is optional, but makes it easier to convert multiple VCFs at once
(with separate jobs). The genotype array is shrunk using `numpy.packbits` and joblib compression can reduce file size
even further.  
As a test, the genotype-only chr22 snps from the 1kgp (2,504 samples x 1,103,547 variants)
<table><tr><th>File</th><th>Size (megabytes)</th><th>Fold Decrease</th>
<tr><td>VCF</td><td>198</td><td>.</td></tr>
<tr><td>without packbits</td><td>64</td><td>3.09</td></tr>
<tr><td>with packbits</td><td>29</td><td>6.82</td></tr>
<tr><td>--compress 9</td><td>26</td><td>7.62</td></tr>
<tr><td>both axis pack*</td><td>12</td><td>16.5</td></tr>
<tr><td>both axis pack* -c 9 </td><td>11</td><td>18</td></tr>
</table>

\* both axis pack is not yet implemented since the overhead it requires slows runtime a little bit. I'll implement it if there's any demand

```
usage: convert [-h] [--no-singleton] [--lowmem] [-B BUFFER] [-c COMPRESS] in_file out_file
```
### Arguments:
* `--no-singleton` exclude variants which are singletons
* `--lowmem` lowers memory usage by converting directly into an intermediate hdf5 file.
* `--buffer` is an integer pased to scikit-allele during conversion and represents how many variants are parsed at
  once. The default value is probably fine for many use-cases, but performance can be tested. Generally, VCFs with 
  extremely high sample counts (100k+) should have a smaller buffer (the default). VCFs without many samples or
  with high amounts of memory available can use a higher buffer in an attempt to speed up conversion.
* `--compress` compression level performed by joblib between 0 (lowest compression, fastest IO) and 10 (high
  compression, slower IO)
* `in_file` an input vcf[.gz]
* `out_file` an output .jl file

### Preparsing
VCFs can be pre-filtered and piped into convert e.g.:
```
bcftools view -i "FILTER == 'PASS'" input.vcf.gz | utmos convert /dev/stdin output.jl
```

## utmos select

Select samples for validation and resequencing.
Works by iteratively selecting samples based on a score and reporting how many variants the sample contains as well as
how many unseen variants are contributed by the sample. Scores by default are variant counts which can optionally be
weighed with `--weights` and/or `--af`

```
usage: select [-h] [-c COUNT] [-o OUT] [--debug] [--af] [--weights WEIGHTS]
	      [--subset SUBSET] [--exclude EXCLUDE] [--lowmem LOWMEM]
	      [--buffer BUFFER] [--maxmem MAXMEM]
	      [in_files ...]
```

### Basic arguments:
#### in_files
  one or more input files and can be a mix of vcf[.gz] or jl files from `utmos convert`. If there was a
  previous run of `utmos select --lowmem example.hdf5`, the a single in_file of the hdf5 file can be specified. This
  saves time by not needing to concatenate results between multiple runs of a dataset.

#### count
Sets how many samples are selected. Can be a count (>= 1) or a percent of samples if < [0, 1). Select all samples (or
until all variants have been used, whichever comes first) with -1

### Scoring arguments:
#### af
Reduce bias towards rare/private alleles by using the allele-frequency weighted matrix for scoring. So, instead of 0/1
for variant presence in each cell of the matrix, the value will be `presence * AF`

#### weight
A tab-delimited file of samples and their weights. Not every sample in the vcf needs to be given a score in the weight file.
Any sample without a provided weight is given a 1.

#### subset
A subset of samples to include in the selection. To help organize your metadata, multiple `--subset` arguments can be
provided and they will be concatenated. e.g. you have `subsetA.txt` and you would like to add sample `HG1234` to that
subset for an utmos run, simply specify `utmos select --subset subsetA.txt --subset HG1234 input.jl`

#### exclude
A list of samples to exclude from the analysis. Similar to `--subset`, you can provide multiple `--exclude` arguments.

### Memory arguments:
By default, utmos will hold all variants in memory. This consumes approximately

`(num_variants * num_samples * dtype_bytes) / 1e9` GB of memory. 

Where `dtype_bytes` is the size per-datapoint (typically 4 bytes).

For datasets with too much data to hold in memory, utmos allows a user to take advantage of hdf5 files and use disk
storage instead of memory. However, this comes at the expense of needing to create intermediate files and increases
IO, thus one should expect an increased runtime.

The `--lowmem` mode works by first concatenating all the input vcf/jl together into an hdf5 file. `--buffer`
rows will be held in-memory before writing/appending the data into the hdf5 file. Then, for each iteration, the dataset
is reduced by removing the used sample's column and all variants the sample contained from the rows before writing the 
reduced matrices to a temporary hdf5 in the `$TMPDIR`. If at any point through the iterations the dataset is estimated 
to have a total size below `--maxmem`, the full dataset is pulled into memory before continuing to the end of the
analysis.

#### lowmem
Name of the hdf5 file which is filled with the concatenated in_files.

#### buffer
Number of variants (rows) to buffer in memory during concatenation before writing to the hdf5 file.

#### maxmem
Maximum memory (in GB) available to utmos. For debugging, taking the 'shortcut' to pull the data into memory if/when
possible can be skipped by setting `--maxmem 0`. However, when actually performing summation, utmos assumes at least 1GB
of memory is available.

## Reusing the `--lowmem file.hdf5` with `select`
If you have a previous run in which you used `--lowmem` to make an hdf5 file, you can reuse it and save
the concatenation/conversion work. Simply provide a single `in_file` of `file.hdf5` or provide no `in_files` and the
parameter `--lowmem file.hdf5`. Note that if you create an hdf5 file with `select --af`, it will hold the AF weighted
matrix and can only be reused with `select --af`.

## Multi-allelic VCFs
When running `select --af`, variant positions are weighed by their allele frequency. For multi-allelic VCFs, the site is
weighed by the largest allele frequency observed. If this is not the desired behavior, split multi-allelics in the VCF 
with `bcftools norm -N -m-any`.

## Performace metrics
Running on a 2013 Mac Pro and using chr22 from 1kgp genotype  
`ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`

2,504 samples x 1,103,547 variants

Utmos runtime:
```
real	24m26.507s
user	11m37.160s
sys	5m12.713s
```

SVCollector runtime: (including 30 seconds to uncompress the VCF)
```
real	61m34.008s
user	49m28.622s
sys	1m24.693s
```

On newer hardware (Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz) and running the docker through singularity:
```
#Utmos
real	6m30.870s
user	5m53.305s
sys	0m36.403s

#SVCollector
real	7m44.110s
user	7m8.516s
sys	0m17.368s
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
