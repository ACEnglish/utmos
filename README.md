1) `utmos cvt --in in.vcf --out out.jl --af --compress 3`
Convert VCF to a joblib matrix of present or not 
--af
Attempts to pull 'AF' 
if not available
>>> ac = g.count_alleles()
>>> ac.to_frequencies()

--in (default stdin)
--out (hdf5 for now, I guess?)

2) utmos calc --out required.txt \*files [--safe]
Calculate the greedy on inputs

--vcf (default = False)
# parse args decides if we're doing vcfs (in-memory or not)
# this makes me belive we can't do hdf5 because of it is directly to a file...
If args.vcf:
	for i in args[:]:
		# param for in memory?
		utmos cvt i (tmpfile)
		args.inputs.append(result)

open each args.inputs and concat (tiny memory hopefully)

--safe
Ensure that the sample names are the same between each of the parts

--count
How many samples to pick (so we can eventually stop, if it is a float (pick a percent) if int, take the count
Interpret 1 as 1%

[maybe]
--af
Use AF for something? IDK
Use the AF data from above

--mode [greedy, topN, random]

3) `utmos plot output.txt`
make the curves with seaborn/pandas (easy enough)
do we want metadata to color it by... --meta [sample\tlabel]
and just throw it in a svg or png in that order

