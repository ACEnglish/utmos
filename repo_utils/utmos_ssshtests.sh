
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest

# Work inside of the repo folder
cd "$( dirname "${BASH_SOURCE[0]}" )"/../
INDIR=repo_utils/test_files
ANSDIR=repo_utils/answer_key
OD=test_results
COVERAGE_RCFILE=.coveragerc

# Reset test results
rm -rf $OD
mkdir -p $OD

ut="coverage run --concurrency=multiprocessing -p -m utmos.__main__"
# ------------------------------------------------------------
#                                 test helpers
# ------------------------------------------------------------
fn_md5() {
    fn=$1
    # simple md5sum checking
    md5sum $fn | cut -f1 -d\  
}

df_check() {
    # check if joblib saved pandas dataframes are equivalent
    test_name=$1
    base_df=$2
    comp_df=$3
    run $test_name python3 -c """
import joblib;
a = joblib.load(\"$base_df\")
b = joblib.load(\"$comp_df\")
assert a.equals(b), \"$base_df != $comp_df\";
"""
}

# ------------------------------------------------------------
#                                 entry help
# ------------------------------------------------------------
run test_help $ut
assert_exit_code 0
assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)

run test_version $ut version
assert_exit_code 0


# ------------------------------------------------------------
#                                 conversion
# ------------------------------------------------------------
run test_conversion_blank $ut convert $INDIR/chunk0.vcf.gz $OD/chunk0.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/chunk0.jl) $(fn_md5 $INDIR/chunk0.jl)

run test_conversion_lowmem $ut convert --lowmem $INDIR/chunk1.vcf.gz $OD/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/chunk1.jl) $(fn_md5 $INDIR/chunk1.jl)

run test_conversion_compress $ut convert -c 1 $INDIR/chunk2.vcf $OD/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/chunk2.jl) $(fn_md5 $INDIR/chunk2.jl)

# ------------------------------------------------------------
#                                 select
# ------------------------------------------------------------
run test_select $ut select $INDIR/chunk2.vcf
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)

run test_select_intcnt $ut select --count 10 $INDIR/chunk1.jl 
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_intcnt.txt)

run test_select_floatcnt $ut select --count 0.01 $INDIR/chunk2.jl 
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_floatcnt.txt)

run test_select_fileout $ut select $INDIR/chunk1.vcf.gz -o $OD/select_fileout.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_fileout.txt) $(fn_md5 $ANSDIR/select_fileout.txt)

run test_select_multi $ut select -o $OD/select_multi.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.vcf
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multi.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_multi2 $ut select -o $OD/select_multi2.txt $INDIR/chunk0.jl $INDIR/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multi2.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_multimix $ut select -o $OD/select_multimix.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multimix.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_exclude $ut select -c 20 --exclude NA21117 $INDIR/chunk[01].jl -o $OD/select_exclude.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_exclude.txt) $(fn_md5 $ANSDIR/select_exclude.txt)

#run test_select_include $ut select -c 20 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include.txt) $(fn_md5 $ANSDIR/select_include.txt)

run test_select_weights $ut select -c 20 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weights.txt) $(fn_md5 $ANSDIR/select_weights.txt)

run test_select_af $ut select -c 20 --af -o $OD/select_af.txt $INDIR/chunk0.jl $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af.txt) $(fn_md5 $ANSDIR/select_af.txt)

run test_select_weightsaf $ut select -c 5 --af --weights $INDIR/weights.txt \
    -o $OD/select_weightsaf.txt $INDIR/chunk0.jl $INDIR/chunk1.jl  
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weightsaf.txt) $(fn_md5 $ANSDIR/select_weightsaf.txt)

#run test_select_exin $ut select -o $OD/select_exin.txt -c 18 --exclude NA21117 \
#    --include HG00096 --exclude HG02332,HG03097 $INDIR/chunk[01].jl --weights $INDIR/weights.txt --af 
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_exin.txt) $(fn_md5 $ANSDIR/select_exin.txt)

#run test_select_exinfile $ut select -o $OD/select_exinfile.txt -c 18 --exclude NA21117 --include HG00096 \
#    --exclude $INDIR/exclude.txt $INDIR/chunk[01].jl --weights $INDIR/weights.txt --af 
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_exinfile.txt) $(fn_md5 $ANSDIR/select_exin.txt)

run test_select_tiny $ut select -c 20 $INDIR/chunk_tiny.vcf -o $OD/select_tiny.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_tiny.txt) $(fn_md5 $ANSDIR/select_tiny.txt)

run test_select_one_af $ut select -c 0.005 --af -o $OD/select_one_af.txt $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_one_af.txt) $(fn_md5 $ANSDIR/select_one_af.txt)

run test_select_badsamps $ut select $INDIR/chunk_tiny.vcf $INDIR/chunk0.jl
assert_exit_code 1

run test_bad_filetype $ut select doesntexist.txt 
assert_exit_code 1

#run test_select_include_topN $ut select --mode topN -c 5 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include_topN.txt) $(fn_md5 $ANSDIR/select_include_topN.txt)

#run test_select_weights_topN $ut select --mode topN -c 5 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_weights_topN.txt) $(fn_md5 $ANSDIR/select_weights_topN.txt)

#run test_select_af_topN $ut select --mode topN -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_af_topN.txt) $(fn_md5 $ANSDIR/select_af_topN.txt)

#run test_select_random $ut select --mode random -c 5 -o $OD/select_random.txt $INDIR/chunk0.jl
#assert_exit_code 0

#run test_select_include_subset $ut select --subset $INDIR/subset.txt -c 5 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include_subset.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include_subset.txt) $(fn_md5 $ANSDIR/select_include_subset.txt)

run test_select_weights_subset $ut select --subset $INDIR/subset.txt -c 5 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weights_subset.txt) $(fn_md5 $ANSDIR/select_weights_subset.txt)

run test_select_af_subset $ut select --subset $INDIR/subset.txt -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af_subset.txt) $(fn_md5 $ANSDIR/select_af_subset.txt)

# ------------------------------------------------------------
#                                 select lowmem
#   I can shorten this later, but to start I want to test EVERYTHING
# ------------------------------------------------------------
run test_select_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 $INDIR/chunk2.vcf
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)

run test_select_intcnt_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 --count 10 $INDIR/chunk1.jl 
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_intcnt.txt)

run test_select_floatcnt_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 --count 0.01 $INDIR/chunk2.jl 
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_floatcnt.txt)

run test_select_fileout_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 $INDIR/chunk1.vcf.gz -o $OD/select_fileout.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_fileout.txt) $(fn_md5 $ANSDIR/select_fileout.txt)

run test_select_multi_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -o $OD/select_multi.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.vcf
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multi.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_multi2_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -o $OD/select_multi2.txt $INDIR/chunk0.jl $INDIR/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multi2.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_multimix_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -o $OD/select_multimix.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multimix.txt) $(fn_md5 $ANSDIR/select_multi.txt)

run test_select_exclude_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 20 --exclude NA21117 $INDIR/chunk[01].jl -o $OD/select_exclude.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_exclude.txt) $(fn_md5 $ANSDIR/select_exclude.txt)

#run test_select_include_lm $ut select --lowmem $OD/tmp.hdf5 -c 20 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include.txt) $(fn_md5 $ANSDIR/select_include.txt)

run test_select_weights_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 20 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weights.txt) $(fn_md5 $ANSDIR/select_weights.txt)

run test_select_af_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 20 --af -o $OD/select_af.txt $INDIR/chunk0.jl $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af.txt) $(fn_md5 $ANSDIR/select_af.txt)

run test_select_weightsaf_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 5 --af --weights $INDIR/weights.txt \
    -o $OD/select_weightsaf.txt $INDIR/chunk0.jl $INDIR/chunk1.jl  
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weightsaf.txt) $(fn_md5 $ANSDIR/select_weightsaf.txt)

#run test_select_exin_lm $ut selec --lowmem $OD/tmp.hdf5t -o $OD/select_exin.txt -c 18 --exclude NA21117 \
#    --include HG00096 --exclude HG02332,HG03097 $INDIR/chunk[01].jl --weights $INDIR/weights.txt --af 
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_exin.txt) $(fn_md5 $ANSDIR/select_exin.txt)

#run test_select_exinfile_lm $ut select --lowmem $OD/tmp.hdf5 -o $OD/select_exinfile.txt -c 18 --exclude NA21117 --include HG00096 \
#    --exclude $INDIR/exclude.txt $INDIR/chunk[01].jl --weights $INDIR/weights.txt --af 
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_exinfile.txt) $(fn_md5 $ANSDIR/select_exin.txt)

run test_select_tiny_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 20 $INDIR/chunk_tiny.vcf -o $OD/select_tiny.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_tiny.txt) $(fn_md5 $ANSDIR/select_tiny.txt)

run test_select_one_af_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 0.005 --af -o $OD/select_one_af.txt $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_one_af.txt) $(fn_md5 $ANSDIR/select_one_af.txt)

run test_select_badsamps_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 $INDIR/chunk_tiny.vcf $INDIR/chunk0.jl
assert_exit_code 1

run test_bad_filetype_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 doesntexist.txt 
assert_exit_code 1

#run test_select_include_topN_lm $ut select --lowmem $OD/tmp.hdf5 --mode topN -c 5 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include_topN.txt) $(fn_md5 $ANSDIR/select_include_topN.txt)

#run test_select_weights_topN_lm $ut select --lowmem $OD/tmp.hdf5 --mode topN -c 5 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_weights_topN.txt) $(fn_md5 $ANSDIR/select_weights_topN.txt)

#run test_select_af_topN_lm $ut select --lowmem $OD/tmp.hdf5 --mode topN -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_topN.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_af_topN.txt) $(fn_md5 $ANSDIR/select_af_topN.txt)

#run test_select_random_lm $ut select --lowmem $OD/tmp.hdf5 --mode random -c 5 -o $OD/select_random.txt $INDIR/chunk0.jl
#assert_exit_code 0

#run test_select_include_subset_lm $ut select --lowmem $OD/tmp.hdf5 --subset $INDIR/subset.txt -c 5 --include HG00096 $INDIR/chunk[01].jl -o $OD/select_include_subset.txt
#assert_exit_code 0
#assert_equal $(fn_md5 $OD/select_include_subset.txt) $(fn_md5 $ANSDIR/select_include_subset.txt)

run test_select_weights_subset_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 --subset $INDIR/subset.txt -c 5 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weights_subset.txt) $(fn_md5 $ANSDIR/select_weights_subset.txt)

run test_select_af_subset_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 --subset $INDIR/subset.txt -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af_subset.txt) $(fn_md5 $ANSDIR/select_af_subset.txt)


# ------------------------------------------------------------
#                                 coverage.py
# ------------------------------------------------------------

printf "\n${BOLD}generating test coverage reports${NC}\n"
coverage combine
coverage report --include=utmos/*
coverage html --include=utmos/* -d $OD/htmlcov/
coverage json --include=utmos/* -o $OD/coverage.json
run test_coverage python3 repo_utils/coverage_maker.py $OD/coverage.json
assert_exit_code 0
