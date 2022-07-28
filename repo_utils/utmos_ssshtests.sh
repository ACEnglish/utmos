
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
run test_conversion_blank $ut convert -B 2000 $INDIR/chunk0.vcf.gz $OD/chunk0.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/chunk0.jl) $(fn_md5 $INDIR/chunk0.jl)

run test_conversion_lowmem $ut convert -B 2000 --lowmem $INDIR/chunk1.vcf.gz $OD/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/chunk1.jl) $(fn_md5 $INDIR/chunk1.jl)

run test_conversion_compress $ut convert -B 2000 -c 1 $INDIR/chunk2.vcf $OD/chunk2.jl
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

run test_select_tiny $ut select -c 20 $INDIR/chunk_tiny.vcf -o $OD/select_tiny.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_tiny.txt) $(fn_md5 $ANSDIR/select_tiny.txt)

run test_select_one_af $ut select -c 0.005 --af -o $OD/select_one_af.txt $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_one_af.txt) $(fn_md5 $ANSDIR/select_one_af.txt)


run test_select_weights_subset $ut select --subset $INDIR/subset.txt -c 5 \
    --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_weights_subset.txt) $(fn_md5 $ANSDIR/select_weights_subset.txt)

run test_select_af_subset $ut select --subset $INDIR/subset.txt -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_subset.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af_subset.txt) $(fn_md5 $ANSDIR/select_af_subset.txt)

# ------------------------------------------------------------
#                                 select param checks
# ------------------------------------------------------------

run test_badfile $ut select doesntexist.txt 
assert_exit_code 1

run test_select_badmultih5 $ut select multi.hdf5 multi.hdf5
assert_exit_code 1

run test_select_badnoin $ut select
assert_exit_code 1

# ------------------------------------------------------------
#                                 select lowmem
# ------------------------------------------------------------
run test_select_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 $INDIR/chunk2.vcf
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)

run test_select_p1_lm $ut select --maxmem 1 --lowmem $INDIR/tiny.hdf5
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)

run test_select_p2_lm $ut select --maxmem 1 $INDIR/tiny.hdf5
assert_exit_code 0
assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)

run test_select_af_lm $ut select --maxmem 0 --lowmem $OD/tmp.hdf5 -c 20 --af \
    -o $OD/select_af.txt $INDIR/chunk0.jl $INDIR/chunk1.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_af.txt) $(fn_md5 $ANSDIR/select_af.txt)

run test_select_big_lm $ut select --maxmem 1 --lowmem $OD/tmp.hdf5 --count 5 $INDIR/big.jl \
    -o $OD/select_big.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_big.txt) $(fn_md5 $ANSDIR/select_big.txt)

run test_select_h5append $ut select --maxmem 0 -c 5 --buffer 500 --af -o $OD/select_multi_lm.txt\
    --lowmem $OD/tmp.hdf5 $INDIR/chunk0.jl $INDIR/chunk2.jl
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_multi_lm.txt) $(fn_md5 $ANSDIR/select_multi_lm.txt)

run test_select_tinyweight $ut select --maxmem 0 --count 20 --lowmem $OD/tmp.hdf5 $INDIR/chunk_tiny.vcf \
    -o $OD/select_tiny_weights.txt --weights $INDIR/weights.txt
assert_exit_code 0
assert_equal $(fn_md5 $OD/select_tiny_weights.txt) $(fn_md5 $ANSDIR/select_tiny_weights.txt)


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
