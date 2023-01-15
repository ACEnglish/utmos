
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest
STOP_ON_FAIL=1
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

jl_check() {
    # check if joblib saved pandas dataframes are equivalent
    test_name=$1
    base_df=$2
    comp_df=$3
    run $test_name python3 -c """
import joblib;
a = joblib.load(\"$base_df\")
b = joblib.load(\"$comp_df\")
assert a['stats'] == b['stats'], \"unequal stats - $base_df != $comp_df\";
for key in [\"GT\", \"AF\", \"samples\"]:
    assert (a[key] == b[key]).all, f\"unequal {key} - $base_df != $comp_df\";
"""
    assert_exit_code 0
}

# ------------------------------------------------------------
#                                 entry help
# ------------------------------------------------------------
run test_help $ut
if [ $test_help ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)
fi

run test_version $ut version
if [ $test_help ]; then
    assert_exit_code 0
fi


# ------------------------------------------------------------
#                                 conversion
# ------------------------------------------------------------
run test_conversion_blank $ut convert -B 2000 $INDIR/chunk0.vcf.gz $OD/chunk0.jl
if [ $test_conversion_blank ]; then
    assert_exit_code 0
    jl_check test_conversion_blank $OD/chunk0.jl $INDIR/chunk0.jl
fi

run test_conversion_lowmem $ut convert -B 2000 --lowmem $INDIR/chunk1.vcf.gz $OD/chunk1.jl
if [ $test_conversion_lowmem ]; then
    assert_exit_code 0
    jl_check test_conversion_lowmem $OD/chunk1.jl $INDIR/chunk1.jl
fi

run test_conversion_compress $ut convert -B 2000 -c 1 $INDIR/chunk2.vcf $OD/chunk2.jl
if [ $test_conversion_compress ]; then
    assert_exit_code 0
    jl_check test_conversion_compress $OD/chunk2.jl $INDIR/chunk2.jl
fi

# ------------------------------------------------------------
#                                 select
# ------------------------------------------------------------
run test_select $ut select $INDIR/chunk2.vcf
if [ $test_select ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)
fi

run test_select_intcnt $ut select --count 10 $INDIR/chunk1.jl 
if [ $test_select_intcnt ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_intcnt.txt)
fi

run test_select_floatcnt $ut select --count 0.01 $INDIR/chunk2.jl 
if [ $test_select_floatcnt ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_floatcnt.txt)
fi

run test_select_fileout $ut select $INDIR/chunk1.vcf.gz -o $OD/select_fileout.txt
if [ $test_select_fileout ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_fileout.txt) $(fn_md5 $ANSDIR/select_fileout.txt)
fi

run test_select_multi $ut select -o $OD/select_multi.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.vcf
if [ $test_select_multi ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_multi.txt) $(fn_md5 $ANSDIR/select_multi.txt)
fi

run test_select_multi2 $ut select -o $OD/select_multi2.txt $INDIR/chunk0.jl $INDIR/chunk2.jl
if [ $test_select_multi2 ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_multi2.txt) $(fn_md5 $ANSDIR/select_multi.txt)
fi

run test_select_multimix $ut select -o $OD/select_multimix.txt $INDIR/chunk0.vcf.gz $INDIR/chunk2.jl
if [ $test_select_multimix ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_multimix.txt) $(fn_md5 $ANSDIR/select_multi.txt)
fi

run test_select_exclude $ut select -c 20 --exclude NA21117 $INDIR/chunk[01].jl -o $OD/select_exclude.txt
if [ $test_select_exclude ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_exclude.txt) $(fn_md5 $ANSDIR/select_exclude.txt)
fi

run test_select_weights $ut select -c 20 --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights.txt
if [ $test_select_weights ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_weights.txt) $(fn_md5 $ANSDIR/select_weights.txt)
fi

run test_select_af $ut select -c 20 --af -o $OD/select_af.txt $INDIR/chunk0.jl $INDIR/chunk1.jl
if [ $test_select_af ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_af.txt) $(fn_md5 $ANSDIR/select_af.txt)
fi

run test_select_weightsaf $ut select -c 5 --af --weights $INDIR/weights.txt \
    -o $OD/select_weightsaf.txt $INDIR/chunk0.jl $INDIR/chunk1.jl  
if [ $test_select_weightsaf ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_weightsaf.txt) $(fn_md5 $ANSDIR/select_weightsaf.txt)
fi

run test_select_tiny $ut select -c 20 $INDIR/chunk_tiny.vcf -o $OD/select_tiny.txt
if [ $test_select_tiny ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_tiny.txt) $(fn_md5 $ANSDIR/select_tiny.txt)
fi

run test_select_one_af $ut select -c 0.005 --af -o $OD/select_one_af.txt $INDIR/chunk1.jl
if [ $test_select_one_af ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_one_af.txt) $(fn_md5 $ANSDIR/select_one_af.txt)
fi


run test_select_weights_subset $ut select --subset $INDIR/subset.txt -c 5 \
    --weights $INDIR/weights.txt $INDIR/chunk0.jl -o $OD/select_weights_subset.txt
if [ $test_select_weights_subset ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_weights_subset.txt) $(fn_md5 $ANSDIR/select_weights_subset.txt)
fi

run test_select_af_subset $ut select --subset $INDIR/subset.txt -c 5 --af  $INDIR/chunk0.jl -o $OD/select_af_subset.txt
if [ $test_select_af_subset ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $OD/select_af_subset.txt) $(fn_md5 $ANSDIR/select_af_subset.txt)
fi

# ------------------------------------------------------------
#                                 select param checks
# ------------------------------------------------------------

run test_badfile $ut select doesntexist.txt 
if [ $test_badfile ]; then
    assert_exit_code 1
fi

run test_select_badmultih5 $ut select multi.hdf5 multi.hdf5
if [ $test_select_badmultih5 ]; then
    assert_exit_code 1
fi

run test_select_badnoin $ut select
if [ $test_select_badnoin ]; then
    assert_exit_code 1
fi

# ------------------------------------------------------------
#                                 select lowmem
# ------------------------------------------------------------

run test_select_lm $ut select --maxmem 0 --lowmem $OD/tiny.hdf5 $INDIR/chunk2.vcf
if [ $test_select_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)
    assert_equal $(fn_md5 $OD/tiny.hdf5) $(fn_md5 $INDIR/tiny.hdf5)
fi

run test_select_p1_lm $ut select --maxmem 1 --lowmem $INDIR/tiny.hdf5
if [ $test_select_p1_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)
fi

run test_select_p2_lm $ut select --maxmem 1 $INDIR/tiny.hdf5
if [ $test_select_p2_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_first.txt)
fi

run test_select_af_lm $ut select --maxmem 0 -c 20 --af --lowmem $OD/tiny.af.hdf5 $INDIR/chunk0.vcf $INDIR/chunk1.vcf.gz
if [ $test_select_af_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_af_h5.txt)
    assert_equal $(fn_md5 $OD/tiny.af.hdf5) $(fn_md5 $INDIR/tiny.af.hdf5)
fi

run test_select_af_p1_lm $ut select --maxmem 1 -c 20 --lowmem $INDIR/tiny.af.hdf5
if [ $test_select_af_p1_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_af_h5.txt)
fi

run test_select_af_p2_lm $ut select --af --maxmem 1 -c 20 $INDIR/tiny.af.hdf5
if [ $test_select_af_p2_lm ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDOUT_FILE) $(fn_md5 $ANSDIR/select_af_h5.txt)
fi


# ------------------------------------------------------------
#                                 coverage.py
# ------------------------------------------------------------
# Don't generate coverage when doing subset of tests
if [ -z "$1" ]; then
    printf "\n${BOLD}generating test coverage reports${NC}\n"
    coverage combine
    coverage report --include=utmos/*
    coverage html --include=utmos/* -d $OD/htmlcov/
    coverage json --include=utmos/* -o $OD/coverage.json
    python3 repo_utils/coverage_maker.py $OD/coverage.json
fi
rm -f .coverage.*
