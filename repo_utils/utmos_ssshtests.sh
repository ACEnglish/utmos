
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest

# Work inside of the repo folder
cd "$( dirname "${BASH_SOURCE[0]}" )"/../
INDIR=repo_utils/test_files
ANSDIR=$INDIR/answer_key
OD=test_results
COVERAGE_RCFILE=.coveragerc

# Reset test results
rm -rf $OD
mkdir -p $OD

ut="coverage run --concurrency=multiprocessing -p -m ut.__main__"
# ------------------------------------------------------------
#                                 test helpers
# ------------------------------------------------------------
fn_md5() {
    fn=$1
    # simple md5sum checking
    md5sum $fn | cut -f1 -d\  
}

# ------------------------------------------------------------
#                                 conversion
# ------------------------------------------------------------
run test_conversion_1 $ut
assert_exit_code 0
assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)
# Simple md5sum works

# ------------------------------------------------------------
#                                 selection
# ------------------------------------------------------------
run test_selection $ut
# --safe (fail) and --keep as float and as int


# ------------------------------------------------------------
#                                 plots
# ------------------------------------------------------------



