#!/bin/bash

# Assumes build has been build with ENABLE_GCOV and DEBUG
# Run from build directory

# Output:
#  cov_base - coverage for the base run(s)
#  cov_unit - coverage for the unit tests (the usual code coverage metric)
#  cov_diff - coverage for the unit tests relative to the base run

# Each directory contains the *.gcov files, one for each source file.
# If the 'gcovr' tool is present ( http://gcovr.com/ ), an HTML report is also
#  produced in cov_detail.html

CURRENT_DIR=`pwd`
#SRC_ROOT=`pwd`/..
# Need parent directory without putting a '..' in the path
SRC_ROOT=`dirname $CURRENT_DIR`

#
#  Run coverage on the base
#

echo "Running base coverage"

ctest -L coverage

raw_base_dir="cov_base_raw"
if [ ! -d $raw_base_dir ]; then
  mkdir $raw_base_dir
fi

base_dir="cov_base"
if [ ! -d $base_dir ]; then
  mkdir $base_dir
fi

cd $raw_base_dir
find `pwd`/.. -name \*.gcda | xargs -i gcov -b -p -l -s ${SRC_ROOT} {}

# Filter out unwanted files
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a process --base-dir .

cd ../$base_dir
# Combine different gcov files corresponding to the same source
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a merge --base-dir ../$raw_base_dir --prefix $SRC_ROOT


# If gcovr is present, create an html report
if hash gcovr 2>/dev/null; then
  echo "Generating HTML report for base coverage"
  # Arguments:
  # -k - keep the *.gcov files, don't delete them
  # -g - operation on existing *.gcov files (by default it would run gcov again)
  # --root - otherwise all the files are filtered out
  #  Run w/o the --html and --html-details options to get a text summary
  gcovr -k -g --root=$SRC_ROOT -o cov_detail.html --html --html-details .
fi

cd ..

${SRC_ROOT}/tests/coverage/clean_gcda.sh

#
#  Unit test coverage
#

echo "Running unit test coverage"
ctest -L unit

raw_unit_dir="cov_unit_raw"
if [ ! -d $raw_unit_dir ]; then
  mkdir $raw_unit_dir
fi

unit_dir="cov_unit"
if [ ! -d $unit_dir ]; then
  mkdir $unit_dir
fi

cd $raw_unit_dir
find `pwd`/.. -name \*.gcda | xargs -i gcov -b -p -l -s ${SRC_ROOT} {}

cd ../$unit_dir

# Combine different gcov files corresponding to the same source
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a merge --base-dir ../$raw_unit_dir --prefix $SRC_ROOT


# Filter out unwanted files
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a process --base-dir .

# If gcovr is present, create an html report
if hash gcovr 2>/dev/null; then
  echo "Creating HTML report for unit tests coverage"
  gcovr -k -g --root=$SRC_ROOT -o cov_detail.html --html --html-details .
fi

cd ..


#
# Compute diff
#

echo "Computing coverage diff"


diff_dir="cov_diff"

if [ ! -d $diff_dir ]; then
  mkdir $diff_dir
fi

# Compute diff
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a compare --base-dir $base_dir --unit-dir $unit_dir --output $diff_dir

cd $diff_dir

# If gcovr is present, create an html report
if hash gcovr 2>/dev/null; then
  echo "Creating HTML report for diff coverage"
  gcovr -k -g --root=$SRC_ROOT -o cov_detail.html --html --html-details .
fi

cd ..
