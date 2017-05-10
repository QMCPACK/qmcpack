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

SRC_ROOT=`pwd`/..

#
#  Run coverage on the base
#

echo "Running base coverage"

ctest -L coverage

base_dir="cov_base"
if [ ! -d $base_dir ]; then
  mkdir $base_dir
fi

cd $base_dir
find `pwd`/.. -name \*.gcda | xargs -i gcov -b -p {}

# Filter out unwanted files
python ${SRC_ROOT}/tests/coverage/compare_gcov.py -a process --base-dir .


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

unit_dir="cov_unit"

if [ ! -d $unit_dir ]; then
  mkdir $unit_dir
fi

cd $unit_dir
find `pwd`/.. -name \*.gcda | xargs -i gcov -b -p {}

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
