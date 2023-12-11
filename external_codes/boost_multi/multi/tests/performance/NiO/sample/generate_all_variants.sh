#!/bin/bash

# this script generates all the input files for GPU and CPU-J3 tests from CPU test input files.

for name in dmc-*-cpu
do

# generate GPU input files
  new_folder=`echo $name | sed "s/cpu/gpu/"`
  natom=`echo $name | sed "s/^.*a//" | sed "s/\-e.*$//"`
  size=$((natom/4))
  mkdir $new_folder
  ref_file=$name/NiO-fcc-S$size-dmc.xml
  new_file=$new_folder/NiO-fcc-S$size-dmc.xml
  ../../adjust_qmcpack_input.py -w 32 $ref_file -o $new_file

  # bring back the commended lines
  grep '^    <!' $ref_file > commented_lines
  sed -i "0,/^    $/{/^    $/d;}" $new_file
  sed -i "/^    $/r commented_lines" $new_file
  sed -i "/^    $/d" $new_file

# generate GPU input files with J3
  if [ -e $name/J123.xml ]
  then
    new_folder=$name-J3
    natom=`echo $name | sed "s/^.*a//" | sed "s/\-e.*$//"`
    size=$((natom/4))
    mkdir $new_folder
    ref_file=$name/NiO-fcc-S$size-dmc.xml
    new_file=$new_folder/NiO-fcc-S$size-dmc.xml
    ../../adjust_qmcpack_input.py -j $name/J123.xml $ref_file -o $new_file

    # bring back the commended lines
    grep '^    <!' $ref_file > commented_lines
    sed -i "0,/^    $/{/^    $/d;}" $new_file
    sed -i "/^    $/r commented_lines" $new_file
    sed -i "/^    $/d" $new_file
  fi

done

rm commented_lines
