#!/bin/bash

if [ ! -d src ]; then
  echo Not in QMCPACK repo top directory. No action taken.
  exit
fi

match="DUMMY"
for folder in `ls src`
do
  if [ -d src/$folder ]; then
    echo processing folder name : $folder
    match=${match}"\|"${folder}
  fi
  for sub_folder in `ls src/$folder`
  do
    if [ -d src/$folder/$sub_folder ]; then
      echo processing sub folder name : $sub_folder
      match=${match}"\|"${sub_folder}
    fi
  done
done
echo $match

grep -R "<\($match\)\/.*>" src

for extension in h hpp cpp c cu cpp.in h.in
do
  echo processing extension $extension
  find src -name "*\.${extension}" -exec sed -i "/include/s/<\(${match}\)\/\(.*\)>/\"\1\/\2\"/" {} \;
done

for extension in cpp c cu
do
  for name in `find src -name "*\.${extension}"`
  do
    file_name_without_extension=`basename $name | sed "s/\.${extension}$"//""`
    full_path_file_header_h=`echo $name | sed "s/\.${extension}$"/.h/""`
    full_path_file_header_hpp=`echo $name | sed "s/\.${extension}$"/.hpp/""`
    if [ -f $full_path_file_header_h ] || [ -f $full_path_file_header_hpp ]; then
      echo $file_name_without_extension $full_path_file_header_h $full_path_file_header_hpp
      sed -i "/include/s/\"${file_name_without_extension}.\(h\|hpp\)\"/\"${file_name_without_extension}.\1\"/" $name
      sed -i "/include/s/\".*\/${file_name_without_extension}.\(h\|hpp\)\"/\"${file_name_without_extension}.\1\"/" $name
      sed -i "/include/s/<${file_name_without_extension}.\(h\|hpp\)>/\"${file_name_without_extension}.\1\"/" $name
      sed -i "/include/s/<.*\/${file_name_without_extension}.\(h\|hpp\)>/\"${file_name_without_extension}.\1\"/" $name
    fi
    #sed -i "s/#include\"/#include \"/" $name
    #sed -i "s/#include</#include </" $name
  done
done
