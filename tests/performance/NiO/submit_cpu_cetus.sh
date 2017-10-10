script=qmcpack-cpu-cetus.sub
nthreads=16
appendix=n256-t${nthreads}-BGQ
subjob=qsub

for size in 8 16 32 64
do
  source_folder=dmc-S$size-cpu
  folder=$source_folder-${appendix}
  echo $folder
  mkdir $folder
  cp sample/$source_folder/* $folder/
  cp $script $folder/
  cd $folder
  sed -i "s/XX/$size/" $script
  $subjob $script
  cd ..
done
