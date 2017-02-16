script=qmcpack.sub
nthreads=16

for size in 8 16 32 64
do
  source_folder=dmc-S$size-cpu
  folder=$source_folder-n256-t${nthreads}-BGQ
  echo $folder
  mkdir $folder
  cp sample/$source_folder/* $folder/
  cp $script $folder/
  cd $folder
  sed -i "s/XX/$size/" $script
  qsub $script
  cd ..
done
