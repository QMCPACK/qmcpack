script=qmcpack-gpu.sub
walker=128

for size in 8 16 32 64
do
  source_folder=dmc-S$size-gpu
  folder=$source_folder-n2-w$walker-K80
  echo $folder
  mkdir $folder
  cp sample/$source_folder/* $folder/
  cp $script $folder/
  cd $folder
  sed -i "s/XX/$size/" $script
  sed -i "/walkers/s/32/$walker/" *.xml
  qsub $script
  cd ..
done
