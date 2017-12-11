script=qmcpack-gpu-cooley.sub
walker=128
appendix=n2-w$walker-K80
subjob=qsub

for size in 8 16 32 64
do
  Natom=$((size * 4))
  Nelec=$((size * 48))
  source_folder=dmc-a$Natom-e$Nelec-gpu
  folder=$source_folder-${appendix}
  echo $folder
  mkdir $folder
  cp sample/$source_folder/* $folder/
  cp $script $folder/
  cd $folder
  sed -i "s/XX/$size/" $script
  sed -i "/walkers/s/32/$walker/" *.xml
  $subjob $script
  cd ..
done
