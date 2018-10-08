MIN_SIZE=16
MAX_SIZE=16
DIS=20

for SIZE in `seq $MIN_SIZE $MAX_SIZE`
do

N_elec_half=$((4 * SIZE))
N_atom=$((3 * SIZE))

eff_out=coeff.out
echo -n > $eff_out
for cell in `seq $SIZE`
do
  for ib in `seq 4`
  do
    for copy in `seq $SIZE`
    do
      if [ $copy -eq $cell ] ; then
        cat band$ib >> $eff_out
      else
        cat band0   >> $eff_out
      fi
    done
  done
done

sed "/COEFF/r $eff_out" H2O.empty.wfs.xml | sed "/COEFF/d" | sed "s/N_ELEC_HALF/$N_elec_half/" > H2O.wfs.xml

up_pos_out=up_pos.out
down_pos_out=down_pos.out
atom_pos_out=atom_pos.out
label_out=label.out
echo -n > $up_pos_out
echo -n > $down_pos_out
echo -n > $atom_pos_out
echo -n > $label_out
for cell in `seq $SIZE`
do
  awk -v dis=$((DIS*(cell-1))) '{print "      ", $1+dis, $2, $3}' elec_up_pos >> $up_pos_out
  awk -v dis=$((DIS*(cell-1))) '{print "      ", $1+dis, $2, $3}' elec_down_pos >> $down_pos_out
  awk -v dis=$((DIS*(cell-1))) '{print "     ", $1+dis, $2, $3}' atom_pos >> $atom_pos_out
  echo "      O H H" >> $label_out
done

sed "/ELEC_UP_POS/r $up_pos_out" H2O.empty.ptcl.xml | sed "/ELEC_UP_POS/d" > H2O.ptcl.xml
sed "/ELEC_DOWN_POS/r $down_pos_out" H2O.ptcl.xml | sed "/ELEC_DOWN_POS/d" > H2O.temp.ptcl.xml
sed "/ATOM_POS/r $atom_pos_out" H2O.temp.ptcl.xml | sed "/ATOM_POS/d" > H2O.ptcl.xml
sed "/ATOM_LABEL/r $label_out" H2O.ptcl.xml | sed "/ATOM_LABEL/d" > H2O.temp.ptcl.xml

sed -i "s/N_ELEC_HALF/$N_elec_half/" H2O.temp.ptcl.xml
sed -i "s/N_ATOM/$N_atom/" H2O.temp.ptcl.xml
mv H2O.temp.ptcl.xml H2O.ptcl.xml

rm $up_pos_out $down_pos_out $atom_pos_out $label_out $eff_out

folder=S$SIZE-no
mkdir $folder
mv H2O.ptcl.xml H2O.wfs.xml $folder
cp simple-H2O.xml $folder

cp qmcpack-theta.sub $folder
cd $folder
qsub qmcpack-theta.sub
cd ..

done
