is_h5=`echo $2 | grep "\.h5" | wc -l`

if [ $is_h5 -gt 0 ]
then
  if [ -e $2 ]
  then
    echo remove $2
    rm $2
  fi
  if [ -e $1 ]
  then
    ln -s ../$1 $2
    echo link ../$1 to $2
  else
    echo "$1 doesn't exist"
    exit 1
  fi
else
  echo Wrong file $2. Do nothing!
  exit 1
fi


