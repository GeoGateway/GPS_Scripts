#!/tools/bin/tcsh

foreach file (nr_*.c)
  sed -e s/=vector/=NR_vector/g $file > out/$file
  mv out/$file $file
  sed -e s/=dvector/=NR_dvector/g $file > out/$file
  mv out/$file $file
  sed -e s/=lvector/=NR_lvector/g $file > out/$file
  mv out/$file $file
  sed -e s/=cvector/=NR_cvector/g $file > out/$file
  mv out/$file $file
  sed -e s/=ivector/=NR_ivector/g $file > out/$file
  mv out/$file $file
  sed -e s/=matrix/=NR_matrix/g $file > out/$file
  mv out/$file $file
  sed -e s/=dmatrix/=NR_dmatrix/g $file > out/$file
  mv out/$file $file
  sed -e s/=cmatrix/=NR_cmatrix/g $file > out/$file
  mv out/$file $file
  sed -e s/=imatrix/=NR_imatrix/g $file > out/$file
  mv out/$file $file
end
