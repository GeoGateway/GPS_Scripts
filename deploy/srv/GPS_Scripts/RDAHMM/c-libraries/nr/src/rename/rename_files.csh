#!/tools/bin/tcsh

foreach func (`cat func_list`)
  foreach file (nr_*.c)
    sed -e s/$func/NR_$func/g $file > out/$file
    mv out/$file $file
  end
end

foreach func (`cat func_list2`)
  foreach file (nr_*.c nr_*.h)
    sed -e s/$func/NR_$func/g $file > out/$file
    mv out/$file $file
  end
end

foreach file (nr_*.c nr_*.h)
  sed -e s/nrerror/NR_error/g $file > out/$file
  mv out/$file $file
  sed -e s/fcomplex/nr_fcomplex/g $file > out/$file
  mv out/$file $file
  sed -e s/nrutil\.h/nr_util\.h/g $file > out/$file
  mv out/$file $file
  sed -e s/complex\.h/nr_complex\.h/g $file out/$file
  mv out/$file $file
end
