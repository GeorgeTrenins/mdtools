#!/bin/bash

ext=$(python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
echo $ext

for name in _rdfs_fort _angle_pdfs_fort ; do
  python -m numpy.f2py -c -m $name $name.f90
  CWD=$(pwd)
  cd ../mdtools/structure
  ln -sf ../../f90/${name}${ext}
  cd $CWD
done
