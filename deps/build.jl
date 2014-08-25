run(`ls`)
@osx_only run(`gfortran -shared -O3 ./doubleshift/src/*.f90 doubleshift/src/turnovers/*.f90 singleshift/src/*.f90 singleshift/src/turnovers/*.f90 -o libamvw.dylib`)
@windows_only run(`gfortran -shared -O3 doubleshift/src/*.f90 doubleshift/src/turnovers/*.f90 singleshift/src/*.f90 singleshift/src/turnovers/*.f90 -o libamvw.dll`)
@unit_only(`gfortran -shared -O3 doubleshift/src/*.f90 doubleshift/src/turnovers/*.f90 singleshift/src/*.f90 singleshift/src/turnovers/*.f90 -o libamvw.so`)
