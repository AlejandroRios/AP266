# compile3.sh
# compiles fortran3.F90 and mymodule.F90

PROG=llt_main
MOD1=llt_module


# Remove previous files
rm *.o
rm *.exe

# Compilation
gfortran -c -fdefault-real-8 $MOD1.F90 $PROG.F90

# Linking
gfortran $MOD1.o $PROG.o -o $PROG.exe

# Execution
./$PROG.exe

#PROG=fortran3

# Remove previous files
#rm *.o
#rm *.exe

# Compilation
#gfortran -c -fdefault-real-8 $PROG.F90

# Linking
#gfortran $PROG.o -o $PROG.exe

# Execution
#./$PROG.exe
