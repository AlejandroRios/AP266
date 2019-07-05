# compile3.sh
# compiles fortran3.F90 and mymodule.F90

PROG=llt_main
MOD1=llt_module
MOD2=llt_module_d
MOD3=llt_module_b


# Remove previous files
rm *.o
rm *.exe

# Compilation
gfortran -c -fdefault-real-8 $MOD1.F90 $MOD2.F90 $MOD3.F90 $PROG.F90

# Compile most used files from the ADFirstAidKit
gcc -c $CFLAGS ADFirstAidKit/adStack.c
gfortran -c $FFLAGS ADFirstAidKit/adBuffer.f

# Linking
gfortran adStack.o adBuffer.o $MOD1.o $MOD2.o $MOD3.o $PROG.o -o $PROG.exe

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
