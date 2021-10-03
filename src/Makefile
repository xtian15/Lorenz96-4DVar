SRC=lorenzN.f90
OBJ=lorenz.so
FLAG="-fdefault-real-8"

${OBJ}:		${SRC}
		f2py --f90flags=${FLAG} -m lorenz -c ${SRC}
clean:
		rm *.so
