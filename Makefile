gravity.so: gravity.f90
	f2py3 -m gravity -c gravity.f90 --f90flags=-fopenmp -lgomp

clean:
	rm -f *.so

