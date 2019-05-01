gravity: gravity.so

clean:
	rm -f gravity.so
	rm -f *.pyc

gravity.so: gravity.f90
	f2py -m gravity -c gravity.f90 --f90flags=-fopenmp -lgomp
