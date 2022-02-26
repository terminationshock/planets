gravity.so: gravity.f90
	f2py3 -m gravity -c gravity.f90

clean:
	rm -f *.so

