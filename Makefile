gravity.so: gravity.f90
	f2py3 -m gravity -c gravity.f90

test: gravity.so
	./test.py

clean:
	rm -f *.so

