FC = mpif90 -ffree-line-length-none

all: clean
	$(FC) -c initialization.f90
	$(FC) -c staticfields.f90 loop.f90
	$(FC) main.f90 initialization.f90 staticfields.f90 loop.f90 -o surya

clean:
	rm -f *.o *.mod surya
