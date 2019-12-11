f2py=python -m numpy.f2py

conda=conda
CONDA_PREFIX=$(shell $(conda) info --base)
activate=$(CONDA_PREFIX)\Scripts\activate

.PHONY: all clean new

all:
	@echo "Use 'make linux' or 'make windows' to build this library"

linux: dependencies interface minpack.f
	$(f2py) minpack.pyf minpack.f -L. -ldependencies -c

windows: dependencies interface minpack.f
	$(conda) create -y -n f2py python=3.7
	cp distutils.cfg $(CONDA_PREFIX)\\envs\\f2py\\Lib\\distutils
	$(activate) f2py && \
	$(conda) install -y numpy libpython m2w64-toolchain && \
	$(f2py) --compiler=mingw32 minpack.pyf minpack.f -L. -ldependencies -c

interface: minpack.f
	$(f2py) -m minpack -h minpack.pyf minpack.f

backup-interface: minpack.pyf
	cp minpack.pyf minpack.pyf.bak

clean:
	rm -fr *.o dependencies/*.o libdependencies.a

clean-interface:
	rm -fr minpack.pyf

clean-libs:
	rm -fr minpack*.so

new: clean all

new-interface: backup-interface clean-interface interface

# Dependencies

dependencies: libdependencies.a

files = $(patsubst %.f,%.o,$(wildcard dependencies/*.f))
libdependencies.a: $(files)
	ar crs libdependencies.a $^

%.o: %.f
	gfortran -c -fPIC $< -o $@
