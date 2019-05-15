TARGET = cindo.exe
FC = gfortran
RM = del

FFLAGS = -g -O
LFLAGS = blas.dll lapack.dll
SRC = factorials.f90 legendre.f90 basegen.f90 \
	binomial.f90 pol2mul.f90 getycoef.f90 getzcoef.f90 \
	aint.f90 bint.f90 rovint.f90 redovp.f90 trnsform.f90 \
	coul_int.f90 core_int.f90 dipint.f90 scf_rhf.f90 \
	scf_uhf.f90 properties.f90 dipind.f90 dipmom_rhf.f90 \
	dipmom_uhf.f90 plot_1d_rhf.f90 plot_2d_rhf.f90 plot_1d_uhf.f90 \
	plot_2d_uhf.f90 spectrum.f90 basfunc.f90 cindo.f90
OBJ = factorials.o legendre.o basegen.o \
	binomial.o pol2mul.o getycoef.o getzcoef.o \
	aint.o bint.o rovint.o redovp.o trnsform.o \
	coul_int.o core_int.o dipint.o scf_rhf.o \
	scf_uhf.o properties.o dipind.o dipmom_rhf.o \
	dipmom_uhf.o plot_1d_rhf.o plot_2d_rhf.o plot_1d_uhf.o \
	plot_2d_uhf.o spectrum.o basfunc.o cindo.o

.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FFLAGS) $<
	
.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) -o $(TARGET) $(LFLAGS) $^

.PHONY: clean
clean:
		-$(RM) $(TARGET)
		-$(RM) $(OBJ)

