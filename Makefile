OBJS = \
	subs.o 
FFILE = set3d
FLAGS = -O3 -fdefault-real-8
STLFILE = cube40.stl

compile: $(OBJS)
	gfortran-mp-4.8 ${FFILE}.f90 ${OBJS} $(FLAGS) -o $(FFILE).exec

file: 
	vim -O Makefile /Desktop/levelSet/Makefile

run:
	./$(FFILE).exec $(STLFILE)

.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	gfortran-mp-4.8 $< $(FLAGS) -c -o $@

clean: 
	rm -f *.mod
	rm -f *.exec 
	rm -f *.o
	rm -f *.txt
	rm -rf output
