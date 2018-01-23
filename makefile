#CFLAGS = -g -Wall
CFLAGS = -Wall
CC = gcc
VGL = prof
OPTIONS = $(OPT)

ifeq ($(strip $(OPT)),$(VGL))

	OPTIONS = -fopenmp
	CC = kinst-ompp gcc

endif

ifeq ($(strip $(OPT)),)

	OPTIONS = -fopenmp -ffast-math -O3
	CC = gcc

endif




	


LDFLAGS = -lm -I../FEM/LinA/include/
LIB = -L../FEM/LinA/lib -lLinAlg
#CFILES = test.c FEM.c FEM2D.c File_stuff.c geometry2D.c linear_algebra.c gmres.c utils.c PCGS.c Multigrid.c
#CFILES = slimer.c FEM.c FEM2D.c File_stuff.c geometry2D.c linear_algebra.c gmres.c utils.c PCGS.c Multigrid.c stack.c misc.c red_black_tree.c GMRES_Newton.c
CFILES = buenorovio.c fft.c
OFILES = $(CFILES:%.c=%.o)
ARCH = $(shell arch)
OLD_ARCH = $(shell cat makefile.log) 

#LIBS_A = ~/C-Programme/FEM/AMD/Lib/libamd.a ~/C-Programme/FEM/UMFPACK/Lib/libumfpack.a
#BLA = -I/home/radszuweit/C-Programme/FEM/AMD/Include -I/home/radszuweit/C-Programme/FEM/UMFPACK/Include \
#-I/home/radszuweit/C-Programme/FEM/UFconfig

#gcc $(CFLAGS) -o test $(CFILES) $(LDFLAGS)
#@echo $(A) wenn echo nicht ausgeben

fast: 
	make look

debug: 

	#CFLAGS = -Wall -g 
	make look

profile:

	#CFLAGS = -Wall -pg -g
	make look

look:
ifneq ($(strip $(ARCH)), $(strip $(OLD_ARCH)))
	@echo new architecture: $(OLD_ARCH) to $(ARCH)
	make clean
else 
	@echo $(ARCH) architecture
endif	
	make compile

compile: LinAlg.a $(OFILES)

	$(CC) $(CFLAGS) $(OPTIONS) -o buenorovio$(ARCH) $(OFILES) $(LDFLAGS) $(LIB) 

%.o : %.c
	$(CC) $(CFLAGS) $(OPTIONS) -c $< $(LIBS_A) $(BLA) $(LDFLAGS)

	@echo $(ARCH) > makefile.log

LinAlg.a:
	cd ../FEM/LinA && make OPT="-g $(OPT)"

clean:
	rm -f *.o buenorovio

