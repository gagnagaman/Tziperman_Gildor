#####################################################################
###  Definitions of compiler and flags for hosts running ncarg:   ###
#####################################################################
BeachDbgFlags=-g -C -c -convert little_endian
BeachRunFlags=-O2 -C -c -convert little_endian
BeachFortran=f77
BeachLoader=f77
BeachRunLoaderFlags=-O2 -C 
BeachDbgLoaderFlags=-g -C 

#####################################################################
###  Set compiler and flags:  					  ###
#####################################################################

FORTRAN=$(BeachFortran)
LOADER=$(BeachLoader)
FLAGS=$(BeachRunFlags)  
LoaderFlags=$(BeachRunLoaderFlags)
LIBS=$(BeachLibs)
 
#####################################################################
###  Programs to be compiled by this makefile:			  ###
#####################################################################

all:: \
       coup51

#####################################################################
###  Compilation instructions for main programs:		  ###
#####################################################################

coup51: couple.o atmosphere.o Asm2dt1.o humidity.o atmosphere.h ocean.o Osm2dt1.o UNESCO.o ocean.h p4grads.o output.o couple.h seaice.o seaice.h landice.o landice.h biochemical.o biochemical.h biochem_coef.o biochem_coef.h states.o statet.o insolation1.o averaging.o zeromean.o
	${FORTRAN} ${LoaderFlags} -o coup51 \
        $(LIBS) \
                couple.o atmosphere.o Asm2dt1.o humidity.o ocean.o Osm2dt1.o UNESCO.o p4grads.o output.o seaice.o landice.o biochemical.o biochem_coef.o states.o statet.o insolation1.o averaging.o zeromean.o

couple.o: couple.f atmosphere.h ocean.h couple.h seaice.h landice.h biochemical.h averaging.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=1950 couple.f 

atmosphere.o: atmosphere.f atmosphere.h couple.h seaice.h ocean.h landice.h biochemical.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=1400 atmosphere.f

Asm2dt1.o: Asm2dt1.f  atmosphere.h couple.h biochemical.h
	${FORTRAN} ${FLAGS} Asm2dt1.f

humidity.o: humidity.f atmosphere.h couple.h
	${FORTRAN} ${FLAGS} humidity.f

ocean.o: ocean.f   ocean.h couple.h atmosphere.h landice.h seaice.h biochemical.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=1500 ocean.f

UNESCO.o: UNESCO.f  ocean.h couple.h
	${FORTRAN} ${FLAGS} UNESCO.f

states.o: states.f  ocean.h couple.h
	${FORTRAN} ${FLAGS} states.f

statet.o: statet.f  ocean.h couple.h
	${FORTRAN} ${FLAGS} statet.f

Osm2dt1.o: Osm2dt1.f  ocean.h couple.h biochemical.h
	${FORTRAN} ${FLAGS} Osm2dt1.f

p4grads.o: p4grads.f couple.h atmosphere.h ocean.h seaice.h
	${FORTRAN} ${FLAGS} p4grads.f

output.o: output.f couple.h atmosphere.h ocean.h seaice.h biochemical.h landice.h
	${FORTRAN} ${FLAGS} output.f

seaice.o: seaice.f couple.h atmosphere.h ocean.h seaice.h
	${FORTRAN} ${FLAGS} seaice.f

landice.o: landice.f couple.h atmosphere.h ocean.h landice.h
	${FORTRAN} ${FLAGS} landice.f

biochemical.o: biochemical.f atmosphere.h ocean.h couple.h biochemical.h biochem_coef.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=1700 biochemical.f

biochem_coef.o: biochem_coef.f biochem_coef.h ocean.h couple.h
	${FORTRAN} ${FLAGS} biochem_coef.f

insolation1.o: insolation1.f atmosphere.h couple.h seaice.h ocean.h landice.h biochemical.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=1400 insolation1.f

averaging.o: averaging.f averaging.h atmosphere.h ocean.h couple.h seaice.h landice.h biochemical.h
	${FORTRAN} ${FLAGS} -OPT:Olimit=2500 averaging.f

zeromean.o: zeromean.f averaging.h 
	${FORTRAN} ${FLAGS} -OPT:Olimit=1800 zeromean.f
