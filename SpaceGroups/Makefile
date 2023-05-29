#######################################
#  FOPT    -  optimization
#  C_LIBS  -  C libraries
#  LDFLAGS -  linker flags
#  EXEC    -  name of out file
#  OBJS    -  name of obj files
#######################################

EXEC    = sgroup
OBJS    = $(EXEC).o math_sg.o sto.o io.o lat.o pgrp.o \
pgrp_op.o  rotb.o sgrp_cub.o sgrp_hex.o sgrp_tet.o  sgrp_ort.o sgrp_trg.o \
sgrp_mon.o sgrp_trc.o 

#FOPT    = -O3 -malign-double -malign-loops=4 -malign-jumps=4
FOPT    = 
LDFLAGS = -L/usr/lib 
C_LIBS  = -lm


# Make executable
$(EXEC) :  $(OBJS)
	cc -o $(EXEC) $(OBJS) $(LDFLAGS) $(C_LIBS)

clean :
	rm -f $(OBJS) $(EXEC)

.c.o:
	cc -c $(FOPT) $<






