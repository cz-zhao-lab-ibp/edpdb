# This makefile is for recompiling the program on the SuSE-linux computer.
# For a RedHat-linux computer, see the embedded comments, 
#  or use the following link
#ln -s /usr/lib/gcc/x86_64-redhat-linux/3.4.6/libg2c.a   /usr/lib/libg2c.a
#/usr/lib/gcc-lib/i386-redhat-linux/2.96/libg2c.a
PROGRAM		= edpdb_v13a

OPT		= 
ENV		= 

FC	   	= g77 
CC		= cc  

FFLAGS		=  ${OPT} $(ENV)  
CFLAGS		= -I/usr/local/include $(ENV)

LDFLAGS		=  
#LIBS		=  -lncurses  -lreadline /usr/lib/libg2c.a
LIBS		=  -lncurses  -lreadline /usr/lib/gcc/x86_64-redhat-linux/3.4.6/libg2c.a
		   
## for Redhat operating system, use the following
#                  /usr/lib/gcc-lib/i386-redhat-linux/2.96/libg2c.a
## to find out where the libg2c.a is located, type
#find /usr/ -name "libg2c*"	   

OBJS	      	= clique.o       rtnout.o \
		abcd.o        harker.o       pair.o  \
	        setw.o \
		atom.o        param.o  \
	        snayb.o \
		axis.o         eular.o  pickr.o  \
		chklib.o       find1.o     polar.o  \
				get_command.o   order.o          edpdb.o \
		listatom.o        needlemen1d.o  volume.o subr-linux-g77.o \
		thread.o smg.o

SRCS	      	= clique.f        rtnout.f \
		abcd.f        harker.f       pair.f  \
	        setw.f \
		atom.f     param.f  \
	        snayb.f \
		axis.f         eular.f  pickr.f  \
		chklib.f       find1.f     polar.f  \
		smg.c \
		get_command.f       order.f  edpdb.f \
  		listatom.f         needlemen1d.f  volume.f subr-linux-g77.f \
		thread.f smg.c

.SUFFIXES :	.o .c .f 
.f.o :
		$(FC) -c $(FFLAGS) $*.f 
.c.o :
		$(CC) -c $(CFLAGS) $*.c

$(PROGRAM):	$(OBJS)
		@echo -n "Loading $(PROGRAM) ... "
		@$(FC) $(FFLAGS)  $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)
		@echo "done"
