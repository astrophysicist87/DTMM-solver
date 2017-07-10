# ============================================================================
#  Makefile CFWR                             Chris Plumberg, September 8, 2015
# ============================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS += -lgsl -lgslcblas -g

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
LIBS            =       -I/home/kapustaj/cplumber/eigen-eigen-5a0156e40feb
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	main.e
endif

SRC		=	main.cpp \
			lib.cpp \
			gauss_quadrature.cpp \
			matrixInterpolation.cpp

INC		= 	lib.h \
			gauss_quadrature.h \
			matrixInterpolation.h

# -------------------------------------------------

OBJDIR		=	.
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	.

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(LIBS) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm $(OBJECTS)
		-rm $(MAIN)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
gauss_quadrature.cpp: gauss_quadrature.h
main.cpp: lib.h gauss_quadrature.h matrixInterpolation.h
lib.cpp: lib.h gauss_quadrature.h matrixInterpolation.h
matrixInterpolation.cpp: lib.h gauss_quadrature.h matrixInterpolation.h


