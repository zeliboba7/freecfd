CC=g++
#OPTS=-Wall -g
OPTS=
EXE=read

INCLUDES= -I$(CGNS_BASE)/include 
LIBS= -L$(CGNS_BASE)/lib 

SRC= \
	grid_read.cc
.SUFFIXES: .cc .o

OBJ=	$(SRC:.cc=.o)

.cc.o:
	$(CC) $(OPTS) $(INCLUDES) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIBS) -lcgns

clean:	
	rm -rf $(OBJ) $(EXE)
