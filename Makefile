VERSION=1.0.2
CC = gcc
CFLAGS = -O3 -funroll-loops -Wall -fopenmp -finline-limit-50000
LFLAGS = -lm -pthread
TMPFILE_DIR=.

DEF = -DTMPFILE_DIR=\"$(TMPFILE_DIR)\" -DVERSION=\"$(VERSION)\"

DEBUG_FLAGS = -g -DDEBUG

PRG = platanus
OBJ = main.o common.o assemble.o scaffold.o gap_close.o bubble_map.o correct.o mcounter.o mgraph.o lcounter.o lgraph.o seqlib.o mapper.o
SRC = main.c common.c assemble.c scaffold.c gap_close.c bubble_map.c correct.c mcounter.c mgraph.c lcounter.c lgraph.c seqlib.c mapper.c

all: $(PRG)

$(PRG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(CFLAGS) $(LFLAGS)

.c.o:
	$(CC) -o $@ -c $< $(CFLAGS) $(DEF)

debug: 
	make all CFLAGS="${CFLAGS} ${DEBUG_FLAGS}"

prof: 
	make all CFLAGS="${CFLAGS} -pg"

clean:
	rm -f $(PRG) $(OBJ)

depend:
	makedepend $(SRC)
