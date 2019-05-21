CC=gcc
CFLAGS=-g
CLIB=-lm

OBJ=pureins.bin
$(OBJ): main.o ins.o insio.o
	$(CC) main.o ins.o insio.o $(CLIB) -o $(OBJ)

UT_OBJ=unittest.bin
ut: unittest.o ins.o
	$(CC) unittest.o ins.o -g -lm -lcriterion -o $(UT_OBJ)
	./$(UT_OBJ)

unittest.o: unittest.c
main.o: main.c
ins.o: ins.c ins.h
insio.o: insio.c ins.h

clean:
	trash *.o $(OBJ) $(UT_OBJ)
test:
	./$(OBJ) > ./data/tmp_output.csv
