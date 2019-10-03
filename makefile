CC 	= gcc
CPP = g++
AR  = ar

CFLAGS = -std=c99 -O2
CPPFLAGS = -std=c++98 -O2
CLINK = -static
CPPLINK = -static

yins_bin: main.o libyins.a
	$(CPP) $(CPPLINK) main.o -o yins_bin -L. -lyins

main.o: main.cpp yins_core/ins.h
	$(CPP) $(CPPFLAGS) -c main.cpp -o main.o

static: libyins.a

libyins.a: ins.o inscmn.o insio.o inskf.o
	$(AR) -r libyins.a ins.o inscmn.o insio.o inskf.o

ins.o: yins_core/ins.c yins_core/ins.h
	$(CC) $(CFLAGS) -c yins_core/ins.c -o ins.o

inscmn.o: yins_core/inscmn.c yins_core/ins.h
	$(CC) $(CFLAGS) -c yins_core/inscmn.c -o inscmn.o

insio.o: yins_core/insio.c yins_core/ins.h
	$(CC) $(CFLAGS) -c yins_core/insio.c -o insio.o

inskf.o: yins_core/inskf.c yins_core/ins.h
	$(CC) $(CFLAGS) -c yins_core/inskf.c -o inskf.o

dist-clean: clean
	rm libyins.a yins_bin

clean:
	rm ins.o inscmn.o insio.o inskf.o main.o
