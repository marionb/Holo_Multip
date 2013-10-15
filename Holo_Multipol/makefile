#make file for holo_multip
#author: Marion Baumgartner

#flags
OBJS = Param.o Data.o Multipole.o DOmega.o  main.o
CC=g++ #compiler to use
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all: holo_multip.exe

#create executable with the first target
holo_multip.exe: $(OBJS)
	$(CC) $(OBJS) -o run

Multipole.o: Multipole.cpp specialVector.h main.o
	$(CC) $(CFLAGS) Multipole.cpp

Data.o: Data.cpp main.o
	$(CC) $(CFLAGS) Data.cpp

Param.o: Param.cpp main.o
	$(CC) $(CFLAGS) Param.cpp

DOmega.o: DOmega.cpp main.o
	$(CC) $(CFLAGS) DOmega.cpp

main.o: main.cpp Param.h Data.h Multipole.h DOmega.h 
	$(CC) $(CFLAGS) main.cpp

clean:
	\rm *.o run

#TODO tar:
#	tar cfv task.tar 

