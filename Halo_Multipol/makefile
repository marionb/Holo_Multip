#make file for holo_multip
#author: Marion Baumgartner

#flags
OBJS = Param.o Data.o DataProcess.o main.o
CC=g++ #compiler to use
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all: holo_multip.exe

#create executable with the first target
holo_multip.exe: $(OBJS)
	$(CC) $(OBJS) -o run
DataProcess.o:DataProcess.cpp main.o
	$(CC) $(CFLAGS) DataProcess.cpp
Data.o:Data.cpp main.o
	$(CC) $(CFLAGS) Data.cpp
Param.o: Param.cpp main.o
	$(CC) $(CFLAGS) Param.cpp

main.o: main.cpp Param.h Data.h
	$(CC) $(CFLAGS) main.cpp

clean:
	\rm *.o run

tar:
	tar cfv task.tar hellomake.c hellofunc.c hellomake.h

