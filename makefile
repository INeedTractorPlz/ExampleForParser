CC = g++
CCFLAGS = -c --std=c++1z -frtti -lm -O2

all : Eskin.o
	$(CC) Eskin.o --std=c++1z -frtti -O2 -o main

Eskin.o: Eskin.cpp
	$(CC) $(CCFLAGS) Eskin.cpp

clean:
	rm -rf *.o main