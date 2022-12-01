CC = g++
CFLAGS = -O3 --fast-math
LDFLAGS = 
a.out: main.o matrix.o
	$(CC) $(LDFLAGS) main.o matrix.o -o a.out
main.o: main.cpp matrix.h
	$(CC) $(CFLAGS) -c main.cpp
matrix.o: matrix.cpp matrix.h
	$(CC) $(CFLAGS) -c matrix.cpp
clean:
	rm -f *.o a.out
