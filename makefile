CC=gcc
LIBS =  -lm

practica1: practica1.c
	$(CC) $(LIBS) $^ -o $@.o

lines: random_line.c
	$(CC) $(LIBS)  $^ -o $@.o

practica2: practica2.c
	$(CC) $(LIBS)  $^ -o $@.o

hermite: hermite_curve.c
	$(CC) $(LIBS)  $^  -o $@.o

bezier: bezierCurve.c
	$(CC) $(LIBS)  $^ -o $@.o

bezierSurface: bezierSurface.c
	$(CC) $(LIBS) -framework opencl $^  -o $@.o 

.PHONY: clean

clean:
	rm -f ./*.o
