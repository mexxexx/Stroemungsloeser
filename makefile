CC = c99
CFLAGS = -Wall

%.o: %.c
	$(CC) $(CFLAGS) -c $< 

OBJ = init.o boundary.o uvp.o visual.o main.o

simulation: $(OBJ)
	$(CC) $(CFLAGS) -o simulation.exe $(OBJ) -lm
	
init.o			:	init.h
boundary.o		:	init.h boundary.h
uvp.o 			:	init.h boundary.h uvp.h 
visual.o		:	init.h visual.h
main.o			:	init.h boundary.h uvp.h visual.h

clean:
	rm -rf *.o *.dat simulation *.prt *.prt_s *.prt_t *.stackdump gmon.out
	
rmdat:
	rm -rf *.dat *.part *.prt_s *.prt_t