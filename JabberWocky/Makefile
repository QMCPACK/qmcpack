

CC=mpicc
CFLAGS=-c -Wall
LDFLAGS=
#SOURCES=JabberWocky.c bp_open.c utils.c
SOURCES=JabberWocky.c process_file.c bp_open.c utils.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=JabberWocky

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o JabberWocky
