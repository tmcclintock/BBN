EXEC = reaction_rates.so

OBJS = reaction_rates.o

CC = gcc

CFLAGS = -fPIC

OFLAGS = -shared -Wl,-soname=$(EXEC)

INCL = 

LIBS = 

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~