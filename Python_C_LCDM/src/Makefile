EXEC = pybbn.so

OBJS = reaction_rates.o rk4.o timestep.o collect_derivs.o abundance_derivs.o dynamics_derivs.o densities.o bessel.o partial_dynamics_derivs.o initialize_constants.o expansion_rate.o

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