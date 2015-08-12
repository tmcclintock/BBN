CC=gcc
CFLAGS=-c -Wall -O2 -I.
DEPS=physicalConstants.h programConstants.h
LDFLAGS=-lm
VPATH=AbundanceDerivs:Bessel:Constants:CriticalDerivatives:Densities:DriverRoutines:InitialConditions:MainFile:OutputRoutines:PartialDerivatives:ReactionRates:TimeSteps
SOURCES=mainBBN.c geAbundanceDerivs.c bessel.c criticalDerivs.c densDerivs.c densities.c driver.c initializeAll.c initialAbundances.c initialDynamics.c outputAbundances.c outputDynamics.c outputHeader.c outputRates.c partials.c reactionRates.c rk4BBN.c timeStep.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=NewBBN

all:$(SOURCES) $(EXECUTABLE) $(DEPS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $^ -o $@

move:
	mv *.o ObjectFiles/

clean:
	rm -f *.o *~ ObjectFiles/* $(EXECUTABLE)