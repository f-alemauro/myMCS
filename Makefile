CC=g++
CFLAGS=-c -I/usr/share/R/include
SOURCES=fmcs_R_wrap.cpp util.cpp MCSCompound.cpp MCSMap.cpp MCS.cpp MCSRingDetector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mcswrap

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)
