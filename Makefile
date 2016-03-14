CC=g++
CFLAGS=-c -I/usr/share/R/include
SOURCES=src/fmcs_R_wrap.cpp src/util.cpp src/MCSCompound.cpp src/MCSMap.cpp src/MCS.cpp src/MCSRingDetector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mcswrap

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@

clean:
	$(RM) src/*.o src/*~ $(MAIN)
