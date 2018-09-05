CC = gcc
CXX = g++
RM = rm -f
CPPFLAGS = -g $(shell root-config --cflags) -I/usr/include
LDFLAGS = -g $(shell root-config --ldflags)
LDLIBS = $(shell root-config --libs) -lmpfr -lgmp -L/usr/lib

SRCS = $(wildcard *.cxx)
OBJS = $(SRCS:.cxx=.o)

all: cHNLdecay

cHNLdecay: $(OBJS)
	$(CXX) $(LDFLAGS) -o cHNLdecay $(OBJS) $(LDLIBS) -ggdb

# To obtain object files
%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJS)
