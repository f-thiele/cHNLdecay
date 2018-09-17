CC = gcc
CXX = g++
RM = rm -f
CPPFLAGS = -g $(shell root-config --cflags) -I/usr/include
LDFLAGS = -g $(shell root-config --ldflags)
LDLIBS = $(shell root-config --libs) -lmpfr -lgmp -L/usr/lib -L/usr/lib64
LDLIBS_DBG = $(shell root-config --libs) -lmpfr -lgmp -L/usr/lib -L/usr/lib64 -lprofiler -ltcmalloc

SRCS = $(wildcard *.cxx)
OBJS = $(SRCS:.cxx=.o)

all: cHNLdecay

debug: $(OBJS)
	$(CXX) $(LDFLAGS) -o cHNLdecay $(OBJS) $(LDLIBS_DBG) -ggdb

cHNLdecay: $(OBJS)
	$(CXX) $(LDFLAGS) -o cHNLdecay $(OBJS) $(LDLIBS)

# To obtain object files
%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJS)
