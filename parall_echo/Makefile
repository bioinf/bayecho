LDFLAGS+= -lboost_program_options -lpthread -ltbb
CXX = g++
CXXFLAGS = -Wall -O3 -w -fopenmp
HEADERS = *.h*
LIB = DNASeq.o NeighborSet.o KmerHashMap.o NeighborSetOldLoader.o
TARGETS = HashMerge NeighborMerge ParallelHashing ParallelNeighborJoin ParallelNeighborJoinParam ParallelVoting 

%.o:	%.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $*.o $(LDFLAGS) 

%:	%.cpp $(LIB) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $* $*.cpp $(LIB) $(LDFLAGS)

all:	$(LIB) $(TARGETS)

clean:
	rm -f $(LIB) $(TARGETS) $(UTIL)

release:
	mkdir release
	mkdir release/output
	mkdir release/tmp
	mkdir release/log
	cp ErrorCorrection.py release/
	cp util.hpp DNASeq.hpp KmerHashMap.hpp MMAP.hpp MMAPReads.hpp NeighborSet.hpp release/
	cp DNASeq.cpp KmerHashMap.cpp hashing.cpp HashMerge.cpp NeighborSet.cpp NeighborJoin.cpp NeighborJoinParam.cpp ParallelNeighborJoin.cpp ParallelNeighborJoinParam.cpp NeighborMerge.cpp Voting.cpp release/
	cp Makefile release/
	cp sample_data.txt release/
	cp sample_answer.txt release/
	cp LICENSE release/
	cp README release/
