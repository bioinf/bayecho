CXX = g++
CXXFLAGS = -Wall -O3 -w
HEADERS = *.hpp
LIB = collection_utils.o distances.o math_func.o hammer.o DNASeq.o NeighborSet.o KmerHashMap.o
TARGETS = hashing NeighborJoin NeighborJoinParam Voting

%.o:	%.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $*.o

%:	%.cpp $(LIB) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $* $*.cpp $(LIB)

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
	cp DNASeq.cpp KmerHashMap.cpp hashing.cpp HashMerge.cpp NeighborSet.cpp NeighborJoin.cpp NeighborJoinParam.cpp NeighborMerge.cpp Voting.cpp release/
	cp Makefile release/
	cp sample_data.txt release/
	cp sample_answer.txt release/
	cp LICENSE release/
	cp README release/
