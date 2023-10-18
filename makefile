mst-backbone: mst-backbone.o
	g++ -o mst-backbone mst-backbone.o -lstdc++fs

UNRESTselector_perEdge: UNRESTselector_perEdge.o
	g++ -o UNRESTselector_perEdge UNRESTselector_perEdge.o -lstdc++fs

all: mst-backbone UNRESTselector_perEdge

mst-backbone.o: mst-backbone.cpp
	g++ -c mst-backbone.cpp -o mst-backbone.o -O3 -Ieigen3 -Iboost -I. -std=c++11

UNRESTselector_perEdge.o: UNRESTselector_perEdge.cpp
	g++ -c UNRESTselector_perEdge.cpp -o UNRESTselector_perEdge.o -O3 -Ieigen3 -Iboost -I. -std=c++11

clean:
	rm mst-backbone mst-backbone.o UNRESTselector_perEdge UNRESTselector_perEdge.o

test: mstbackbone
	./mstbackbone --seq test.fas --constraint_size 10 --distance_measure LogDet --out test --root_supertree no