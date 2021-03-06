all: mst-backbone UNRESTselector_perEdge

mst-backbone: mst-backbone.o
	g++ -o mst-backbone mst-backbone.o -lstdc++fs

UNRESTselector_perEdge: UNRESTselector_perEdge.o
	g++ -o UNRESTselector_perEdge UNRESTselector_perEdge.o -lstdc++fs

mst-backbone.o: mst-backbone.cpp
	g++ -c mst-backbone.cpp -o mst-backbone.o -O3 -Ieigen3 -Iboost -I. -std=c++17

UNRESTselector_perEdge.o: UNRESTselector_perEdge.cpp
	g++ -c UNRESTselector_perEdge.cpp -o UNRESTselector_perEdge.o -O3 -Ieigen3 -Iboost -I. -std=c++17

clean:
	rm mst-backbone mst-backbone.o UNRESTselector_perEdge UNRESTselector_perEdge.o
