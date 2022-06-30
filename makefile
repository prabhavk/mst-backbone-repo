all: mst-backbone UNRESTselector

mst-backbone: mst-backbone.o
	g++ -o mst-backbone mst-backbone.o -lstdc++fs

UNRESTselector: UNRESTselector.o
	g++ -o UNRESTselector UNRESTselector.o -lstdc++fs

mst-backbone.o: mst-backbone.cpp
	g++ -c mst-backbone.cpp -o mst-backbone.o -O3 -Ieigen3 -Iboost -I. -std=c++17

UNRESTselector.o: UNRESTselector.cpp
	g++ -c UNRESTselector.cpp -o UNRESTselector.o -O3 -Ieigen3 -Iboost -I. -std=c++17

clean:
	rm mst-backbone mst-backbone.o UNRESTselector UNRESTselector.o
