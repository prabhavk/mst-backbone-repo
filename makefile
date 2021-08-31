mst-backbone: mst-backbone.o
        g++ -o mst-backbone mst-backbone.o -lstdc++fs

mst-backbone.o: main.cpp
        g++ -c main.cpp -o mst-backbone.o -O3 -Ieigen3 -Iboost -I. -std=c++17

clean:
        rm mst-backbone mst-backbone.o
