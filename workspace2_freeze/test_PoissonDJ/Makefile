CXX = pgc++

execute: main.o
	$(CXX) -o execute main.o -Minfo
main.o: main.cpp
	$(CXX) -c -acc -Minfo -std=c++11 main.cpp

clean:
	rm -rf *.o execute 

