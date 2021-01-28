CXX = pgc++

execute: main.o
	$(CXX) -o execute main.o -Minfo -Mcuda -ta=tesla
main.o: main.cpp
	$(CXX) -c -acc -Minfo -Mcuda -std=c++11 main.cpp

clean:
	rm -rf *.o execute 

