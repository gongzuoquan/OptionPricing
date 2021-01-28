CXX = pgc++

execute: main.o
	$(CXX) -o execute main.o -ta=tesla -Mcuda -Minfo
main.o: main.cpp
	$(CXX) -c -acc -Minfo -Mcuda=nollvm -std=c++11 main.cpp

clean:
	rm -rf *.o execute *.prof

