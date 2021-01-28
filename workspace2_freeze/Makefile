CXX = pgc++

execute: main.o
	$(CXX) -o execute main.o -Minfo
main.o: main.cpp
	$(CXX) -c -acc -Minfo -std=c++11 main.cpp -I ~/workspace/NLPwork/eigen

clean:
	rm -rf *.o execute 

