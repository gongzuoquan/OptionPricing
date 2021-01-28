CXX = pgc++

execute: main.o option.o util.o
	$(CXX) -o execute main.o option.o util.o -Minfo
main.o: main.cpp
	$(CXX) -c -acc -Minfo -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen 
option.o: option.cpp
	$(CXX) -c -acc -Minfo -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen
util.o: util.cpp
	$(CXX) -c -acc -Minfo -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen
clean:
	rm -rf *.o execute 

