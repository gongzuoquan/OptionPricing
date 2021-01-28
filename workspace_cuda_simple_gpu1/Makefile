CXX = pgc++
NVCC = /usr/local/cuda/bin/nvcc

execute: main.o option.o util.o
	$(CXX) -o execute main.o option.o util.o -Minfo -Mcuda
main.o: main.cpp
	$(CXX) -c -acc -Minfo -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen 
option.o: option.cu
	$(NVCC) -c -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen
util.o: util.cpp
	$(CXX) -c -acc -Minfo -std=c++11 $< -I/home/zqgong/workspace/NLPwork/eigen
clean:
	rm -rf *.o execute 

