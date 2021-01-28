CXX = pgc++
NVCC = /usr/local/cuda/bin/nvcc

execute: main.o option.o
	$(CXX) -o execute main.o option.o -Minfo -Mcuda
main.o: main.cpp
	$(CXX) -c -Minfo main.cpp
option.o: option.cu
	$(NVCC) -c option.cu

clean:
	rm -rf *.o execute 

