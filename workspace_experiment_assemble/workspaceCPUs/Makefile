CXX = pgc++

execute: main.o
	$(CXX) -o execute main.o -ta=tesla -Mcuda -Minfo -L/home/zqgong/include/trng-4.17/src/.libs /home/zqgong/include/trng-4.17/src/.libs/libtrng4.a
main.o: main.cpp
	$(CXX) -c -acc -mp -Minfo -Mcuda=nollvm -std=c++11 main.cpp -I/home/zqgong/include/trng-4.17

clean:
	rm -rf *.o execute 

