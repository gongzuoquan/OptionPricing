#include<iostream>
#include<cstdlib>
#include"tree.hpp"
#include<cmath>
#include<cstdlib>
#include<chrono>
#include<openacc_curand.h>
#include<ctime>
using namespace std;
using namespace chrono;

int random_seed=123456;
string param_file="params.txt";

int main(int argc, char **argv)
{
    //srand(random_seed);
    srand(time(NULL));

    Params params(param_file);
    params.display();

    AmericanOption ap(params);
    ap.display();
    ap.MC_pricing(params.sim_num,params.time_step,params.rho,true);

    return 0;
}


