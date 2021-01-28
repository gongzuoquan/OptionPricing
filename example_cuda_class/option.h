#ifndef OPTION_CLASS_H_
#define OPTION_CLASS_H_

#define dtype float
class Option
{
    public:
        dtype S0;
        dtype K;
        dtype T;
        dtype R;
        dtype V;

        dtype price;
    public:
        Option(dtype in_S0,dtype in_K,dtype in_T,dtype in_R,dtype in_V):S0(in_S0),K(in_K),T(in_T),R(in_R),V(in_V){}
        ~Option(){}
    public:
        void display();
        void MC_pricing_gpu(const int n,const int m,unsigned int seed);
};

#endif
