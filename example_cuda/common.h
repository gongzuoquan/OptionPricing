#ifndef _COMMON_H_
#define _COMMON_H_
#include<stdio.h>
#define dtype float
typedef struct TOption
{
    dtype S0;
    dtype K;
    dtype T;
    dtype R;
    dtype V;
    public:
        TOption(dtype in_S0,dtype in_K,dtype in_T,dtype in_R,dtype in_V):S0(in_S0),K(in_K),T(in_T),R(in_R),V(in_V)
        {
            return;
        }
        ~TOption(){}
    public:
        void display()
        {
            printf("S0: %.2f\n",S0);
            printf("K: %.2f\n",K);
            printf("T: %.2f\n",T);
            printf("R: %.2f\n",R);
            printf("V: %.2f\n",V);
            printf("\n");
        }
}TOption;

#endif
