#include<iostream>
#include<vector>
using namespace std;
class MonteCarloOptionPricing
{
    public:
        double r;
        double S0;
        double K;
        double T;
        double mue;
        double sigma;
        double div_yield;
        int simulation_rounds;
        int time_steps;
        bool fix_random_seed;

        double *r_mat;
        double *sigma_mat;
        double dt;

        double prices;
    public:
        MonteCarloOptionPricing(
                double in_r,
                double in_S0,
                double in_K,
                double in_T,
                double in_mue,
                double in_sigma,
                double in_div_yield=0.0,
                int in_simulation_rounds=10000,
                int in_time_steps=4,
                bool in_fix_random_seed=false):r(in_r),S0(in_S0),K(in_K),T(in_T),mue(in_mue),sigma(in_sigma),div_yield(in_div_yield),simulation(in_simulation_rounds),time_steps(in_time_steps),fix_random_seed(in_fix_random_seed)
        /*
         :param S0: 当前标的资产价格（例如：股票）
         :param K: 执行价格（敲定价格）
         :param T: 期权成熟期（单位为年，可以为小数）
         :param r: 利率（这里我们假定为常数利率模型）
         :param sigma: 资产年收益率的波动率
         :param div_yield: 标的资产的年分红
         :param simulation_rounds: 通常，蒙特卡洛算法需要大量的迭代次数
         :param time_steps: 时间轴上的采样次数，在0到T之间的切片次数，例如252
         :param fix_random_seed: 一个bool值，True或者False
        */
        {
            assert(sigma>=0);
            assert(S0>=0);
            assert(T>=0);
            assert(div_yield>=0);
            assert(time_steps>=0);
            assert(simulation_rounds>=0);
            int sample_num=simulation_rounds*time_steps;
            r_mat=new double [sample_num];
            for(int i=0;i<sample_num;i++)
            {
                r_mat[i]=r/(T*time_steps);
            }

            sigma_mat=new double [sample_num];
            for(int i=0;i<sample_num;i++)
            {
                sigma_mat[i]=sigma;
            }
            dt=T/time_steps;

            if(fix_random_seed)
                srand(15000);
        }

        ~MonteCarloOptionPricing()
        {
            delete r_mat;
            delete sigma_mat;
        }
    public:
        void vasicek_model(double r0, double alpha, double b, double interest_vol)
        /*
         用于利率模拟的vasicek模型
         :param r0:
         :param alpha:
         :param b:
         :param interest_vol:
         :return:
         */
        {
            //TODO
            return;
        }
        void Cox_Ingersoll_Ross_model(double r0, double alpha, double b, double interest_vol)
        {
            return ;
        }
};




















