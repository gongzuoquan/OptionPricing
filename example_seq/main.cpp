//2020.12.27
//版本：单个cpu+option模板
//规模：2^20 * 200
//期权类型：欧式看跌
#include<iostream>
#include<ctime>
#include<cmath>
using namespace std;

#define Time_Step 200 //每条期权的时间采样点个数
#define N (int)pow(2,20) //期权样本数
#define dtype double
typedef struct TOption
{
    dtype S0 = 100;
    dtype K = 100;
    dtype R = 0.05;
    dtype V = 0.2;
    dtype T = 1.0;

    public:
        TOption(){
            S0=0.0;
            K=0.0;
            R=0.0;
            V=0.0;
            T=0.0;
        }
        TOption(dtype in_S0,dtype in_K,dtype in_R,dtype in_V,dtype in_T):S0(in_S0),K(in_K),R(in_R),V(in_V),T(in_T)
        {}
    public:
        void display()
        {
            cout<<"S0: "<<S0<<endl;
            cout<<"K : "<<K<<endl;
            cout<<"R : "<<R<<endl;
            cout<<"V : "<<V<<endl;
            cout<<"T : "<<T<<endl;
        }
}TOption; //期权协议模板（template of option agreement）-结构体

dtype gaussrand(){
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    return sqrt(-2*log(v))*sin(2*M_PI*u);
}

int main()
{
    clock_t start=clock();
    srand(time(NULL));

    TOption option(100.0,100.0,0.05,0.2,1.0);
    option.display();
    cout<<N<<" * "<<Time_Step<<endl;

    int i, j;
    dtype sum=0.0;
    dtype aver;
    dtype discount_factor;
    dtype dt=option.T/Time_Step;
    dtype P[N];

    // MC 方法生成N条样本路径
    for(i=0; i<N; i++){
        dtype asset_price=option.S0;
        // 一条样本路径的生成
        for(j=1; j<Time_Step; j++){
          asset_price = asset_price * exp((option.R-0.5*option.V*option.V)*dt+option.V*gaussrand()*sqrt(dt));
        }

        // 计算一条样本的payoff
        dtype payoff=0.0;
        //if(asset_price-option.K>0)
        if(option.K-asset_price>0)
             //payoff = asset_price - option.K;
             payoff = option.K-asset_price;
        else  payoff = 0.0;

        P[i] = payoff;
        sum += payoff;
    }

    // 计算N条路径的payoff的均值
    aver = sum / N;
    //cout<<"aver: "<<aver<<endl;
    // 对payoff均值进行折现，获得当前时间的期权价格
    discount_factor=exp(-option.R*option.T);
    dtype price = aver * discount_factor;

    dtype total=0.0;
    for(i=0; i<N; i++){
        dtype tmp = (P[i] - aver)*(P[i] - aver);
        total+=tmp;
    }
    dtype variance=total/N;
    dtype error=sqrt(variance/N)*exp(-2*option.R*option.T);

    printf("the price = %lf\n",price);
    printf("the error = %lf\n",error);
    cout<<"95%置信区间为："<<"["<<price-error<<","<<price+error<<"]"<<endl;

    clock_t end=clock();
    dtype time=((end-start)/(dtype)CLOCKS_PER_SEC);
    printf("time: %f\n",time);
    return 0;
}
