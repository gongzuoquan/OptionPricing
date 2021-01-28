//2020.1.28
//版本：单个gpu+option模板
//规模：100000 * 50
//期权类型：美式看跌
//方法：最小二乘 Monte Carlo 方法 （LSMC）
//注：验证数据参考LSMC原始论文（price=5.622）
#include<iostream>
#include<ctime>
#include<cmath>
#include"util.hpp"
using namespace std;
//第一步：按照BS模型生成价格路径，并保存在二维矩阵
//第二步：从最后时刻开始回溯：
//1.每个时刻都需要比较两个payoff：当前行权的payoff和当前持有的payoff
//2.当前行权的payoff直接用当前价格即可算出
//3.持有的payoff由最小二乘计算出
//4.二者比较决定是否在当前行权
//第三步：得出每条期权路径的最佳执行时刻，及相应的最佳payoff
//第四步：计算期望，得到定价结果

//#define Time_Step 200 //每条期权的时间采样点个数
//#define N (int)pow(2,20) //期权样本数
#define Time_Step 50 //每条期权的时间采样点个数
#define N 100000 //期权样本数
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

// 正态分布随机数
dtype gaussrand(){
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    return sqrt(-2*log(v))*sin(2*M_PI*u);
}

int main()
{
    clock_t start=clock();
    srand(time(NULL));

    TOption option(44,40,0.06,0.4,2.0);
    option.display();
    cout<<N<<" * "<<Time_Step<<endl;

    int i, j;
    dtype sum=0.0;
    dtype price;
    dtype dt=option.T/Time_Step;

    dtype asset_mat[N][Time_Step];
    dtype payoff_mat[N][Time_Step];

    // ***第一个循环（二重循环）：MC 方法生成N条样本路径
    for(i=0; i<N; i++){
        dtype asset_price=option.S0;
        // 一条样本路径的生成
        for(j=1; j<Time_Step; j++){
            asset_price = asset_price * exp((option.R-0.5*option.V*option.V)*dt+option.V*gaussrand()*sqrt(dt));
            asset_mat[i][j]=asset_price;  //存储资产价格在二维数组中

            // 计算一条样本的payoff
            dtype payoff=0.0;
            if(option.K-asset_price>0)
                //payoff = asset_price - option.K;
                payoff = option.K-asset_price;
            else  payoff = 0.0;
            payoff_mat[i][j]=payoff;
        }
    }
    //asset_mat[i][j] 中存储着标的资产的价格
    //payoff_mat[i][j] 中存储着在每个时刻行权的payoff

    // ***第二个循环：创建数组记录执行时刻，默认执行时刻为最后时刻
    int exec_time[N];
    for(int e=0;e<N;e++)
    {
        exec_time[e]=Time_Step-1;
    }

    dtype ls_x[N];
    dtype ls_y[N];
    // ***第三个循环：最重要的循环，通过最小二乘计算持有payoff，并将其与执行payoff比较，以决策是否在当前时刻执行期权
    for(int t=Time_Step-2;t>=0;t--) //最后一个时刻是Time_Step-1，但最后时刻没有持有价值，所以不需要比较
    {
        int non_zero_num=0;
        //第一个循环是为了最小二乘计算持有价值做准备
        //需要与其pk的当前执行payoff已经在第一个循环里计算完毕了
        //这里，如果当前执行payoff为0，那么可以自动跳过，即没有收益就不用比了，自动认为继续持有胜出
        for(int i=0;i<N;i++)
        {
            if(payoff_mat[i][t]!=0)
            {
                //x——当前标的资产价格，y——当前最佳收益的折现
                ls_x[non_zero_num]=asset_mat[i][t];
                ls_y[non_zero_num]=payoff_mat[i][exec_time[i]]*exp(-option.R*dt*(exec_time[i]-t));
                non_zero_num+=1;
            }
        }

        // *获得到进行最小二乘的x和y数组后，开始进行最小二乘
        LeastSquareEigen ls(ls_x,ls_y,non_zero_num);

        //最小二乘完成，可以获得持续payoff，然后我们进行两个payoff的比较，以确定是否在当前时刻行权
        for(int i=0;i<N;i++)
        {
            if(payoff_mat[i][t]!=0)
            {
                dtype continuation_payoff=ls.get_Y(asset_mat[i][t]);
                if(payoff_mat[i][t]>continuation_payoff)
                {
                    exec_time[i]=t;
                }
            }
        }

    }

    // ***第四个循环：计算每个期权路径下的最佳payoff
    //每个时间的最佳payoff已经得出
    //当开始计算payoff均值作为price
    dtype payoff_array[N];
    for(int i=0;i<N;i++)
    {
        payoff_array[i]=payoff_mat[i][exec_time[i]]*exp(-option.R*dt*exec_time[i]);
        sum+=payoff_array[i];
    }

    // 计算N条路径的payoff的均值
    price = sum / N;

    // ***第五个循环：计算error（标准误差）
    dtype total=0.0;
    for(i=0; i<N; i++){
        dtype tmp = (payoff_array[i] - price)*(payoff_array[i] - price);
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
