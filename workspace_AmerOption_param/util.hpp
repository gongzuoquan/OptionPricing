#ifndef UTIL_H_
#define UTIL_H_
#include<cstdlib>
#include<cmath>
#include"Eigen/Dense"
#include"Eigen/Cholesky"
#include<iostream>
#include<string>
#include<vector>
using namespace std;
using namespace Eigen;

#define dtype double
#define SND_NUM 0.3989422804
#define MT_N 624
//工具类，用于记录一些实用函数

//MT19937 待验证
class MT19937{ 
    bool isInit;
    unsigned int MT[MT_N];
    int index;
    void srand(int seed)
    {
        index=0;
        isInit=1;
        MT[0]=seed;
        for (int i=1;i<MT_N;i++)
        {
            unsigned int t=1812433253*(MT[i-1])^(MT[i-1]>>30)+i;
            MT[i]=t&0xffffffff;
        }
    }
    void generate()
    {
        for(int i=0;i<MT_N;i++)
        {
            unsigned int y=(MT[i]&0x80000000)+(MT[(i+1)%MT_N]&0x7fffffff);
            MT[i]=MT[(i+397)%MT_N]^(y>>1);
            if(y&1)
                MT[i]^=2567483615;
        }
    }
    unsigned int Rand()
    {
        if (!isInit)
            srand((int)time(NULL));
        if (index==0)
            generate();
        unsigned int y=MT[index];
        y=y^(y>>11);
        y=y^((y<<7)&2636928640);
        y=y^((y<<15)&4022730752);
        y=y^(y>>18);
        index=(index+1)%MT_N;
        return y;
    }
};

dtype gaussianBoxMuller()
{
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    return sqrt(-2*log(v))*sin(2*M_PI*u);
}

dtype gaussianPolarMethod()
{
    dtype x = 0.0;
    dtype y = 0.0;
    dtype euclidSq = 0.0;
    do{
        x = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        y = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        euclidSq = x*x+y*y;
    }while(euclidSq >= 1.0);

    return x*sqrt(-2*log(euclidSq)/euclidSq);
}

//用于同时返回两个返回值（gaussianBoxMuller2，gaussianPolarMethod2）
//其实使用引用参数也可以

pair<dtype,dtype> gaussianBoxMuller2()
{
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    dtype u=sqrt(-2*log(v))*cos(2*M_PI*u);
    dtype v=sqrt(-2*log(v))*sin(2*M_PI*u);
    return {u,v};
}

pair<dtype,dtype> gaussianPolarMethod2()
{
    dtype x = 0.0;
    dtype y = 0.0;
    dtype euclidSq = 0.0;
    do{
        x = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        y = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        euclidSq = x*x+y*y;
    }while(euclidSq >= 1.0);

    dtype u= x*sqrt(-2*log(euclidSq)/euclidSq);
    dtype v= y*sqrt(-2*log(euclidSq)/euclidSq);
    return {u,v};
}

//标准正态分布
dtype std_normal(dtype input)
{
    dtype output=(SND_NUM)*exp(-input*input/2);
    return output;
}

//用于生成正态分布的随机数
dtype MoroInvCND(dtype P){
    const dtype a1 = 2.50662823884;
    const dtype a2 = -18.61500062529;
    const dtype a3 = 41.39119773534;
    const dtype a4 = -25.44106049637;
    const dtype b1 = -8.4735109309;
    const dtype b2 = 23.08336743743;
    const dtype b3 = -21.06224101826;
    const dtype b4 = 3.13082909833;
    const dtype c1 = 0.337475482272615;
    const dtype c2 = 0.976169019091719;
    const dtype c3 = 0.160797971491821;
    const dtype c4 = 2.76438810333863E-02;
    const dtype c5 = 3.8405729373609E-03;
    const dtype c6 = 3.951896511919E-04;
    const dtype c7 = 3.21767881768E-05;
    const dtype c8 = 2.888167364E-07;
    const dtype c9 = 3.960315187E-07;
    dtype y, z;
    y = P - 0.5;  //偏移到中点为原点
    //Moro算法以0.08和0.92为分界点，中间(0.08<x<0.92)用Beasley算法，两边用截断的切比雪夫序列
    if(fabs(y) < 0.42){          
        z = y * y;
        z = y * (((a4 * z + a3) * z + a2) * z + a1) / ((((b4 * z + b3) * z + b2) * z + b1) * z + 1);  //嵌套乘法
    }
    else
    {
        if(y > 0)
            z = log(-log(1.0 - P));
        else
            z = log(-log(P));
        z = c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * (c6 + z * (c7 + z * (c8 + z * c9)))))));
        if(y < 0) z = -z;
    }

    return z;
}

//标准正态分布的累积函数
dtype normalCFD(dtype value)
{
    return 0.5*erfc(-value*M_SQRT2);
}

inline dtype max(dtype value1,dtype value2)
{
    if(value1>value2)
        return value1;
    else
        return value2;
}

inline dtype min(dtype value1,dtype value2)
{
    if(value1<value2)
        return value1;
    else
        return value2;
}

int poisson(int lambda)
//产生一个泊松分布的随机数，lamda为总体平均数
{
    int x=-1;
    dtype u;
    dtype log1,log2;
    log1=0;
    log2=-lambda;
    do
    {
        u=rand()/static_cast<dtype>(RAND_MAX-1);
        log1+=log(u);
        x++;
    }while(log1>=log2);
    return x;
}

//简单的最小二乘拟合，拟合多项式为一阶多项式（y=ax+b）
class SimpleLeastSquare
{
    public:
        dtype a,b;
    public:
        SimpleLeastSquare(const dtype *x,const dtype *y, const int length)
        {
            dtype s1=0,s2=0,s3=0,s4=0;
            for(int i=0;i<length;i++)
            {
                s1+=x[i]*y[i];
                s2+=x[i];
                s3+=y[i];
                s4+=x[i]*x[i];
            }
            a=(s1*length-s2*s3)/(s4*length-s2*s2);
            b=(s4*s3-s2*s1)/(s4*length-s2*s2);

            return ;
        }
        dtype get_Y(const dtype x)
        {
            return a*x+b;
        }

    ////测试该类的代码如下（正确拟合系数为a:-0.9，b:0.2）：
    //dtype x[5]={-1.0,0.0,1.0, 2.0};
    //dtype y[5]={ 1.0,0.0,0.0,-2.0};
    //SimpleLeastSquare ls(x,y,4);
    //cout<<"a: "<<ls.a<<"   b: "<<ls.b<<endl;
    //cout<<ls.get_Y(8)<<endl;
};









inline dtype laguerre(dtype x, int order)
{
    //cout<<x<<" "<<order<<endl;
    dtype LaguerrePolynome =0;
    if(order==0)     
    {
        LaguerrePolynome=1;

    }
    else 
    {
        if(order==1)     
        {
            LaguerrePolynome=(1-x);
        }
        else 
        {
            if(order==2) 
                LaguerrePolynome=(1-2*x+0.5*x*x);
        }
    }
    
    return LaguerrePolynome;
    //return (LaguerrePolynome*exp(-x/2)); //这个指数部分不能带
}




//利用Eigen库实现 任意阶多项式 最小二乘拟合（默认二阶）
class LeastSquareEigen{
    public:
        dtype *coef; //1.系数个数随方程阶数变动
        int order; //2.多项式阶数
    public:
        LeastSquareEigen(
                dtype *x,  //自变量值
                dtype *y,  //因变量值，要使用Eigen的Map，不能加const关键字
                const int length, //数据点个数
                const int order=2  //多项式阶数，默认为二阶多项式
                )
        {
            this->order=order;
            int row=length; //4
            int col=order+1; //3
            //1.先要把我们的C++数组装进eigen的矩阵和向量中去
            //2.使用householder-QR方法求解最小二乘解
            coef=new dtype[col];
            dtype *A=new dtype[row*col];

            for(int i=0;i<col;i++)
            {
                for(int j=0;j<row;j++)
                {
                    //A[i*row+j]=pow(x[j],i);
                    A[i*row+j]=laguerre(x[j],i);
                }
            }
            //cout<<row<<" "<<col<<endl;
            //for(int i=0;i<row;i++)
            //{
            //    for(int j=0;j<col;j++)
            //    {
            //        cout<<A[j*row+i]<<"   ";
            //    }
            //    cout<<endl;
            //}
            Map<MatrixXd> mat_A(A,row,col);
            //cout<<mat_A<<endl;
            Map<MatrixXd> vec_y(y,length,1);
            //cout<<vec_y<<endl;

            //VectorXd out=mat_A.HouseholderQr().solve(vec_y); //快速但不稳定
            //VectorXd out=mat_A.colPivHouseholderQr().solve(vec_y); //速度和稳定性折中
            VectorXd out=mat_A.fullPivHouseholderQr().solve(vec_y); //慢但稳定

            //cout<<out<<endl;
            //cout<<"参数为"<<endl;
            for(int i=0;i<col;i++)
            {
                coef[i]=out[i];
                //cout<<coef[i]<<",";
            }
            //cout<<endl;
            
            return;
        }

        ~LeastSquareEigen()
        {
            delete[] coef;
        }

    public:
        dtype get_Y(const dtype x)
        {
            dtype sum=0.0;
            int number_of_items=order+1;
            for(int i=0;i<number_of_items;i++)
            {
                sum+=coef[i]*laguerre(x,i);
                //cout<<laguerre(x,i)<<" "<<coef[i]<<endl;
                //sum+=coef[i]*pow(x,i);
            }
            return sum;
        }

};

//2020.11.8
//为了对American option LSMC 程序debug，增加了对价格路径保存和加载的方法
//由于保存格式为csv，加载时涉及字符串切分，固自建了split函数
//来源： //https://blog.csdn.net/luoyexuge/article/details/81902346?utm_medium=distribute.pc_relevant_bbs_down.none-task-blog-baidujs-1.nonecase&depth_1-utm_source=distribute.pc_relevant_bbs_down.none-task-blog-baidujs-1.nonecase
vector<string> split(const string &s,const string &seperator)
{
    vector<string> result;
    typedef string::size_type string_size;
    string_size i=0;

    while( i!= s.size())
    {
        int flag=0;
        while(i!=s.size() && flag==0)
        {
            flag=1;
            for(string_size x=0;x<seperator.size();++x)
                if(s[i]==seperator[x])
                {
                    ++i;
                    flag=0;
                    break;
                }
        }

        flag=0;
        string_size j=i;
        while(j!=s.size()&&flag==0)
        {
            for(string_size x=0;x<seperator.size();++x)
                if(s[j]==seperator[x])
                {
                    flag=1;
                    break;
                }
            if(flag==0)
                ++j;
        }
        if(i!=j)
        {
            result.push_back(s.substr(i,j-i));
            i=j;
        }
    }
    return result;
}
void cholesky_decomposition(dtype *corr,dtype *L,int row,int col)
{
    Map<MatrixXd> mat_corr(corr,row,col);
    MatrixXd ml=mat_corr.llt().matrixL();
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
            L[i*col+j]=ml(i,j);
    }
    //cout<<ml<<endl;
    
    return ;
}
#endif
