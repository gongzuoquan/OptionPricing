#ifndef UTIL_H_
#define UTIL_H_

#include<vector>
#include<string>
#include<iostream>
using namespace std;
#define dtype double
#define SND_NUM 0.3989422804
//工具类，用于记录一些实用函数

extern dtype gaussianBoxMuller();

extern dtype gaussianPolarMethod();

//用于同时返回两个返回值（gaussianBoxMuller2，gaussianPolarMethod2）
//其实使用引用参数也可以
struct mytuple
{
    dtype u;
    dtype v;
};

extern mytuple gaussianBoxMuller2();

extern mytuple gaussianPolarMethod2();

//标准正态分布
extern dtype std_normal(dtype input);

//用于生成正态分布的随机数
extern dtype MoroInvCND(dtype P);

//标准正态分布的累积函数
extern dtype normalCFD(dtype value);

/*
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
*/

extern int poisson(int lambda);
//产生一个泊松分布的随机数，lamda为总体平均数

//简单的最小二乘拟合，拟合多项式为一阶多项式（y=ax+b）
class SimpleLeastSquare
{
    public:
        dtype a,b;
    public:
        SimpleLeastSquare(const dtype *x,const dtype *y, const int length);
        dtype get_Y(const dtype x);

    ////测试该类的代码如下（正确拟合系数为a:-0.9，b:0.2）：
    //dtype x[5]={-1.0,0.0,1.0, 2.0};
    //dtype y[5]={ 1.0,0.0,0.0,-2.0};
    //SimpleLeastSquare ls(x,y,4);
    //cout<<"a: "<<ls.a<<"   b: "<<ls.b<<endl;
    //cout<<ls.get_Y(8)<<endl;
};









inline dtype laguerre(dtype x, int order);




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
                );
        
        ~LeastSquareEigen();
        
    public:
        dtype get_Y(const dtype x);
        
};

//2020.11.8
//为了对American option LSMC 程序debug，增加了对价格路径保存和加载的方法
//由于保存格式为csv，加载时涉及字符串切分，固自建了split函数
//来源： //https://blog.csdn.net/luoyexuge/article/details/81902346?utm_medium=distribute.pc_relevant_bbs_down.none-task-blog-baidujs-1.nonecase&depth_1-utm_source=distribute.pc_relevant_bbs_down.none-task-blog-baidujs-1.nonecase
extern vector<string> split(const string &s,const string &seperator);

#endif
