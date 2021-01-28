#ifndef TREE_H_
#define TREE_H_
#include"option.hpp"
#include<new>
#include<cmath>
#include<iomanip>
class BitreeOption:public Option
{
    public:
        int node_num;
        dtype *bitree;
        dtype up_range;
        dtype down_range;
    public:
        //无参数构造函数需复制到Option子类
        BitreeOption():Option()
        {
            node_num=0;
            bitree=NULL;
            up_range=0.0;
            down_range=0.0;
        }
        BitreeOption(
                const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V,
                const dtype in_up_range=0.0,
                const dtype in_down_range=0.0
                ):Option(in_S0,in_K,in_T,in_R,in_V)

        {
            node_num=0;
            bitree=NULL;
            up_range=in_up_range;
            down_range=in_down_range;
            return;
        }
        ~BitreeOption()
        {
            if(bitree==NULL)
                return;
            delete[] bitree;
            bitree=NULL;
        }
    public:
        /*
        void pricing(const int time_step,const bool is_input=false)
        {

            node_num=pow(2,time_step);
            bitree=new(nothrow) dtype[node_num];
            if(bitree==NULL)
            {
                cerr<<"存储空间不足，分配错误"<<endl;
                return;
            }
            bitree[0]=0;
            bitree[1]=S0;
            //计算股票价格
            for(int i=2;i<node_num;i++)
            {
                if(i%2==0)
                {
                    bitree[i]=bitree[i/2]*up_range;
                }
                else
                {
                    bitree[i]=bitree[i/2]*down_range;
                }
            }
            dtype payoff=0.0;
            dtype St_num=pow(2,time_step-1);
            for(int i=node_num-1;i>=pow(2,time_step-1);i--)
            {
                if(is_input)
                    payoff+=payoff_put(bitree[i]);
                else
                    payoff+=payoff_call(bitree[i]);
                cout<<payoff_put(bitree[i])<<endl;
            }
            payoff/=St_num;
            price=payoff;
            //delete[] bitree;
            return;
        }
        */
        
        void pricing(const int time_step,const bool is_input=false)
        {
            node_num=time_step+1;
            bitree=new(nothrow) dtype[node_num];
            if(bitree==NULL)
            {
                cerr<<"存储空间不足，分配错误"<<endl;
                return;
            }
            dtype dt=T/time_step;

            discount=exp(-R*dt);  //折现率
            dtype discount_inv=1.0/discount; //逆折现率 

            up_range=exp(V*sqrt(dt));
            down_range=1/up_range;
            dtype up_range_square=up_range*up_range;
            bitree[0]=S0*pow(down_range,time_step);

            dtype p=(discount_inv-down_range)/(up_range-down_range);
            dtype q=1-p;

            for(int i=1;i<node_num;i++)
            {
                bitree[i]=bitree[i-1]*up_range_square;
                //cout<<bitree[i]<<" "<<endl;
            }
            //cout<<endl;
            for(int i=0;i<node_num;i++)
            {
                if(is_input)
                    bitree[i]=payoff_put(bitree[i]);
                else
                    bitree[i]=payoff_call(bitree[i]);
            }
            for(int t=time_step-1;t>=0;t--)
            {
                for(int i=0;i<=t;i++)
                {
                    bitree[i]=(bitree[i]*q+bitree[i+1]*p)*discount;
                }
            }
            price=bitree[0];
            return;
        }

        /*
        void pricing(const int time_step,const bool is_input=false)
        {
            node_num=time_step+1;
            bitree=new(nothrow) dtype[node_num];
            if(bitree==NULL)
            {
                cerr<<"存储空间不足，分配错误"<<endl;
                return;
            }
            dtype dt=T/time_step;

            discount=exp(-R*dt);  //折现率
            dtype discount_inv=1.0/discount; //逆折现率 

            //up_range=exp(V*sqrt(dt));
            //down_range=1/up_range;
            //dtype up_range_square=up_range*up_range;
            //bitree[0]=S0*pow(down_range,time_step);

            dtype p=(discount_inv-down_range)/(up_range-down_range);
            dtype q=1-p;
            //cout<<discount_inv-down_range<<endl;
            //cout<<exp(R*dt)-down_range<<endl;

            for(int i=0;i<node_num;i++)
            {
                bitree[i]=S0*pow(up_range,node_num-1-i)*pow(down_range,i);
                //cout<<i<<" "<<bitree[i]<<endl;
            }
            for(int i=0;i<node_num;i++)
            {
                if(is_input)
                    bitree[i]=payoff_put(bitree[i]);
                else
                    bitree[i]=payoff_call(bitree[i]);
            }
            for(int t=time_step-1;t>=0;t--)
            {
                for(int i=0;i<=t;i++)
                {
                    bitree[i]=(bitree[i+1]*q+bitree[i]*p)*discount;
                }
            }

            price=bitree[0];
            return;
        }
        */
};

#endif
