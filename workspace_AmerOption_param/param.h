#ifndef OPTION_PARAM_H
#define OPTION_PARAM_H
//参数类
class Params
{
    public:
        dtype S0,K,T,V,R;
        int sim_num,time_step;
        dtype rho,kappa,theta,epslon;
    public:
        Params(string param_file)
        {
            ifstream fin(param_file,ios::in);
            string temp;
            int param_num=0;
            while(getline(fin,temp))
            {
                vector<string> result=split(temp,"=");
                if(result[0]=="S0")
                    S0=stod(result[1]);
                else if(result[0]=="K")
                    K=stod(result[1]);
                else if(result[0]=="T")
                    T=stod(result[1]);
                else if(result[0]=="R")
                    R=stod(result[1]);
                else if(result[0]=="V")
                    V=stod(result[1]);
                else if(result[0]=="sim_num")
                    sim_num=stoi(result[1]);
                else if(result[0]=="time_step")
                    time_step=stoi(result[1]);
                else if(result[0]=="rho")
                    rho=stod(result[1]);
                else if(result[0]=="kappa")
                    kappa=stod(result[1]);
                else if(result[0]=="theta")
                    theta=stod(result[1]);
                else if(result[0]=="epslon")
                    epslon=stod(result[1]);
                param_num++;
            }
            fin.close();
        }

        void display()
        {
            cout<<"参数列表如下："<<endl;
            cout<<"S0: "<<S0<<endl;
            cout<<"K: "<<K<<endl;
            cout<<"T: "<<T<<endl;
            cout<<"V: "<<V<<endl;
            cout<<"R: "<<R<<endl;
            cout<<"rho: "<<rho<<endl;
            cout<<"kappa: "<<kappa<<endl;
            cout<<"theta: "<<theta<<endl;
            cout<<"epslon: "<<epslon<<endl;
            cout<<endl;
        }

};


