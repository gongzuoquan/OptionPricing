from sys import argv
_,path=argv
count=0
simN=0
params=''
with open(path,"r") as fin,open(path+".csv","w",encoding="gbk") as fout:
    fout.write("模拟规模,期权种类,期权参数,期权价格,期权误差,95%误差上限,95%误差下限,运行时间(s)\n")
    sim_write_flag=True
    param_write_flag=False
    for line in fin.readlines():
        if count>3 :
            fout.write("\n")
            count=0
        line=line.strip()

        if param_write_flag and line=="":
            param_write_flag=False

        if line=="":
            continue

        if line[0]=="@":
            simN=line.split(" ")[1]
            sim_write_flag=True
            fout.write(str(simN)+",")  

        if line[:2]=="S0":
            param_write_flag=True
        if param_write_flag:
            params+=line+';'
        
        if line[0]==">":
            if count==0:
                if sim_write_flag:
                    sim_write_flag=False
                else:
                    fout.write(",")
                option_type=line.split("价格为")[0][1:]
                fout.write(option_type+",")
                fout.write(params+',')
                params=''

            if count==2:
                line=line.split("[")[1][:-1].strip()
                borders=line.split(";")
                up=borders[0].strip()
                down=borders[1].strip()
                fout.write(up+",")
                fout.write(down+",")
            else:
                line=line.split("：")[1]
                fout.write(line+",")
            count+=1


