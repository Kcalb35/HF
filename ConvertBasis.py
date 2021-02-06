import re

def process(li,number,file):
    if len(li)>0:
        name = split(li[0][0])[0]
        file.write(str(number)+" "+name)
        
        s =[]
        p =[]
        d=[]
        for shell in li:
            obital = split(shell[0])[1]
            if (obital=="S"):
                tmp = []
                for line in shell[1:]:
                    componet = [float(x) for x in split(line)]
                    tmp.append(componet)
                s.append(tmp)
            elif (obital=="SP"):
                tmps = []
                tmpp =[]
                for line in shell[1:]:
                    component = split(line)
                    tmps.append([float(component[0]),float(component[1])])
                    tmpp.append([float(component[0]),float(component[2])])
                s.append(tmps)
                p.append(tmpp)
            elif(obital=="D"):
                raise NotImplementedError
        
        file.write(" "+str(len(s)+len(p)+len(d))+"\n")

        print(s)
        print(p)
        file.write(seralizeOneL(s,0))
        file.write(seralizeOneL(p,1))



def seralize(li):
    li = [str(x) for x in li]
    return " ".join(li) 

def seralizeOneL(li,l):
    s = ""
    n = l
    for ob in li:
        n +=1
        total = len(ob)
        s += str(n) + " "+str(l)+" "+str(total)+"\n"
        for gto in ob:
            s += seralize(gto)+"\n"
    return s

def split(s):
   return [ x for x in s.split(" ") if x]


with open("baseset/STO-3G.txt") as f:
    with open("STO-3G.txt",'w') as f2:
        line = f.readline()
        li =[]
        i = -1
        number=-1
        NumOfAtom = -1 
        while True:
            line = line.replace('\n','')
            if re.search("^#",line) is not None or line=="":
                number = number+1
                process(li,number,f2)
                NumOfAtom += 1
                li = []
                i=-1
                if line=="":
                    break 
            else :
                # add to buffer
                if re.search("^[A-z]",line) is not None :
                    i = i+1
                    li.append([])
                li[i].append(line)
            line = f.readline()