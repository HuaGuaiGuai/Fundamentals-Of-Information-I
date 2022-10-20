import numpy as np
from pyldpc import make_ldpc, encode, decode, get_message
import math
import random as random
def p_random(arr1,arr2):
    assert len(arr1) == len(arr2), "Length does not match."
    assert sum(arr2) == 1 , "Total rate is not 1."

    sup_list = [len(str(i).split(".")[-1]) for i in arr2]
    top = 10 ** max(sup_list)
    new_rate = [int(i*top) for i in arr2]
    rate_arr = []
    for i in range(1,len(new_rate)+1):
        rate_arr.append(sum(new_rate[:i]))
    rand = random.randint(1,top)
    data = None
    for i in range(len(rate_arr)):
        if rand <= rate_arr[i]:
            data = arr1[i]
            break
    return data
F=0
A=1000
for ab in range(A):
    n =1024
    d_v = 3
    d_c = 8
    p =0.06
    H, G = make_ldpc(n, d_v, d_c, systematic=True, sparse=True)     #生成H和G
    k = H.shape[0]
    r=np.zeros([2], dtype=float)
    P=np.zeros([2],dtype=float)
    R=np.zeros([n],dtype=float)
    r[1]=1
    P[0]=1-p
    P[1]=p
    for i in range(n):
        R[i]=p_random(r,P)
    for i in range(n):
        if R[i]==0:
            R[i]=math.log((1-p)/(p),2)                #求取信道传输的r
        else:
            R[i]=math.log(p/(1-p),2)

    '''--------------------解码过程------------------'''
    M=np.zeros([k,n], dtype=float)
    E=np.zeros([k,n], dtype=float)
    '''--------------------迭代过程--------------------'''
    for i in range(20):           #最大迭代次数20次
        for j in range(n):
            for a in range(k):          #循环经过每个元素
                if H[a,j]!=0:             #如果遇到非0元素
                    X=E.sum(axis=0)         #先对E矩阵进行按列求和，获得一个一维矩阵
                    M[a,j]=X[j]-E[a,j]+R[j]     #对元素进行更新，H=E对应列求和-E+信道输入

        for a in range(k):
                Y=1
                for b in range(n):
                    if M[a, b] != 0:
                        Y = Y * math.tanh((M[a, b]) / 2)  # 对M矩阵进行处理
                for j in range(n):  # 循环经过每个元素
                    if M[a,j]!=0:
                        E[a,j]=2*math.atanh(Y/math.tanh(M[a,j]/2))            #更新结果

        Y=E.sum(axis=0)                     #将M矩阵按列求和进行译码

        for b in range(n):                      #如果大于0，则译为0，小于0，则译为1
            Y[b]=Y[b]+R[b]
            if Y[b]>0:
                Y[b]=0
            else:
                Y[b]=1
        N=0

        for b in range(n):                      #比对原始码字和译码码字，找出不同的个数
            if Y[b]!=0:
                N=N+1
        if N==0:
            print('迭代次数:',i)
            break

    print(N)
    if N>0:
        F=F+1
    F=float(F)
print(F/A)

