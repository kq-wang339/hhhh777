from numpy import *
import matplotlib.pyplot as plt
from io import BytesIO

data_am15 = open("AM15G.dat",'r', encoding='utf-8') # solar spectrum
data_alpha = open("absorption.dat", 'r',encoding='utf-8') #光吸收系数
Eg = 1.40                              #带隙，单位eV
L = 3                                 #输入厚度，单位为μm
f = 1                               #f，直接带隙为1，间接带隙需要修改
L_max = 3                           #厚度最大值，单位微米
T = 300                             #温度，单位K
e = ev = 1.60217648740E-19
h = 6.626068E-34
c = 299792458
h_ev = h / ev
c_nm = c * 1E9
Pi = 3.1415926
k = 1.3806505E-23                   #单位J/K
k_ev = k / ev
Pin = 1000                          #太阳光输入效率


#将太阳光谱数据转换成二维列表
am15 = []
for line in data_am15:
    s = line.strip().split('\t')
    s1 = ' '.join(s) + '\n'
    s2 = s1.split()
    if s2 != []:
        am15.append([s2[0],s2[1]])
data_am15.close()

#将光吸收系数变为二维列表
alpha = []
for line in data_alpha:
    s = line.strip().split('\t')
    s1 = ' '.join(s) + '\n'
    s2 = s1.split()
    if s2 != []:
        alpha.append([float(s2[0]), float(s2[1])])


# preparing the data for calculating slme
# 差值过程，思路就是将光吸收与SLME的横坐标对标
data_in = []

for l in range(1, len(am15)) :          #am15为太阳光谱
#    x = am15[l].split()
    hv  = float(am15[l][0])                   #波长，nm
    nhv = float(am15[l][1])                   #入射能量
    for ll in range(len(alpha)-1) :
        if alpha[ll][0] <= hv and alpha[ll+1][0] >= hv :
            fact = (hv - alpha[ll][0])/(alpha[ll+1][0] - alpha[ll][0])
            tmp1 = alpha[ll][1]*(1-fact) + fact*alpha[ll+1][1]
            data_in.append([hv, nhv, tmp1])
            #数据内容分别为波长，太阳光入射能量，tmp1为光吸收系数
            break

dat = open('data_in_1.dat','w',encoding='utf-8')
for i in range(len(data_in)):
    string = str(data_in[i][0]) + '\t' + str(data_in[i][1]) + '\t' + str(data_in[i][2]) +'\n'
#   print(string)
    dat.write(string)
dat.close()

def get_I(l,f=1,data_in=data_in):
#产生短路电流和暗电流的函数，需要修改的参数有：l，厚度，单位微米；f，直接带隙为1，间接带隙需要修改

    Isc = 0.0
    I0 = 0.0
    L = l * 1E-4  # 厚度，单位微米

    for l in range(len(data_in) - 1):
        hv0 = data_in[l][0]  # 积分单元矩阵左横坐标
        hv1 = data_in[l + 1][0]  # 积分单元矩阵右横坐标
        #
        des1 = hv1 - hv0
        #
        aE0 = 1.0 - exp(-2.0 * L * data_in[l][2])
        aE1 = 1.0 - exp(-2.0 * L * data_in[l + 1][2])

        is0 = data_in[l][1] * (hv0 / h / c_nm) * aE0
        is1 = data_in[l + 1][1] * (hv1 / h / c_nm) * aE1

        Isc = Isc + e * (is0 + is1) * des1 / 2.0

        hv_0 = 1240 / hv0
        hv_1 = 1240 / hv1
        des2 = hv_0 - hv_1

        irb0 = 2 * Pi * hv_0 ** 2 / h_ev ** 3 / c ** 2 * (exp(-1 * hv_0 / k_ev / T)) * aE0
        irb1 = 2 * Pi * hv_1 ** 2 / h_ev ** 3 / c ** 2 * (exp(-1 * hv_1 / k_ev / T)) * aE1

        I0 = I0 + e * Pi / f * (irb0 + irb1) * des2 / 2.0

    return Isc, I0

def get_JVcurve(Isc, I0, Eg):
#产生JV曲线的函数，需要用到get_I输出的参数，Eg为带隙，单位为eV
    I = []
    V = []
    npts = int(Eg / 0.001)
    for ll in range(npts):
        Vap = ll * 0.001
        i = Isc - I0 * (exp(Vap / k_ev / T) - 1)
        #    print(I)
        I.append(i)
        V.append(Vap)
        if i <= 0:
            break

    plt.plot(V,I,'r', label='J-V curve')
    plt.ylim(0,Isc+50)              # xlim、ylim：分别设置X、Y轴的显示范围
    plt.xlim(0,Vap+0.05)
    plt.title("JV curve")           # title：设置子图的标题
    plt.savefig('JV-curve.png')
    plt.show()

    dat = open('JV-curve.dat', 'w', encoding='utf-8')
    for i in range(len(I)):
        string = str(V[i]) + '\t' + str(I[i]) + '\t' + str(I[i]*V[i]) +'\n'
 #       print(string)
        dat.write(string)
    dat.close()

    print('JV-curve中的信息：')
    print('开路电压 = ', Vap)
    print('短路电流 = ', Isc)
    print('SLME =' + str(get_slme(Isc,I0)) + '\t' + '厚度 = ' + str(L) + 'μm')

    return 0

def get_slme(Isc,I0):
#计算SLME的函数，会同时打印出短路电流，开路电压和SLME数据，需要用到get_I的输出参数
    npts = int(Eg / 0.001)
    maxIV = 0
    IVtmp = 0
    for ll in range(npts):
        Vap = ll * 0.001
        I = Isc - I0 * (exp(Vap / k_ev / T) - 1)
        IVtmp = Vap * I
        #    print(I)
        if IVtmp >= maxIV:
            maxIV = IVtmp
        elif I <= 0:
            break
#    print("短路电流 = ", Isc, "A/m2")
#    print("开路电压 = ", Vap, "V")
    slme = maxIV / Pin
#    print("SLME = ", slme)
    return slme


#主函数部分
#第一部分是画给定厚度的JV曲线，同时给出开路电压，短路电流和SLME


Isc,I0 = get_I(l=L, f=f)
#print(I0)
get_JVcurve(Isc, I0, Eg)
get_slme(Isc,I0)


#第二部分是画SLME随厚度变化曲线，需要输入曲线中厚度最大值和曲线撒点数

n = 100                  #曲线撒的点

npts = int(L_max*n)
Y = []
X = []
dat = open('SLME-curve.dat', 'w', encoding='utf-8')
slme = 0
slme_max = 0
for i in range(npts+1):
    l = i / n
    Isc, I0 = get_I(l=l)
#    print("厚度 =", l,"μm")
    slme = get_slme(Isc, I0)
    Y.append(slme)
    X.append(l)
    dat.write(str(l) + '\t' + str(slme) + '\n')
    if slme >= slme_max:
        slme_max = slme
        l_max = l
dat.close()
print('SLME-curve内信息：')
print('SLME_max = ' + str(slme_max) + '\t' + '厚度 = ' + str(l_max) + 'μm')

plt.plot(X,Y)
plt.ylim(0,Y[-1]+0.025)              # xlim、ylim：分别设置X、Y轴的显示范围
plt.xlim(0,L_max)
plt.title("SLME curve")           # title：设置子图的标题
plt.savefig('SLME-curve.png')
plt.show()

