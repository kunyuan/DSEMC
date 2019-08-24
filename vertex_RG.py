import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob
import os
import re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

rs = 1.0
Lambda = 2.0
Beta = 20
# XType = "Tau"
XType = "Mom"
# XType = "Angle"

OrderByOrder = False

# 0: I, 1: T, 2: U, 3: S
Channel = [1, ]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, 1, 2, 3]

folder = "./Beta{0}_rs{1}_lambda{2}/".format(Beta, rs, Lambda)

Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)
TauBin = None
AngleBin = None
ExtMomBin = None
AngleBinSize = None
TauBinSize = None
ExtMomBinSize = None

##############   2D    ##################################
###### Bare Green's function    #########################
kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble = 0.15916/2  # 2D, Beta=10, rs=1
Bubble = 0.0795775  # 2D, Beta=20, rs=1


def AngleIntegation(Data, l):
    # l: angular momentum
    shape = Data.shape[1:]
    Result = np.zeros(shape)
    for x in range(AngleBinSize):
        Result += Data[x, ...] * \
            np.cos(l*AngleBin[x])*2.0*np.pi/AngleBinSize
    return Result/2.0/np.pi


for order in Order:
    for chan in Channel:

        files = os.listdir(folder)
        Num = 0
        Data0 = None
        if(order == 0):
            FileName = "vertex{0}_pid[0-9]+.dat".format(chan)
        else:
            FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)

        for f in files:
            if re.match(FileName, f):
                print "Loading ", f
                with open(folder+f, "r") as file:
                    line1 = file.readline()
                    line2 = file.readline()
                    if TauBin is None:
                        TauBin = np.fromstring(line2.split(":")[1], sep=' ')
                        TauBinSize = len(TauBin)
                    line3 = file.readline()
                    if AngleBin is None:
                        AngleBin = np.fromstring(line3.split(":")[1], sep=' ')
                        AngleBinSize = len(AngleBin)
                    line4 = file.readline()
                    if ExtMomBin is None:
                        ExtMomBin = np.fromstring(line4.split(":")[1], sep=' ')
                        ExtMomBinSize = len(ExtMomBin)
                        ExtMomBin /= kF
                Num += 1
                d = np.loadtxt(folder+f)
                if Data0 is None:
                    Data0 = d
                else:
                    Data0 += d
        Data0 /= Num
        if(chan == 1):
            Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize, TauBinSize))
        else:
            Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize))

        DataWithAngle[(order, chan)] = Data0

        # average the angle distribution
        Data[(order, chan)] = AngleIntegation(Data0, 0)


def ErrorPlot(p, x, d, color, marker, label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


w = 1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy', 'y']
ColorList = ColorList*40

if(XType == "Scale"):
    for i in range(ExtMomBinSize/4):
        index = 4*i
        ErrorPlot(ax, ScaleBin[:-2], diffData[:-2, index],
                  ColorList[i], 's', "Q {0}".format(ExtMomBin[index]))
    ax.set_xlim([0.0, ScaleBin[-2]])
    ax.set_xlabel("$Scale$", size=size)
elif(XType == "Tau"):
    chan = 1
    for i in range(ExtMomBinSize/8):
        index = 8*i
        for order in Order:
            ErrorPlot(ax, TauBin, Data[(order, chan)][index, :],
                      ColorList[2*i], 's', "Loop {0}, Q {1}".format(order, ExtMomBin[index]))
    ax.set_xlim([0.0, TauBin[-1]+0.1])
    ax.set_xlabel("$Tau$", size=size)
elif (XType == "Mom"):
    i = 0
    for order in Order[1:]:
        for chan in Channel:
            i += 1
            if(chan == 1):
                qData = np.sum(Data[(order, chan)], axis=1) * \
                    Beta/kF**2/TauBinSize
            else:
                qData = Data[(order, chan)]

            # qData = np.sum(qData, axis=1)*Beta/kF**2/TauBinSize
            # qData0 = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData0
            # qData=8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData
            ErrorPlot(ax, ExtMomBin, qData,
                      ColorList[i], 's', "Loop {0}, Chan {1}".format(order, ChanName[chan]))

    i = 0

    for chan in Channel:
        i += 1
        if(chan == 1):
            qData = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)
            qData -= np.sum(Data[(0, chan)], axis=1) * \
                Beta/kF**2/TauBinSize
            # qData -= np.sum(Data[(1, chan)], axis=1) * \
            #     Beta/kF**2/TauBinSize
            # qData -= np.sum(Data[(2, chan)], axis=1) * \
            #     Beta/kF**2/TauBinSize
            # qData -= np.sum(Data[(3, chan)], axis=1) * \
            #     Beta/kF**2/TauBinSize
        else:
            qData = Data[(0, chan)]

        # qData = np.sum(qData, axis=1)*Beta/kF**2/TauBinSize
        # qData0 = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData0
        # qData=8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData
        ErrorPlot(ax, ExtMomBin, qData,
                  ColorList[i], 'o', "Chan {1}".format(0, ChanName[chan]))

    x = np.arange(0, 3.0, 0.001)
    y = x*0.0+Bubble
    for i in range(len(x)):
        if x[i] > 2.0:
            y[i] = Bubble*(1-np.sqrt(1-4/x[i]**2))
    y0 = 8.0*np.pi/(x*x*kF*kF+Lambda)
    # ym=y0-y0*y0*y
    yphy = 8.0*np.pi/(x*x*kF*kF+Lambda+y*8.0*np.pi)

    # ax.plot(x, yphy, 'k-', lw=2, label="physical")
    ax.plot(x, y0, 'k-', lw=2, label="original")

    # ax.plot(x, y0*y0*y, 'r-', lw=2, label="wrong")

    ax.set_xlim([0.0, ExtMomBin[-1]])
    # ax.set_xscale("log")
    ax.set_xlabel("$q/k_F$", size=size)

elif(XType == "Angle"):
    for i in range(ScaleBinSize/8):
        # print i, index
        # print ScaleBin[index]
        index = 8*i
        ErrorPlot(ax, AngleBin, Data[index, :, 5],
                  ColorList[i], 's', "Q {0}".format(ScaleBin[index]))
    ax.set_xlim([0.0, AngleBin[-1]])
    ax.set_xlabel("$Angle$", size=size)
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-\Gamma_4(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
