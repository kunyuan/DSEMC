import os
import sys
import re
import glob
import time
import numpy as np

SleepTime = 10

rs = None
Lambda = None
Beta = None
TotalStep = None
BetaStr = None
rsStr = None
LambdaStr = None

with open("inlist", "r") as file:
    line = file.readline()
    para = line.split(" ")
    BetaStr = para[0]
    Beta = float(BetaStr)
    rsStr = para[1]
    rs = float(rsStr)
    LambdaStr = para[2]
    Lambda = float(LambdaStr)
    TotalStep = float(para[4])

print rs, Beta, Lambda, TotalStep

# 0: I, 1: T, 2: U, 3: S
Channel = [1, ]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

folder = "./Beta{0}_rs{1}_lambda{2}/".format(BetaStr, rsStr, LambdaStr)

TauBin = None
AngleBin = None
ExtMomBin = None
AngleBinSize = None
TauBinSize = None
ExtMomBinSize = None
Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)

##############   2D    ##################################
###### Bare Green's function    #########################
# kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble = 0.15916/2  # 2D, Beta=10, rs=1
# Bubble = 0.0795775  # 2D, Beta=20, rs=1

#############  3D  ######################################
kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs
Bubble = 0.0971916  # 3D, Beta=10, rs=1
Step = None


def AngleIntegation(Data, l):
    # l: angular momentum
    shape = Data.shape[1:]
    Result = np.zeros(shape)
    for x in range(AngleBinSize):
        # Result += Data[x, ...] * \
        #     np.cos(l*AngleBin[x])*2.0*np.pi/AngleBinSize
        Result += Data[x, ...]*2.0/AngleBinSize
    return Result/2.0
    # return Result


def TauIntegration(Data):
    return np.sum(Data, axis=-1) * \
        Beta/kF**2/TauBinSize


while True:

    time.sleep(SleepTime)

    for order in Order:
        for chan in Channel:

            files = os.listdir(folder)
            Num = 0
            Norm = 0
            Data0 = None
            FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)

            for f in files:
                if re.match(FileName, f):
                    print "Loading ", f
                    Norm0 = -1
                    d = None
                    try:
                        with open(folder+f, "r") as file:
                            line0 = file.readline()
                            Step = int(line0.split(":")[-1])/1000000
                            # print "Step:", Step
                            line1 = file.readline()
                            # print line1
                            Norm0 = float(line1.split(":")[-1])
                            # print "Norm: ", Norm0
                            line2 = file.readline()
                            if TauBin is None:
                                TauBin = np.fromstring(
                                    line2.split(":")[1], sep=' ')
                                TauBinSize = len(TauBin)
                            line3 = file.readline()
                            if AngleBin is None:
                                AngleBin = np.fromstring(
                                    line3.split(":")[1], sep=' ')
                                AngleBinSize = len(AngleBin)
                            line4 = file.readline()
                            if ExtMomBin is None:
                                ExtMomBin = np.fromstring(
                                    line4.split(":")[1], sep=' ')
                                ExtMomBinSize = len(ExtMomBin)
                                ExtMomBin /= kF
                        # Num += 1
                        # print "Load data..."
                        d = np.loadtxt(folder+f)

                        if d is not None and Norm0 > 0:
                            if Data0 is None:
                                Data0 = d
                            else:
                                Data0 += d

                            Norm += Norm0

                    # print "Norm", Norm

                    except:
                        print "fail to load ", folder+f

            if Norm > 0 and Data0 is not None:
                print "Total Weight: ", Data0[0]
                Data0 /= Norm
                if(chan == 1):
                    Data0 = Data0.reshape(
                        (AngleBinSize, ExtMomBinSize, TauBinSize))
                else:
                    Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize))

                if DataWithAngle.has_key((order, chan)):
                    DataWithAngle[(order, chan)] = DataWithAngle[(
                        order, chan)]*0.0+Data0*1.0
                else:
                    DataWithAngle[(order, chan)] = Data0
                Data[(order, chan)] = AngleIntegation(
                    DataWithAngle[(order, chan)], 0)

    if len(DataWithAngle) > 0:
        print "Write Weight file."
        with open("weight1.data", "w") as file:
            for angle in range(AngleBinSize):
                for qidx in range(ExtMomBinSize):
                    for tidx in range(TauBinSize):
                        file.write("{0} ".format(
                            DataWithAngle[(0, 1)][angle, qidx, tidx]))

        qData = np.sum(Data[(0, 1)], axis=1) * \
            Beta/kF**2/TauBinSize

        qData = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData

        with open("data.data", "a") as file:
            file.write("{0}\n".format(qData[0]))

        print qData

    if Step >= TotalStep:
        print "End of Simulation!"
        sys.exit(0)
