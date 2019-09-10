import os
import sys
import re
import glob
import time
import numpy as np

rs = 1.0
Lambda = 2.0
Beta = 20
SleepTime = 5
TotalStep = 101

# 0: I, 1: T, 2: U, 3: S
Channel = [1, ]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

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
# kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble = 0.15916/2  # 2D, Beta=10, rs=1
# Bubble = 0.0795775  # 2D, Beta=20, rs=1

#############  3D  ######################################
kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs
Bubble = 0.0971916  # 3D, Beta=10, rs=1
Step = None

while True:

    time.sleep(SleepTime)

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
                    try:
                        with open(folder+f, "r") as file:
                            line1 = file.readline()
                            Step = int(line1.split(":")[-1])/1000000

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
                        Num += 1
                        d = np.loadtxt(folder+f)
                        if Data0 is None:
                            Data0 = d
                        else:
                            Data0 += d
                    except:
                        print "fail to load ", folder+f

            if Num > 0 and Data0 is not None:
                Data0 /= Num
                if(chan == 1):
                    Data0 = Data0.reshape(
                        (AngleBinSize, ExtMomBinSize, TauBinSize))
                else:
                    Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize))

                DataWithAngle[(order, chan)] = Data0

            # average the angle distribution
            # Data[(order, chan)] = AngleIntegation(Data0, 0)
    if len(DataWithAngle) > 0:
        print "Write Weight file."
        with open(folder+"weight1.dat", "w") as file:
            for angle in range(AngleBinSize):
                for qidx in range(ExtMomBinSize):
                    for tidx in range(TauBinSize):
                        file.write("{0} ".format(
                            DataWithAngle[0, 1][angle, qidx, tidx]))

    if Step >= TotalStep:
        print "End of Simulation!"
        sys.exit(0)
