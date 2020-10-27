from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import kRed
from ROOT import Fit
from ROOT import Math
from ROOT import TMath
from ROOT import TRandom
from ROOT import cout
from array import array
import math
import sys

from scipy.optimize import minimize
import numpy as np
    

class NormFunct(TF1):
    def __init__(self, func_to_wrap):
        self.myfunct = func_to_wrap
   
    def operator(self, x, p):
        xmin, xmax = 0,0 #(or None, None)

        myfunct.GetRange(xmin,xmax)
        return myfunct.EvalPar(x,p)/myfunct.Integral(xmin, xmax)    




r_num = TF1("r_num","gaus",-1,1)
r_num.SetParameter(0,0)
r_num.SetParameter(1,0.1)

#data = [1,2,2.1,3.3,3.1,3.0]
data = []
random = TRandom()
for ii in range(0,361):
    y = 3*TMath.Sin(ii * 3.1415926/180.0) + random.Gaus(0,0.6)
    print(y)
    data.append(y)


opt = Fit.DataOptions()
d_range = Fit.DataRange(0,360)
unbin = Fit.UnBinData(len(data))

for ii in data:
    unbin.Add(ii)

fasy=TF1("fasy","[0]*sin(x*3.1415926/180.0)",0.0,360.0)
fasy.SetParameter(0,2)
#my_test_f = TF1("test",NormFunct(f), 0 , 10)

norm_fasy = TF1("norm_fasy",NormFunct(fasy),0.0,360.0)
# ROOT::Math::normal_pdf(x,[1],[0])", 0.0, 5.0) #sin(x*3.1415926/180.0)",0.0,361.0)

#fasy.SetParameter(1,1)

fit_funct = Math.WrappedMultiTF1(norm_fasy, norm_fasy.GetNdim())
fitter = Fit.Fitter()

fitter.SetFunction(fit_funct,False)
#fitter.Config().ParSettings(0).SetLimits(0,10000)
#fitter.Config().ParSettings(1).SetLimits(0,10000)
#fitter.Config().SetUpdateAfterFit()

fitter.LikelihoodFit(unbin)


my_fit_result = fitter.Result()
my_fit_result.Print(cout)


ctest = TCanvas("ctest","ctest",900,900)
ctest.cd()
fit_pars = my_fit_result.GetParams()
print(fit_pars)
fasy.SetParameter(0, fit_pars[0])
#fasy.SetParameter(1, fit_pars[1])
fasy.Draw()
ctest.SaveAs("my_likelihoodtest.pdf")

