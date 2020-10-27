from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue
import numpy as np

from array import array
import math
import sys,os

fd = TFile(sys.argv[1])
datatype = sys.argv[2]
tag=sys.argv[3]
base=os.path.basename(fd.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_hadron_kin_'+tag+'.pdf'
can = TCanvas('can','can',1200,3200)
can.Print('{}['.format(mon_out_file_name))
can.Divide(3,8)




