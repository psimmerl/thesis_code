from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from collections import OrderedDict

import numpy as np
import os

from array import array
import math
import sys

ff = TFile(sys.argv[1])
datatype = sys.argv[2]
tag=sys.argv[3]

base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

def compute_if_absent(d,k,default):
    if k in d:
        return d[k]
    else:
        d[k] = default(k)
        return d[k]

histos = OrderedDict()# {}

def buildP(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,0.0,beam.E())
def buildTheta(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,0.0,60)
def buildPTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title, 200, 0.0, beam.E(), 200, 0.0, 60.0)
def buildPhiTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,200,-180.0,180.0,200, 0.0, 60.0)

builders={
    'p':buildP,
    'theta':buildTheta,
    'ptheta':buildPTheta,
    'phitheta':buildPhiTheta,
}
