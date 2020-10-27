from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
import numpy as np

from array import array
import math
import sys

ff = TFile(sys.argv[1])

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

for ev in ff.clas12:
    el = TLorentzVector(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)
    phimass = (kp+km).M()

    hel = ev.hel

    v3l = beam.Vect().Cross(el.Vect())
    v3h = pro.Vect().Cross((beam-el).Vect()) 
    trento = np.degrees(TMath.ACos(v3l.Dot(v3h)/(v3l.Mag()*v3h.Mag())))
    if(v3l.Dot(pro.Vect())<0): trento=360-trento

    # look at the BSA for different phi mass cuts
    # how is the background influencing the bsa
    for cc in range(0,10):
        if( phimass < 1.025 + cc*0.02 ):
            if( hel > 0 ):
                # fill pos hist
                
            else if( hel < 0 ):
                # fill neg hist

        



    


    
