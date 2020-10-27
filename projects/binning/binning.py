from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
import numpy as np

from array import array
import math
import sys

ff = TFile(sys.argv[1])
datatype = sys.argv[2]

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

for ev in ff.clas12:
    el = TLorentzVector(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)

    lmda = pro+kp 
    vphi = kp+km 
    q = beam-el  
    q2 = -q.M2()  
    xb = q2/(2*target.M()*q.E())
    t = 2*target.M()*(pro.E()-target.M()) 
    nu = q2/(2*target.M()*xb)
    y = nu/beam.E() 

    v3l = beam.Vect().Cross(el.Vect())
    v3h = pro.Vect().Cross((beam-el).Vect()) 
    trento = np.degrees(TMath.ACos(v3l.Dot(v3h)/(v3l.Mag()*v3h.Mag())))
    if(v3l.Dot(pro.Vect())<0): trento=360-trento

    boost_phicm = vphi.BoostVector()
    boost_phicm = -boost_phicm
    boost_kp = TLorentzVector(kp)
    boost_km = TLorentzVector(km)
    boost_kp.Boost(boost_phicm)
    boost_km.Boost(boost_phicm)
    cos_theta_cm_kp = TMath.Cos(boost_kp.Vect().Theta())

    
    
