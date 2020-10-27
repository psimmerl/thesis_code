from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend
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
mc_data = sys.argv[4]


base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution

sig=6.0
epkpkmxe_min = 0.0222 - 0.0262*sig
epkpkmxe_max = 0.0222 + 0.0262*sig

epkpX_min = 0.234 - 0.0348*sig
epkpX_max = 0.234 + 0.0348*sig

epkmX_min = 0.228 - 0.0284*sig
epkmX_max = 0.228 + 0.0284*sig

ekpkmX_min = 0.891 - 0.0285*sig
ekpkmX_max = 0.891 + 0.0285*sig

def passEvent( epkpkmX, epkpX, epkmX, ekpkmX ):
    exclusivity_cut_results = (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min , epkpX.M2() < epkpX_max and epkpX.M2() > epkpX_min, ekpkmX.M2() < ekpkmX_max and ekpkmX.M2() > ekpkmX_min, epkmX.M2() < epkmX_max and epkmX.M2() > epkmX_min)    
    return exclusivity_cut_results

def passEPKPKMX(epkpkmX):
    return ( (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min ))

def passEPKPX(epkpX):
    return epkpX.M2() < epkpX_max and epkpX.M2() > epkpX_min

def passEPKMX(epkmX):
    return epkmX.M2() < epkmX_max and epkmX.M2() > epkmX_min

def passEKPKMX(ekpkmX):
    ekpkmX.M2() < ekpkmX_max and ekpkmX.M2() > ekpkmX_min
            

def returnListResult(index_to_not_check, list_in):
    
    final_result=-1
    cc=0
    for ii in list_in:
        if cc != index_to_not_check:
            if( not ii ):
                return False
        cc+=1
    return True

def compute_if_absent(d,k,default):
    if k in d:
        return d[k]
    else:
        d[k] = default(k)
        return d[k]

histos = OrderedDict()

def buildP(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,0.0,beam.E())
def buildTheta(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,0.0,60)
def buildPhi(title):
    print('building {}'.format(title))
    return TH1F(title,title,400,-180.0,180.0)
def buildPTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title, 200, 0.0, beam.E(), 200, 0.0, 60.0)
def buildPPhi(title):
    print('building {}'.format(title))
    return TH2F(title,title, 300, 0.0, beam.E(), 300, -180.0, 180.0)
def buildPhiTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,200,-180.0, 180.0,200, 0.0, 60.0)
def buildMM2(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,-1.50,1.5)
def buildMM2Large(title):
    print('building {}'.format(title))
    return TH1F(title,title,1000,-15,15)
def buildMass(title):
    print('building {}'.format(title))
    return TH1F(title,title,150,0.8,1.5)
def buildMassLarge(title):
    print('building {}'.format(title))
    return TH1F(title,title,250,0.8,2.5)
def buildMass2D(title):
    print('building {}'.format(title))
    return TH2F(title,title,75,1.4,2.0,75,0.8,1.5)
def buildMass2DLarge(title):
    print('building {}'.format(title))
    return TH2F(title,title,1000,-5,5,1000,-5,5)
def buildMassKin(title):
    print('building {}'.format(title))
    return TH2F(title,title,350,0,beam.E(),350,-4,4)

def buildMMThetaEPX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, 0.60, 1.50)
def buildMMThetaEPKPX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250,0.0,60.0,250, -0.2, 0.6)
def buildMMThetaEPKMX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -0.2, 0.6)
def buildMMThetaEKPKMX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, 0.50, 1.30)

def buildMMpEPX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 5.0, 250, 0.60, 1.50)
def buildMMpEPKPX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 5.0, 250, -0.20, 0.60)
def buildMMpEPKMX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 5.0, 250, -0.20, 0.60)
def buildMMpEKPKMX(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 4.0, 250, 0.50, 1.30)

def buildDeltaPPrSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -0.5, 0.5)
def buildDeltaPPrBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -1, 1)

def buildDeltaPKpSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -0.5, 0.5)
def buildDeltaPKpMed(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -5, 5)
def buildDeltaPKpMed2(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -5, 5)
def buildDeltaPKpBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -1, 1)

def buildDeltaPKmSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -0.5, 0.5)
def buildDeltaPKmMed(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -5, 5)
def buildDeltaPKmMed2(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -5, 5)
def buildDeltaPKmBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -1, 1)

def buildDeltaThetaPrSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -10.0, 10.0)
def buildDeltaThetaPrMed(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -0.50, 0.50)
def buildDeltaThetaPrBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -10.0, 10.0)

def buildDeltaThetaKpSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -10.0, 10.0)
def buildDeltaThetaKpMed(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -0.50, 0.50)
def buildDeltaThetaKpMed2(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -5.0, 5.0)
def buildDeltaThetaKpBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -10.0, 10.0)
def buildDeltaThetaKmSmall(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 3.0, 250, -10.0, 10.0)
def buildDeltaThetaKmMed(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -0.50, 0.50)
def buildDeltaThetaKmMed2(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -5.0, 5.0)
def buildDeltaThetaKmBig(title):
    print('building {}'.format(title))
    return TH2F(title,title, 250, 0.0, 60.0, 250, -10.0, 10.0)

        

#################################################################################
#################################################################################
## resolution from simulation ( generated - reconstructed )
def buildMCResMntmEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResThetaEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)
def buildMCResPxEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPyEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPzEl(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPxEl1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.20)
def buildMCResPyEl1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPzEl1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)


def buildMCResMntmPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)
def buildMCResPxPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 5.0, 250, -0.2, 0.2)
def buildMCResPyPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 5.0, 250, -0.2, 0.2)
def buildMCResPzPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPxPr1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPyPr1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPzPr1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)




def buildMCResMntmKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)
def buildMCResPxKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0, 4.0, 250, -0.2, 0.2)
def buildMCResPyKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0, 4.0, 250, -0.2, 0.2)
def buildMCResPzKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPxKp1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPyKp1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.20, 0.20)
def buildMCResPzKp1D(title):
    print('building {}'.format(title))
    return TH1F(title,title, 250, -0.2, 0.2)


def buildMCResMntmKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)
def buildMCResPxKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0, 4.0, 250, -0.2, 0.2)
def buildMCResPyKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.2, 0.2)
def buildMCResPzKm(title):
    print('buildg {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildMCResPxKm1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPyKm1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)
def buildMCResPzKm1D(title):
    print('building {}'.format(title))
    return TH1F(title,title,250, -0.2, 0.2)


def buildQ2(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0,beam.E())
def buildx(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0.0,1.05)
def buildt(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0,5.0)

def buildQ2x(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,1,100,0.0,beam.E())
def buildQ2t(title):
    print('building {}'.format(title))
    return TH2F(title,title,100, 0, beam.E(), 100, 0.0, 6)
def buildtx(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,6,100,0.0,1.0)
def buildBetaP(title):
    print('building {}'.format(title))
    return TH2F(title,title, 200, 0.0, beam.E(), 200, 0.0, 1.25)

def buildxElTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0.0,1.05, 100,0.0, 60.0 )
def buildxPrTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0.0,1.05, 100, 0.0, 60.0)
def buildxKpTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0.0,1.05, 100, 0.0, 60.0)
def buildxKmTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0.0,1.05, 100, 0.0, 60.0)


def buildDeltaQ2(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, beam.E(), 250, -0.2, 0.2)
def buildDeltaXb(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 1.1, 250, -0.2, 0.2)
def buildDeltaT(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 6.0, 250, -0.2, 0.2)
def buildDeltaCosTheta(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -0.2, 0.2)

def buildXbMM2(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -0.25, 1.25, 250, 0.0, 1.1)


builders={
    'p':buildP,
    'theta':buildTheta,
    'phi':buildPhi,
    'ptheta':buildPTheta,
    'phitheta':buildPhiTheta,
    'pphi':buildPPhi,
    'mm2': buildMM2,
    'mass': buildMass,
    'mm2L': buildMM2Large,
    'massL': buildMassLarge,
    'mass2D': buildMass2D,
    'mass2DL': buildMass2DLarge,
    'masskin': buildMassKin,
    'mmthetaEPX':buildMMThetaEPX,
    'mmthetaEPKPX':buildMMThetaEPKPX,
    'mmthetaEPKMX':buildMMThetaEPKMX,
    'mmthetaEKPKMX':buildMMThetaEKPKMX,
    'mmpEPX':buildMMpEPX,
    'mmpEPKPX':buildMMpEPKPX,
    'mmpEPKMX':buildMMpEPKMX,
    'mmpEKPKMX':buildMMpEKPKMX,
    'delpPrSmall':buildDeltaPPrSmall,
    'delpPrBig':buildDeltaPPrBig,
    'delpKpSmall':buildDeltaPKpSmall,
    'delpKpMed':buildDeltaPKpMed,
    'delpKpMed2':buildDeltaPKpMed2,
    'delpKpBig':buildDeltaPKpBig,
    'delpKmSmall':buildDeltaPKmSmall,
    'delpKmMed':buildDeltaPKmMed,
    'delpKmMed2':buildDeltaPKmMed2,
    'delpKmBig':buildDeltaPKmBig,

    'delthetPrSmall':buildDeltaThetaPrSmall,
    'delthetPrMed':buildDeltaThetaPrMed,
    'delthetPrBig':buildDeltaThetaPrBig,

    'delthetKpSmall':buildDeltaThetaKpSmall,
    'delthetKpMed':buildDeltaThetaKpMed,
    'delthetKpMed2':buildDeltaThetaKpMed2,
    'delthetKpBig':buildDeltaThetaKpBig,

    'delthetKmSmall':buildDeltaThetaKmSmall,
    'delthetKmMed':buildDeltaThetaKmMed,
    'delthetKmMed2':buildDeltaThetaKmMed2,
    'delthetKmBig':buildDeltaThetaKmBig,
    ###############################################
    ###############################################
    'mcresmntmel': buildMCResMntmEl,
    'mcreselpx': buildMCResPxEl,
    'mcreselpy': buildMCResPyEl,
    'mcreselpz': buildMCResPzEl,
    'mcresthetael': buildMCResThetaEl,
    'mcresphiel': buildMCResPhiEl,
    'mcreselpx1d': buildMCResPxEl1D,
    'mcreselpy1d': buildMCResPyEl1D,
    'mcreselpz1d': buildMCResPzEl1D,

    'mcresmntmpr': buildMCResMntmPr,
    'mcresprpx': buildMCResPxPr,
    'mcresprpy': buildMCResPyPr,
    'mcresprpz': buildMCResPzPr,
    'mcresthetapr': buildMCResThetaPr,
    'mcresphipr': buildMCResPhiPr,
    'mcresprpx1d': buildMCResPxPr1D,
    'mcresprpy1d': buildMCResPyPr1D,
    'mcresprpz1d': buildMCResPzPr1D,

    'mcresmntmkp': buildMCResMntmKp,
    'mcreskppx': buildMCResPxKp,
    'mcreskppy': buildMCResPyKp,
    'mcreskppz': buildMCResPzKp,
    'mcresthetakp': buildMCResThetaKp,
    'mcresphikp': buildMCResPhiKp,
    'mcreskppx1d': buildMCResPxKp1D,
    'mcreskppy1d': buildMCResPyKp1D,
    'mcreskppz1d': buildMCResPzKp1D,

    'mcresmntmkm': buildMCResMntmKm,
    'mcreskmpx': buildMCResPxKm,
    'mcreskmpy': buildMCResPyKm,
    'mcreskmpz': buildMCResPzKm,
    'mcresthetakm': buildMCResThetaKm,
    'mcresphikm': buildMCResPhiKm,
    'mcreskmpx1d': buildMCResPxKm1D,
    'mcreskmpy1d': buildMCResPyKm1D,
    'mcreskmpz1d': buildMCResPzKm1D,

    'q2':buildQ2,
    't':buildt,
    'x':buildx,
    'q2x':buildQ2x,
    'tx':buildtx,
    'q2t':buildQ2t,
    'deltaq2': buildDeltaQ2,
    'deltaxb': buildDeltaXb,
    'deltat': buildDeltaT,
    'deltacosthet': buildDeltaCosTheta,
    'betap': buildBetaP,

    'xeltheta': buildxElTheta,
    'xprtheta': buildxPrTheta,
    'xkptheta': buildxKpTheta,   
    'xkmtheta': buildxKmTheta,
    'xbvsmm2': buildXbMM2

}

#def plotHistos(el, pro, kp, km, mc_el, mc_pr, mc_kp, mc_km ):

    
for ev in ff.clas12:

    el = TLorentzVector()
    el.SetPxPyPzE(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector()
    pro.SetPxPyPzE(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector()
    kp.SetPxPyPzE(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector()
    km.SetPxPyPzE(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)

    mc_el=TLorentzVector()
    mc_el.SetPxPyPzE(ev.mc_el_px,ev.mc_el_py,ev.mc_el_pz,ev.mc_el_e)
    mc_pr=TLorentzVector()
    mc_pr.SetPxPyPzE(ev.mc_pr_px,ev.mc_pr_py,ev.mc_pr_pz,ev.mc_pr_e)
    mc_kp=TLorentzVector()
    mc_kp.SetPxPyPzE(ev.mc_kp_px,ev.mc_kp_py,ev.mc_kp_pz,ev.mc_kp_e)
    mc_km=TLorentzVector()
    mc_km.SetPxPyPzE(ev.mc_km_px,ev.mc_km_py,ev.mc_km_pz,ev.mc_km_e)
    
    mc_res_el_p =el.P() -  mc_el.P() 
    mc_res_el_px =el.Px() -  mc_el.Px() 
    mc_res_el_py =el.Py() -  mc_el.Py() 
    mc_res_el_pz =el.Pz() -  mc_el.Pz() 

    mc_res_el_theta = np.degrees(el.Theta() - mc_el.Theta())
    mc_res_el_phi = np.degrees(el.Phi() - mc_el.Phi())

    mc_res_pr_p =  pro.P() - mc_pr.P()
    mc_res_pr_px = pro.Px() -  mc_pr.Px() 
    mc_res_pr_py = pro.Py() -  mc_pr.Py() 
    mc_res_pr_pz = pro.Pz() -  mc_pr.Pz() 
    mc_res_pr_theta = np.degrees(pro.Theta() - mc_pr.Theta())
    mc_res_pr_phi = np.degrees(pro.Phi() - mc_pr.Phi())

    mc_res_kp_p = kp.P() - mc_kp.P()
    mc_res_kp_px = kp.Px() -  mc_kp.Px() 
    mc_res_kp_py = kp.Py() -  mc_kp.Py() 
    mc_res_kp_pz = kp.Pz() -  mc_kp.Pz() 
    mc_res_kp_theta = np.degrees( kp.Theta() - mc_kp.Theta())
    mc_res_kp_phi = np.degrees(kp.Phi() - mc_kp.Phi())

    mc_res_km_p = km.P() - mc_km.P()
    mc_res_km_px = km.Px() -  mc_km.Px() 
    mc_res_km_py = km.Py() -  mc_km.Py() 
    mc_res_km_pz = km.Pz() -  mc_km.Pz() 
    mc_res_km_theta = np.degrees(km.Theta() - mc_km.Theta())
    mc_res_km_phi = np.degrees(km.Phi() - mc_km.Phi())
    
    pr_status = ev.pr_status
    kp_status = ev.kp_status
    km_status = ev.km_status
    
    prb = ev.prb
    kpb = ev.kpb
    kmb = ev.kmb

    all_final_state_status_FD=False
        
    q = beam-el  
    q2 = -q.M2()  
    xb = q2/(2*target.M()*q.E())
    t = 2*target.M()*(pro.E()-target.M()) 

    mc_q = beam-mc_el  
    mc_q2 = -mc_q.M2()  
    mc_xb = mc_q2/(2*target.M()*mc_q.E())
    mc_t = 2*target.M()*(mc_pr.E()-target.M()) 

    delta_q2 = q2 - mc_q2
    delta_t = t - mc_t
    delta_xb = xb - mc_xb


    eX = beam + target - el
    epX = beam + target - el - pro
    epkpkmX = beam + target - el - pro - kp -km
    epkpX = beam + target - el - pro - kp
    epkmX = beam + target - el - pro - km
    ekpkmX = beam + target - el - kp - km

    cplpro = np.degrees(ekpkmX.Vect().Angle(pro.Vect()))
    cplkm = np.degrees(epkpX.Vect().Angle(km.Vect()))
    cplkp = np.degrees(epkmX.Vect().Angle(kp.Vect()))

    pt = math.sqrt(epkpkmX.Px()**2 + epkpkmX.Py()**2)

    lmda = pro+km
    vphi = kp+km 

    select_phi = vphi.M() < max_phi_mass_signal_limit and vphi.M() > min_phi_mass_signal_limit 

    ##########################################################
    ## monte carlo generated distributions of candidate events
    ## 
    compute_if_absent(histos,'gen_prior_cuts_el_p',builders['p']).Fill(mc_el.P())
    compute_if_absent(histos,'gen_prior_cuts_pr_p',builders['p']).Fill(mc_pr.P())
    compute_if_absent(histos,'gen_prior_cuts_kp_p',builders['p']).Fill(mc_kp.P())
    compute_if_absent(histos,'gen_prior_cuts_km_p',builders['p']).Fill(mc_km.P())
    
    compute_if_absent(histos,'gen_prior_cuts_el_theta',builders['theta']).Fill(np.degrees(mc_el.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_pr_theta',builders['theta']).Fill(np.degrees(mc_pr.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_kp_theta',builders['theta']).Fill(np.degrees(mc_kp.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_km_theta',builders['theta']).Fill(np.degrees(mc_km.Theta()))
    
    compute_if_absent(histos,'gen_prior_cuts_el_phi',builders['phi']).Fill(np.degrees(mc_el.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_pr_phi',builders['phi']).Fill(np.degrees(mc_pr.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_kp_phi',builders['phi']).Fill(np.degrees(mc_kp.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_km_phi',builders['phi']).Fill(np.degrees(mc_km.Phi()))
    
    compute_if_absent(histos,'gen_prior_cuts_el_ptheta',builders['ptheta']).Fill(mc_el.P(),np.degrees(mc_el.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_pr_ptheta',builders['ptheta']).Fill(mc_pr.P(),np.degrees(mc_pr.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_kp_ptheta',builders['ptheta']).Fill(mc_kp.P(),np.degrees(mc_kp.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_km_ptheta',builders['ptheta']).Fill(mc_km.P(),np.degrees(mc_km.Theta()))
    
    compute_if_absent(histos,'gen_prior_cuts_el_thetaphi',builders['phitheta']).Fill(np.degrees(mc_el.Phi()),np.degrees(mc_el.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_pr_thetaphi',builders['phitheta']).Fill(np.degrees(mc_pr.Phi()),np.degrees(mc_pr.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_kp_thetaphi',builders['phitheta']).Fill(np.degrees(mc_kp.Phi()),np.degrees(mc_kp.Theta()))
    compute_if_absent(histos,'gen_prior_cuts_km_thetaphi',builders['phitheta']).Fill(np.degrees(mc_km.Phi()),np.degrees(mc_km.Theta()))
    
    compute_if_absent(histos,'gen_prior_cuts_el_pphi',builders['pphi']).Fill(mc_el.P(),np.degrees(mc_el.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_pr_pphi',builders['pphi']).Fill(mc_pr.P(),np.degrees(mc_pr.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_kp_pphi',builders['pphi']).Fill(mc_kp.P(),np.degrees(mc_kp.Phi()))
    compute_if_absent(histos,'gen_prior_cuts_km_pphi',builders['pphi']).Fill(mc_km.P(),np.degrees(mc_km.Phi()))

    compute_if_absent(histos,'q2_prior_cuts',builders['q2']).Fill(q2)
    compute_if_absent(histos,'t_prior_cuts',builders['t']).Fill(t)
    compute_if_absent(histos,'xb_prior_cuts',builders['x']).Fill(xb)
    
    compute_if_absent(histos,'xq2_prior_cuts',builders['q2x']).Fill(xb,q2)
    compute_if_absent(histos,'txb_prior_cuts',builders['tx']).Fill(t,xb)
    compute_if_absent(histos,'q2t_prior_cuts',builders['q2t']).Fill(q2,t)

    compute_if_absent(histos,'xb_el_theta_prior_cuts',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
    compute_if_absent(histos,'xb_pr_theta_prior_cuts',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
    compute_if_absent(histos,'xb_kp_theta_prior_cuts',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
    compute_if_absent(histos,'xb_km_theta_prior_cuts',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))

    compute_if_absent(histos,'xb_epkpkmXe_prior_cuts',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
    compute_if_absent(histos,'xb_ekpkmXe_prior_cuts',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
    compute_if_absent(histos,'xb_epkpXe_prior_cuts',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
    compute_if_absent(histos,'xb_epkmXe_prior_cuts',builders['xbvsmm2']).Fill(epkmX.M2(), xb)
    

    if( (pr_status >= 2000 and pr_status < 4000 and kp_status >= 2000 and kp_status<4000 and km_status >=2000 and km_status < 4000) ): 
        all_final_state_status_FD=True

        compute_if_absent(histos,'pro_beta_vs_p_raw',builders['betap']).Fill(pro.P(),prb)
        compute_if_absent(histos,'kp_beta_vs_p_raw',builders['betap']).Fill(kp.P(),kpb)
        compute_if_absent(histos,'km_beta_vs_p_raw',builders['betap']).Fill(km.P(),kmb)

        compute_if_absent(histos,'q2_prior_cuts_FD',builders['q2']).Fill(q2)
        compute_if_absent(histos,'t_prior_cuts_FD',builders['t']).Fill(t)
        compute_if_absent(histos,'xb_prior_cuts_FD',builders['x']).Fill(xb)
        
        compute_if_absent(histos,'xq2_prior_cuts_FD',builders['q2x']).Fill(xb,q2)
        compute_if_absent(histos,'txb_prior_cuts_FD',builders['tx']).Fill(t,xb)
        compute_if_absent(histos,'q2t_prior_cuts_FD',builders['q2t']).Fill(q2,t)

        compute_if_absent(histos,'xb_el_theta_prior_cuts_FD',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
        compute_if_absent(histos,'xb_pr_theta_prior_cuts_FD',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
        compute_if_absent(histos,'xb_kp_theta_prior_cuts_FD',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
        compute_if_absent(histos,'xb_km_theta_prior_cuts_FD',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))

        if passEPKPKMX(epkpkmX) :
            compute_if_absent(histos,'xb_el_theta_epkpkmXcut_FD',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
            compute_if_absent(histos,'xb_pr_theta_epkpkmXcut_FD',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
            compute_if_absent(histos,'xb_kp_theta_epkpkmXcut_FD',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
            compute_if_absent(histos,'xb_km_theta_epkpkmXcut_FD',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))
            compute_if_absent(histos,'xb_epkpkmXe_epkpkmXcut_FD',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
            compute_if_absent(histos,'xb_ekpkmX_epkpkmXcut_FD',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
            compute_if_absent(histos,'xb_epkpX_epkpkmXcut_FD',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
            compute_if_absent(histos,'xb_epkmX_epkpkmXcut_FD',builders['xbvsmm2']).Fill(epkmX.M2(), xb)


        if passEPKPX(epkpX) :
            compute_if_absent(histos,'xb_el_theta_epkpXcut_FD',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
            compute_if_absent(histos,'xb_pr_theta_epkpXcut_FD',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
            compute_if_absent(histos,'xb_kp_theta_epkpXcut_FD',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
            compute_if_absent(histos,'xb_km_theta_epkpXcut_FD',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))
            compute_if_absent(histos,'xb_epkpkmXe_epkpXcut_FD',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
            compute_if_absent(histos,'xb_ekpkmX_epkpXcut_FD',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
            compute_if_absent(histos,'xb_epkpX_epkpXcut_FD',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
            compute_if_absent(histos,'xb_epkmX_epkpXcut_FD',builders['xbvsmm2']).Fill(epkmX.M2(), xb)

        if passEPKMX(epkmX) :
            compute_if_absent(histos,'xb_el_theta_epkmXcut_FD',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
            compute_if_absent(histos,'xb_pr_theta_epkmXcut_FD',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
            compute_if_absent(histos,'xb_kp_theta_epkmXcut_FD',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
            compute_if_absent(histos,'xb_km_theta_epkmXcut_FD',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))
            compute_if_absent(histos,'xb_epkpkmXe_epkmXcut_FD',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
            compute_if_absent(histos,'xb_ekpkmX_epkmXcut_FD',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
            compute_if_absent(histos,'xb_epkpX_epkmXcut_FD',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
            compute_if_absent(histos,'xb_epkmX_epkmXcut_FD',builders['xbvsmm2']).Fill(epkmX.M2(), xb)

        if passEKPKMX(ekpkmX) :
            print('here')
            compute_if_absent(histos,'xb_el_theta_ekpkmXcut_FD',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
            compute_if_absent(histos,'xb_pr_theta_ekpkmXcut_FD',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
            compute_if_absent(histos,'xb_kp_theta_ekpkmXcut_FD',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
            compute_if_absent(histos,'xb_km_theta_ekpkmXcut_FD',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))
            compute_if_absent(histos,'xb_epkpkmXe_epkpkmXcut_FD',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
            compute_if_absent(histos,'xb_ekpkmX_ekpkmXcut_FD',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
            compute_if_absent(histos,'xb_epkpX_ekpkmXcut_FD',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
            compute_if_absent(histos,'xb_epkmX_ekpkmXcut_FD',builders['xbvsmm2']).Fill(epkmX.M2(), xb)



        if( all( passEvent(epkpkmX, epkpX, epkmX, ekpkmX) ) ): # and select_phi ):

        ################
        ## kin plots
            compute_if_absent(histos,'el_p',builders['p']).Fill(el.P())
            compute_if_absent(histos,'pr_p',builders['p']).Fill(pro.P())
            compute_if_absent(histos,'kp_p',builders['p']).Fill(kp.P())
            compute_if_absent(histos,'km_p',builders['p']).Fill(km.P())
            
            compute_if_absent(histos,'el_theta',builders['theta']).Fill(np.degrees(el.Theta()))
            compute_if_absent(histos,'pr_theta',builders['theta']).Fill(np.degrees(pro.Theta()))
            compute_if_absent(histos,'kp_theta',builders['theta']).Fill(np.degrees(kp.Theta()))
            compute_if_absent(histos,'km_theta',builders['theta']).Fill(np.degrees(km.Theta()))
            
            compute_if_absent(histos,'el_phi',builders['phi']).Fill(np.degrees(el.Phi()))
            compute_if_absent(histos,'pr_phi',builders['phi']).Fill(np.degrees(pro.Phi()))
            compute_if_absent(histos,'kp_phi',builders['phi']).Fill(np.degrees(kp.Phi()))
            compute_if_absent(histos,'km_phi',builders['phi']).Fill(np.degrees(km.Phi()))
            
            compute_if_absent(histos,'el_ptheta',builders['ptheta']).Fill(el.P(),np.degrees(el.Theta()))
            compute_if_absent(histos,'pr_ptheta',builders['ptheta']).Fill(pro.P(),np.degrees(pro.Theta()))
            compute_if_absent(histos,'kp_ptheta',builders['ptheta']).Fill(kp.P(),np.degrees(kp.Theta()))
            compute_if_absent(histos,'km_ptheta',builders['ptheta']).Fill(km.P(),np.degrees(km.Theta()))
            
            compute_if_absent(histos,'el_thetaphi',builders['phitheta']).Fill(np.degrees(el.Phi()),np.degrees(el.Theta()))
            compute_if_absent(histos,'pr_thetaphi',builders['phitheta']).Fill(np.degrees(pro.Phi()),np.degrees(pro.Theta()))
            compute_if_absent(histos,'kp_thetaphi',builders['phitheta']).Fill(np.degrees(kp.Phi()),np.degrees(kp.Theta()))
            compute_if_absent(histos,'km_thetaphi',builders['phitheta']).Fill(np.degrees(km.Phi()),np.degrees(km.Theta()))

            compute_if_absent(histos,'el_pphi',builders['pphi']).Fill(el.P(),np.degrees(el.Phi()))
            compute_if_absent(histos,'pr_pphi',builders['pphi']).Fill(pro.P(),np.degrees(pro.Phi()))
            compute_if_absent(histos,'kp_pphi',builders['pphi']).Fill(kp.P(),np.degrees(kp.Phi()))
            compute_if_absent(histos,'km_pphi',builders['pphi']).Fill(km.P(),np.degrees(km.Phi()))

            compute_if_absent(histos,'q2',builders['q2']).Fill(q2)
            compute_if_absent(histos,'t',builders['t']).Fill(t)
            compute_if_absent(histos,'xb',builders['x']).Fill(xb)
            
            compute_if_absent(histos,'xq2',builders['q2x']).Fill(xb,q2)
            compute_if_absent(histos,'txb',builders['tx']).Fill(t,xb)
            compute_if_absent(histos,'q2t',builders['q2t']).Fill(q2,t)
            
            compute_if_absent(histos,'xb_el_theta',builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
            compute_if_absent(histos,'xb_pr_theta',builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
            compute_if_absent(histos,'xb_kp_theta',builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
            compute_if_absent(histos,'xb_km_theta',builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))

            compute_if_absent(histos,'xb_epkpkmXe',builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
            compute_if_absent(histos,'xb_ekpkmX',builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
            compute_if_absent(histos,'xb_epkpX',builders['xbvsmm2']).Fill(epkpX.M2(), xb)
            compute_if_absent(histos,'xb_epkmX',builders['xbvsmm2']).Fill(epkmX.M2(), xb)


            compute_if_absent(histos,'eX',builders['mm2']).Fill( eX.M2() )
            compute_if_absent(histos,'epX',builders['mm2']).Fill( epX.M2() )
            compute_if_absent(histos,'epkpX',builders['mm2']).Fill( epkpX.M2() )
            compute_if_absent(histos,'epkmX',builders['mm2']).Fill( epkmX.M2() )
            compute_if_absent(histos,'ekpkmX',builders['mm2']).Fill( ekpkmX.M2() )
            compute_if_absent(histos,'epkpkmX',builders['mm2']).Fill( epkpkmX.M2() )
            compute_if_absent(histos,'epkpkmXe',builders['mm2']).Fill( epkpkmX.E() )
            compute_if_absent(histos,'vphimass',builders['mm2']).Fill( vphi.M() )
                        
            compute_if_absent(histos,'LeX',builders['mm2L']).Fill( eX.M2() )
            compute_if_absent(histos,'LepX',builders['mm2L']).Fill( epX.M2() )
            compute_if_absent(histos,'LepkpX',builders['mm2L']).Fill( epkpX.M2() )
            compute_if_absent(histos,'LepkmX',builders['mm2L']).Fill( epkmX.M2() )
            compute_if_absent(histos,'LekpkmX',builders['mm2L']).Fill( ekpkmX.M2() )
            compute_if_absent(histos,'LepkpkmX',builders['mm2L']).Fill( epkpkmX.M2() )
            compute_if_absent(histos,'LepkpkmXe',builders['mm2L']).Fill( epkpkmX.E() )
            compute_if_absent(histos,'Lvphimass',builders['mm2L']).Fill( vphi.M() )
            compute_if_absent(histos,'pt',builders['mm2']).Fill(pt)
                        
            ################
            ## epX 
            compute_if_absent(histos,'epX_vs_theta_pr',builders['mmthetaEPX']).Fill(np.degrees(pro.Theta()), epX.M2())
            compute_if_absent(histos,'epX_vs_p_pr',builders['mmpEPX']).Fill(pro.P(), epX.M2())
            compute_if_absent(histos,'kpkm_vs_theta_pr',builders['mmthetaEPX']).Fill(np.degrees(pro.Theta()), vphi.M())
            compute_if_absent(histos,'kpkm_vs_p_pr',builders['mmpEPX']).Fill(pro.P(), vphi.M())
                        
            ##epkpX
            compute_if_absent(histos,'epkpX_vs_theta_km',builders['mmthetaEPKPX']).Fill(np.degrees(km.Theta()), epkpX.M2())
            compute_if_absent(histos,'epkpX_vs_p_km',builders['mmpEPKPX']).Fill(km.P(), epkpX.M2())
            
            ##epkmX
            compute_if_absent(histos,'epkmX_vs_theta_kp',builders['mmthetaEPKMX']).Fill(np.degrees(kp.Theta()), epkmX.M2())
            compute_if_absent(histos,'epkmX_vs_p_kp',builders['mmpEPKMX']).Fill(kp.P(), epkmX.M2())
            
            ##ekpkmX
            compute_if_absent(histos,'ekpkmX_vs_theta_pr',builders['mmthetaEKPKMX']).Fill(np.degrees(pro.Theta()), ekpkmX.M2())
            compute_if_absent(histos,'ekpkmX_vs_p_pr',builders['mmpEKPKMX']).Fill(pro.P(), ekpkmX.M2())
            
            ##ep->epPhi
            ##delta p and delta theta
            delta_p_pr = pro.P() - ekpkmX.P() 
            delta_theta_pr = np.degrees(pro.Theta() - ekpkmX.Theta())
            
            delta_p_km = km.P() - epkpX.P()
            delta_theta_km = np.degrees(km.Theta() - epkpX.Theta())
            
            delta_p_kp = kp.P() - epkmX.P() 
            delta_theta_kp = np.degrees( kp.Theta() - epkmX.Theta())

            compute_if_absent(histos,'deltaQ2',builders['deltaq2']).Fill(mc_q2, delta_q2)
            compute_if_absent(histos,'deltaT',builders['deltat']).Fill(mc_t, delta_t)
            compute_if_absent(histos,'deltaXb',builders['deltaxb']).Fill(mc_xb, delta_xb)

            
            compute_if_absent(histos,'delta_theta_vs_theta_pr',builders['delthetPrBig']).Fill(np.degrees(pro.Theta()), delta_theta_pr)
            compute_if_absent(histos,'delta_theta_vs_p_pr',builders['delthetPrSmall']).Fill(pro.P(), delta_theta_pr)
            
            compute_if_absent(histos,'delta_p_vs_theta_pr',builders['delthetPrMed']).Fill(np.degrees(pro.Theta()), delta_p_pr)
            compute_if_absent(histos,'delta_p_vs_p_pr',builders['delpPrSmall']).Fill(pro.P(), delta_p_pr)

            ##################################
            ## km
            compute_if_absent(histos,'delta_theta_vs_theta_km',builders['delthetKmBig']).Fill(np.degrees(km.Theta()), delta_theta_km)
            compute_if_absent(histos,'delta_theta_vs_p_km',builders['delthetKmSmall']).Fill(km.P(), delta_theta_km)
            
            compute_if_absent(histos,'delta_p_vs_theta_km',builders['delthetKmMed']).Fill(np.degrees(km.Theta()), delta_p_km)
            compute_if_absent(histos,'delta_p_vs_p_km',builders['delpKmSmall']).Fill(km.P(), delta_p_km)
            
            ##################################
            ## kp
            compute_if_absent(histos,'delta_theta_vs_theta_kp',builders['delthetKpBig']).Fill(np.degrees(kp.Theta()), delta_theta_kp)
            compute_if_absent(histos,'delta_theta_vs_p_kp',builders['delthetKpSmall']).Fill(kp.P(), delta_theta_kp)
            
            compute_if_absent(histos,'delta_p_vs_theta_kp',builders['delthetKpMed']).Fill(np.degrees(kp.Theta()), delta_p_kp)
            compute_if_absent(histos,'delta_p_vs_p_kp',builders['delpKpSmall']).Fill(kp.P(), delta_p_kp)
            
            ##################################
            ## monte carlo simulation resolution in forward detector
            ## generated - reconstructed
            compute_if_absent(histos,'gen_el_p',builders['p']).Fill(mc_el.P())
            compute_if_absent(histos,'gen_pr_p',builders['p']).Fill(mc_pr.P())
            compute_if_absent(histos,'gen_kp_p',builders['p']).Fill(mc_kp.P())
            compute_if_absent(histos,'gen_km_p',builders['p']).Fill(mc_km.P())

            compute_if_absent(histos,'gen_el_theta',builders['theta']).Fill(np.degrees(mc_el.Theta()))
            compute_if_absent(histos,'gen_pr_theta',builders['theta']).Fill(np.degrees(mc_pr.Theta()))
            compute_if_absent(histos,'gen_kp_theta',builders['theta']).Fill(np.degrees(mc_kp.Theta()))
            compute_if_absent(histos,'gen_km_theta',builders['theta']).Fill(np.degrees(mc_km.Theta()))
            
            compute_if_absent(histos,'gen_el_phi',builders['phi']).Fill(np.degrees(mc_el.Phi()))
            compute_if_absent(histos,'gen_pr_phi',builders['phi']).Fill(np.degrees(mc_pr.Phi()))
            compute_if_absent(histos,'gen_kp_phi',builders['phi']).Fill(np.degrees(mc_kp.Phi()))
            compute_if_absent(histos,'gen_km_phi',builders['phi']).Fill(np.degrees(mc_km.Phi()))
            
            compute_if_absent(histos,'gen_el_ptheta',builders['ptheta']).Fill(mc_el.P(),np.degrees(mc_el.Theta()))
            compute_if_absent(histos,'gen_pr_ptheta',builders['ptheta']).Fill(mc_pr.P(),np.degrees(mc_pr.Theta()))
            compute_if_absent(histos,'gen_kp_ptheta',builders['ptheta']).Fill(mc_kp.P(),np.degrees(mc_kp.Theta()))
            compute_if_absent(histos,'gen_km_ptheta',builders['ptheta']).Fill(mc_km.P(),np.degrees(mc_km.Theta()))

            compute_if_absent(histos,'gen_el_thetaphi',builders['phitheta']).Fill(np.degrees(mc_el.Phi()),np.degrees(mc_el.Theta()))
            compute_if_absent(histos,'gen_pr_thetaphi',builders['phitheta']).Fill(np.degrees(mc_pr.Phi()),np.degrees(mc_pr.Theta()))
            compute_if_absent(histos,'gen_kp_thetaphi',builders['phitheta']).Fill(np.degrees(mc_kp.Phi()),np.degrees(mc_kp.Theta()))
            compute_if_absent(histos,'gen_km_thetaphi',builders['phitheta']).Fill(np.degrees(mc_km.Phi()),np.degrees(mc_km.Theta()))
            
            compute_if_absent(histos,'gen_el_pphi',builders['pphi']).Fill(mc_el.P(),np.degrees(mc_el.Phi()))
            compute_if_absent(histos,'gen_pr_pphi',builders['pphi']).Fill(mc_el.P(),np.degrees(mc_pr.Phi()))
            compute_if_absent(histos,'gen_kp_pphi',builders['pphi']).Fill(mc_el.P(),np.degrees(mc_kp.Phi()))
            compute_if_absent(histos,'gen_km_pphi',builders['pphi']).Fill(mc_el.P(),np.degrees(mc_km.Phi()))
                        
            mc_eX = beam + target - mc_el
            mc_epX = beam + target - mc_el - mc_pr
            mc_epkpkmX = beam + target - mc_el - mc_pr - mc_kp - mc_km
            mc_epkpX = beam + target - mc_el - mc_pr - mc_kp
            mc_epkmX = beam + target - mc_el - mc_pr - mc_km
            mc_ekpkmX = beam + target - mc_el - mc_kp - mc_km
            mc_pt = math.sqrt(mc_epkpkmX.Px()**2 + mc_epkpkmX.Py()**2)
            
            compute_if_absent(histos,'gen_eX',builders['mm2']).Fill( mc_eX.M2() )
            compute_if_absent(histos,'gen_epX',builders['mm2']).Fill( mc_epX.M2() )
            compute_if_absent(histos,'gen_epkpX',builders['mm2']).Fill( mc_epkpX.M2() )
            compute_if_absent(histos,'gen_epkmX',builders['mm2']).Fill( mc_epkmX.M2() )
            compute_if_absent(histos,'gen_ekpkmX',builders['mm2']).Fill( mc_ekpkmX.M2() )
            compute_if_absent(histos,'gen_epkpkmX',builders['mm2']).Fill( mc_epkpkmX.M2() )
            compute_if_absent(histos,'gen_epkpkmXe',builders['mm2']).Fill( mc_epkpkmX.E() )
                        
            compute_if_absent(histos,'gen_LeX',builders['mm2L']).Fill( mc_eX.M2() )
            compute_if_absent(histos,'gen_LepX',builders['mm2L']).Fill( mc_epX.M2() )
            compute_if_absent(histos,'gen_LepkpX',builders['mm2L']).Fill( mc_epkpX.M2() )
            compute_if_absent(histos,'gen_LepkmX',builders['mm2L']).Fill( mc_epkmX.M2() )
            compute_if_absent(histos,'gen_LekpkmX',builders['mm2L']).Fill( mc_ekpkmX.M2() )
            compute_if_absent(histos,'gen_LepkpkmX',builders['mm2L']).Fill( mc_epkpkmX.M2() )
            compute_if_absent(histos,'gen_LepkpkmXe',builders['mm2L']).Fill( mc_epkpkmX.E() )
            compute_if_absent(histos,'gen_pt',builders['mm2']).Fill(mc_pt)
            
            compute_if_absent(histos,'mc_res_el_p_FD',builders['mcresmntmel']).Fill(mc_el.P(), mc_res_el_p)
            compute_if_absent(histos,'mc_res_el_px_1d',builders['mcreselpx1d']).Fill(mc_res_el_px)
            compute_if_absent(histos,'mc_res_el_py_1d',builders['mcreselpy1d']).Fill(mc_res_el_py)
            compute_if_absent(histos,'mc_res_el_pz_1d',builders['mcreselpz1d']).Fill(mc_res_el_px)        
            compute_if_absent(histos,'mc_res_el_px_FD',builders['mcreselpx']).Fill(mc_el.Px(), mc_res_el_px)
            compute_if_absent(histos,'mc_res_el_py_FD',builders['mcreselpy']).Fill(mc_el.Py(), mc_res_el_py)
            compute_if_absent(histos,'mc_res_el_pz_FD',builders['mcreselpz']).Fill(mc_el.Pz(), mc_res_el_pz)            
            compute_if_absent(histos,'mc_res_el_theta_FD',builders['mcresthetael']).Fill(np.degrees(mc_el.Theta()), mc_res_el_theta)
            compute_if_absent(histos,'mc_res_el_phi_FD',builders['mcresphiel']).Fill(np.degrees(mc_el.Phi()), mc_res_el_phi)
            
            compute_if_absent(histos,'mc_res_pr_p_FD',builders['mcresmntmpr']).Fill(mc_pr.P(), mc_res_pr_p)
            compute_if_absent(histos,'mc_res_pr_px_1d',builders['mcresprpx1d']).Fill(mc_res_pr_px)
            compute_if_absent(histos,'mc_res_pr_py_1d',builders['mcresprpy1d']).Fill(mc_res_pr_py)
            compute_if_absent(histos,'mc_res_pr_pz_1d',builders['mcresprpz1d']).Fill(mc_res_pr_px)        
            compute_if_absent(histos,'mc_res_pr_px_FD',builders['mcresprpx']).Fill(mc_pr.Px(), mc_res_pr_px)
            compute_if_absent(histos,'mc_res_pr_py_FD',builders['mcresprpy']).Fill(mc_pr.Py(), mc_res_pr_py)
            compute_if_absent(histos,'mc_res_pr_pz_FD',builders['mcresprpz']).Fill(mc_pr.Pz(), mc_res_pr_pz)
            compute_if_absent(histos,'mc_res_pr_theta_FD',builders['mcresthetapr']).Fill(np.degrees(mc_pr.Theta()), mc_res_pr_theta)
            compute_if_absent(histos,'mc_res_pr_phi_FD',builders['mcresphipr']).Fill(np.degrees(mc_pr.Phi()), mc_res_pr_phi)
            low_phi = -180
            for ii in range(1,7):
                max_phi = low_phi + ii*60.0
                min_phi = low_phi + (ii-1)*60.0
                if( np.degrees(pro.Phi()) > min_phi and np.degrees(pro.Phi()) < max_phi ):
                    compute_if_absent(histos,'mc_res_pr_p_FD_minphi'+str(min_phi)+'_maxphi'+str(max_phi),builders['mcresmntmpr']).Fill(mc_pr.P(), mc_res_pr_p)

            low_theta=10
            for ii in range(1,8):
                max_theta = low_theta + ii*2.5
                min_theta = low_theta + (ii-1)*2.5
                if( np.degrees(pro.Theta()) > min_theta and np.degrees(pro.Theta()) < max_theta):
                    compute_if_absent(histos,'mc_res_pr_p_FD_mintheta'+str(min_theta)+'_maxtheta'+str(max_theta),builders['mcresmntmpr']).Fill(mc_pr.P(), mc_res_pr_p)
                

            compute_if_absent(histos,'mc_res_kp_p_FD',builders['mcresmntmkp']).Fill(mc_kp.P(), mc_res_kp_p)
            compute_if_absent(histos,'mc_res_kp_px_1d',builders['mcreskppx1d']).Fill(mc_res_kp_px)
            compute_if_absent(histos,'mc_res_kp_py_1d',builders['mcreskppy1d']).Fill(mc_res_kp_py)
            compute_if_absent(histos,'mc_res_kp_pz_1d',builders['mcreskppz1d']).Fill(mc_res_kp_px)        
            compute_if_absent(histos,'mc_res_kp_px_FD',builders['mcreskppx']).Fill(mc_kp.Px(), mc_res_kp_px)
            compute_if_absent(histos,'mc_res_kp_py_FD',builders['mcreskppy']).Fill(mc_kp.Py(), mc_res_kp_py)
            compute_if_absent(histos,'mc_res_kp_pz_FD',builders['mcreskppz']).Fill(mc_kp.Pz(), mc_res_kp_pz)
            compute_if_absent(histos,'mc_res_kp_theta_FD',builders['mcresthetakp']).Fill(np.degrees(mc_kp.Theta()), mc_res_kp_theta)
            compute_if_absent(histos,'mc_res_kp_phi_FD',builders['mcresphikp']).Fill(np.degrees(mc_kp.Phi()), mc_res_kp_phi)
            
            compute_if_absent(histos,'mc_res_km_p_FD',builders['mcresmntmkm']).Fill(mc_km.P(), mc_res_km_p)
            compute_if_absent(histos,'mc_res_km_px_1d',builders['mcreskmpx1d']).Fill(mc_res_km_px)
            compute_if_absent(histos,'mc_res_km_py_1d',builders['mcreskmpy1d']).Fill(mc_res_km_py)
            compute_if_absent(histos,'mc_res_km_pz_1d',builders['mcreskmpz1d']).Fill(mc_res_km_px)        
            compute_if_absent(histos,'mc_res_km_px_FD',builders['mcreskmpx']).Fill(mc_km.Px(), mc_res_km_px)
            compute_if_absent(histos,'mc_res_km_py_FD',builders['mcreskmpy']).Fill(mc_km.Py(), mc_res_km_py)
            compute_if_absent(histos,'mc_res_km_pz_FD',builders['mcreskmpz']).Fill(mc_km.Pz(), mc_res_km_pz)
            compute_if_absent(histos,'mc_res_km_theta_FD',builders['mcresthetakm']).Fill(np.degrees(mc_km.Theta()), mc_res_km_theta)
            compute_if_absent(histos,'mc_res_km_phi_FD',builders['mcresphikm']).Fill(np.degrees(mc_km.Phi()), mc_res_km_phi)
            
            compute_if_absent(histos,'rec_pr_ptheta_FD',builders['ptheta']).Fill((pro.P()), np.degrees(pro.Theta()))
            compute_if_absent(histos,'rec_pr_ptheta_FD_wres',builders['ptheta']).Fill((pro.P()), np.degrees(pro.Theta()), mc_res_pr_p) 
            
                    
        
fout = TFile('mc_res_phi_kin_ebc_heb_'+tag+'.root','recreate')
for _, hist in histos.items():
    hist.Write()

fout.Close()



    
    
