from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from collections import OrderedDict

import numpy as np
import os

from array import array
import math

#include my Exclusive Cuts class to select the events
import sys
sys.path.append('/work/clas12/bclary/CLAS12/analysis_code_fork0/projects/compareDataSim/src')
from ExclusiveCuts import ExclusiveCuts

ff = TFile(sys.argv[1])
datatype = sys.argv[2]
tag=sys.argv[3]

base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution

epkpkmxe_min=0
epkpkmxe_max=0
epkpX_min=0
epkpX_max=0
epkmX_min=0
epkmX_max=0
ekpkmX_min=0
ekpkmX_max=0
sig=5

excl_cuts = ExclusiveCuts()
excl_cuts.setCuts(datatype,tag,sig)

if datatype == 'sim' and 'inb' in tag or 'outb' in tag :
    print('init simulation cuts')
    epkpkmxe_min = 0.0222 - 0.0262*sig
    epkpkmxe_max = 0.0222 + 0.0262*sig
    
    epkpX_min = 0.234 - 0.0348*sig
    epkpX_max = 0.234 + 0.0348*sig
    
    epkmX_min = 0.228 - 0.0284*sig
    epkmX_max = 0.228 + 0.0284*sig
    
    ekpkmX_min = 0.891 - 0.0285*sig
    ekpkmX_max = 0.891 + 0.0285*sig

if datatype == 'data' and 'inb' in tag:
    epkpkmxe_min = 0.08 - 0.0398*sig
    epkpkmxe_max = 0.08 + 0.0398*sig
    
    epkpX_min = 0.248 - 0.059*sig
    epkpX_max = 0.248 + 0.059*sig
    
    epkmX_min = 0.248 - 0.055*sig
    epkmX_max = 0.248 + 0.055*sig
    
    ekpkmX_min = 0.904 - 0.0719*sig
    ekpkmX_max = 0.904 + 0.0719*sig

if datatype == 'data' and 'inb' in tag:
    epkpkmxe_min = 0.08 - 0.0398*sig
    epkpkmxe_max = 0.08 + 0.0398*sig
    
    epkpX_min = 0.248 - 0.059*sig
    epkpX_max = 0.248 + 0.059*sig
    
    epkmX_min = 0.248 - 0.055*sig
    epkmX_max = 0.248 + 0.055*sig
    
    ekpkmX_min = 0.904 - 0.0719*sig
    ekpkmX_max = 0.904 + 0.0719*sig


def passEPKPKMX(epkpkmX):
    return ( (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min ))

def passEPKPX(epkpX):
    return epkpX.M2() < epkpX_max and epkpX.M2() > epkpX_min

def passEPKMX(epkmX):
    return epkmX.M2() < epkmX_max and epkmX.M2() > epkmX_min

def passEKPKMX(ekpkmX):
    return ekpkmX.M2() < ekpkmX_max and ekpkmX.M2() > ekpkmX_min
    
 
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
def buildCpl(title):
    print('building {}'.format(title))
    return TH1F(title,title,200,0.0,60.0)
def buildKin2D(title):
    print('building {}'.format(title))
    return TH2F(title,title,400,0.0,beam.E(), 400, 0.0, 6.0)
def buildCos(title):
    print('building {}'.format(title))
    return TH1F(title,title,20,-1,1)
def buildAsy(title):
    print('building {}'.format(title))
    return TH1F(title,title,40,-180.0,180.0)
def buildHel(title):
    print('building {}'.format(title))
    return TH1F(title,title,4,-2,2)
def buildQ2x(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,1,100,0.0,beam.E())
def buildtx(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,6,100,0.0,1.0)
def buildQ2t(title):
    print('building {}'.format(title))
    return TH2F(title,title,100, 0, beam.E(), 100, 0.0, 6.0)
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
def buildXbMM2(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -0.25, 1.25, 250, 0.0, 1.1)

    
builders={
    'p':buildP,
    'theta':buildTheta,
    'ptheta':buildPTheta,
    'phitheta':buildPhiTheta,
    'mm2': buildMM2,
    'mass': buildMass,
    'mm2L': buildMM2Large,
    'massL': buildMassLarge,
    'mass2D': buildMass2D,
    'mass2DL': buildMass2DLarge,
    'masskin': buildMassKin,
    'cpl': buildCpl,
    'kin2d':buildKin2D,
    'cos':buildCos,
    'asy':buildAsy,
    'hel':buildHel,
    'q2x':buildQ2x,
    'q2t':buildQ2t,
    'tx':buildtx,
    'betap': buildBetaP,
    'xeltheta': buildxElTheta,
    'xprtheta': buildxPrTheta,
    'xkptheta': buildxKpTheta,   
    'xkmtheta': buildxKmTheta,
    'xbvsmm2': buildXbMM2


}


def plotHistos( el, pro, kp, km, hel, cut_name ):
    eX = beam + target - el
    epX = beam + target - el - pro
    ekpX = beam + target - el - kp    
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
    q = beam-el  
    q2 = -q.M2()  
    xb = q2/(2*target.M()*q.E())
    t = 2*target.M()*(pro.E()-target.M()) 
    nu = q2/(2*target.M()*xb)
    y = nu/beam.E() 

    v3l = beam.Vect().Cross(el.Vect())
    v3h = pro.Vect().Cross((beam-el).Vect()) 
    trento = np.degrees(TMath.ACos(v3l.Dot(v3h)/(v3l.Mag()*v3h.Mag())))
    if(v3l.Dot(pro.Vect())<0): trento=-trento

    #v3l = beam.Vect().Cross(el.Vect())
    v3hphi = vphi.Vect().Cross((beam-el).Vect()) 
    trentophi = np.degrees(TMath.ACos(v3l.Dot(v3h)/(v3l.Mag()*v3h.Mag())))
    if(v3l.Dot(vphi.Vect())<0): trentophi=-trentophi
    
    boost_phicm = vphi.BoostVector()
    boost_phicm = -boost_phicm
    boost_kp = TLorentzVector(kp)
    boost_km = TLorentzVector(km)
    boost_kp.Boost(boost_phicm)
    boost_km.Boost(boost_phicm)
    cos_theta_cm_kp = TMath.Cos(boost_kp.Vect().Theta())

    kp_as_pip = TLorentzVector(kp.Px(), kp.Py(), kp.Pz(), np.sqrt(kp.P()**2 + 0.1396**2))
    vphi_kp_as_pip = kp_as_pip + km

    compute_if_absent(histos,'el_p_'+cut_name,builders['p']).Fill( el.P() )
    compute_if_absent(histos,'pr_p_'+cut_name,builders['p']).Fill( pro.P() )
    compute_if_absent(histos,'kp_p_'+cut_name,builders['p']).Fill( kp.P() )
    compute_if_absent(histos,'km_p_'+cut_name,builders['p']).Fill( km.P() )

    compute_if_absent(histos,'el_theta_'+cut_name,builders['theta']).Fill( np.degrees(el.Theta()) )
    compute_if_absent(histos,'pr_theta_'+cut_name,builders['theta']).Fill( np.degrees(pro.Theta()) )
    compute_if_absent(histos,'kp_theta_'+cut_name,builders['theta']).Fill( np.degrees(kp.Theta()) )
    compute_if_absent(histos,'km_theta_'+cut_name,builders['theta']).Fill( np.degrees(km.Theta()) )
    
    compute_if_absent(histos,'el_ptheta_'+cut_name,builders['ptheta']).Fill( el.P(),  np.degrees(el.Theta()) )
    compute_if_absent(histos,'pr_ptheta_'+cut_name,builders['ptheta']).Fill( pro.P(),  np.degrees(pro.Theta()) )
    compute_if_absent(histos,'kp_ptheta_'+cut_name,builders['ptheta']).Fill( kp.P(),  np.degrees(kp.Theta()) )
    compute_if_absent(histos,'km_ptheta_'+cut_name,builders['ptheta']).Fill( km.P(),  np.degrees(km.Theta()) )
    
    compute_if_absent(histos,'el_phitheta_'+cut_name,builders['phitheta']).Fill( np.degrees(el.Phi()),  np.degrees(el.Theta()) )
    compute_if_absent(histos,'pr_phitheta_'+cut_name,builders['phitheta']).Fill(np.degrees(pro.Phi()),  np.degrees(pro.Theta()) )
    compute_if_absent(histos,'kp_phitheta_'+cut_name,builders['phitheta']).Fill( np.degrees(kp.Phi()),  np.degrees(kp.Theta()) )
    compute_if_absent(histos,'km_phitheta_'+cut_name,builders['phitheta']).Fill( np.degrees(km.Phi()),  np.degrees(km.Theta()) )
    
    compute_if_absent(histos,'eX_'+cut_name,builders['mm2']).Fill( eX.M2() )
    compute_if_absent(histos,'epX_'+cut_name,builders['mm2']).Fill( epX.M2() )
    compute_if_absent(histos,'ekpX_'+cut_name,builders['mm2']).Fill( ekpX.M2() )
    compute_if_absent(histos,'epkpX_'+cut_name,builders['mm2']).Fill( epkpX.M2() )
    compute_if_absent(histos,'epkmX_'+cut_name,builders['mm2']).Fill( epkmX.M2() )
    compute_if_absent(histos,'ekpkmX_'+cut_name,builders['mm2']).Fill( ekpkmX.M2() )
    compute_if_absent(histos,'epkpkmX_'+cut_name,builders['mm2']).Fill( epkpkmX.M2() )
    compute_if_absent(histos,'epkpkmXe_'+cut_name,builders['mm2']).Fill( epkpkmX.E() )
    

    compute_if_absent(histos,'LeX_'+cut_name,builders['mm2L']).Fill( eX.M2() )
    compute_if_absent(histos,'LepX_'+cut_name,builders['mm2L']).Fill( epX.M2() )
    compute_if_absent(histos,'LekpX_'+cut_name,builders['mm2L']).Fill( ekpX.M2() )
    compute_if_absent(histos,'LepkpX_'+cut_name,builders['mm2L']).Fill( epkpX.M2() )
    compute_if_absent(histos,'LepkmX_'+cut_name,builders['mm2L']).Fill( epkmX.M2() )
    compute_if_absent(histos,'LekpkmX_'+cut_name,builders['mm2L']).Fill( ekpkmX.M2() )
    compute_if_absent(histos,'LepkpkmX_'+cut_name,builders['mm2L']).Fill( epkpkmX.M2() )
    compute_if_absent(histos,'LepkpkmXe_'+cut_name,builders['mm2L']).Fill( epkpkmX.E() )

    compute_if_absent(histos,'mntm_el_epkpkmX_'+cut_name,builders['masskin']).Fill(el.P(), epkpkmX.M2() )
    compute_if_absent(histos,'mntm_pr_ekpkmX_'+cut_name,builders['masskin']).Fill(pro.P(), ekpkmX.M2() )
    compute_if_absent(histos,'mntm_kp_epkpkmX_'+cut_name,builders['masskin']).Fill(kp.P(), epkmX.M2() )
    compute_if_absent(histos,'mntm_km_epkpkmX_'+cut_name,builders['masskin']).Fill(km.P(), epkpX.M2() )

    compute_if_absent(histos,'kp_as_pip_'+cut_name,builders['masskin']).Fill(kp_as_pip.P(), vphi_kp_as_pip.M())
    
    compute_if_absent(histos,'vphimass_'+cut_name,builders['mass']).Fill(vphi.M())
    compute_if_absent(histos,'lambdamass_'+cut_name,builders['massL']).Fill(lmda.M())
    compute_if_absent(histos,'lambdaphimass_'+cut_name,builders['mass2D']).Fill(lmda.M(), vphi.M())

    compute_if_absent(histos,'kp_mntm_phimass_'+cut_name,builders['mass2DL']).Fill(kp.P(), vphi.M())
    compute_if_absent(histos,'km_mntm_phimass_'+cut_name,builders['mass2DL']).Fill(km.P(), vphi.M())
    compute_if_absent(histos,'pro_mntm_phimass_'+cut_name,builders['mass2DL']).Fill(pro.P(), vphi.M())
    compute_if_absent(histos,'phi_mntm_phimass_'+cut_name,builders['mass2DL']).Fill(vphi.P(), vphi.M())

    
    compute_if_absent(histos,'LepXvskpkm_'+cut_name,builders['mass2DL']).Fill(epX.M(), vphi.M())
    compute_if_absent(histos,'LepXvslmda_'+cut_name,builders['mass2DL']).Fill(epX.M(), lmda.M())
    compute_if_absent(histos,'LekpXvsprkm'+cut_name,builders['mass2DL']).Fill(ekpX.M(), lmda.M())
    compute_if_absent(histos,'LekpXvskpkm'+cut_name,builders['mass2DL']).Fill(ekpX.M(), vphi.M())
    compute_if_absent(histos,'LepkpXvsepkm'+cut_name,builders['mass2DL']).Fill(epkpX.M(), epkmX.M())
    
    compute_if_absent(histos,'cplpro_'+cut_name,builders['cpl']).Fill( cplpro )
    compute_if_absent(histos,'cplkp_'+cut_name,builders['cpl']).Fill( cplkp )
    compute_if_absent(histos,'cplkm_'+cut_name,builders['cpl']).Fill( cplkm )
    compute_if_absent(histos,'cplproVskpkm_'+cut_name,builders['kin2d']).Fill(vphi.M(), cplpro)
    compute_if_absent(histos,'cplkpVskpkm_'+cut_name,builders['kin2d']).Fill(vphi.M(), cplkp)
    compute_if_absent(histos,'cplkmVskpkm_'+cut_name,builders['kin2d']).Fill(vphi.M(), cplkm)

    compute_if_absent(histos,'pt_'+cut_name,builders['mm2']).Fill(pt)
    compute_if_absent(histos,'xq2_'+cut_name,builders['q2x']).Fill(xb,q2)
    compute_if_absent(histos,'txb_'+cut_name,builders['tx']).Fill(t,xb)
    compute_if_absent(histos,'q2t_'+cut_name,builders['kin2d']).Fill(q2,t)
    compute_if_absent(histos,'q2t_'+cut_name,builders['q2t']).Fill(q2,t)

    compute_if_absent(histos,'xb_el_theta_'+cut_name,builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
    compute_if_absent(histos,'xb_pr_theta_'+cut_name,builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
    compute_if_absent(histos,'xb_kp_theta_'+cut_name,builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
    compute_if_absent(histos,'xb_km_theta_'+cut_name,builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))

    compute_if_absent(histos,'xb_epkpkmXe_'+cut_name,builders['xbvsmm2']).Fill(epkpkmX.E(), xb)
    compute_if_absent(histos,'xb_ekpkmXe_'+cut_name,builders['xbvsmm2']).Fill(ekpkmX.M2(), xb)
    compute_if_absent(histos,'xb_epkpXe_'+cut_name,builders['xbvsmm2']).Fill(epkpX.M2(), xb)
    compute_if_absent(histos,'xb_epkmXe_'+cut_name,builders['xbvsmm2']).Fill(epkmX.M2(), xb)

    
    if( hel > 0 ):
        compute_if_absent(histos,'vphimass_pos_helicity_'+cut_name,builders['mass']).Fill(vphi.M())
    if( hel < 0 ):
        compute_if_absent(histos,'vphimass_neg_helicity_'+cut_name,builders['mass']).Fill(vphi.M())


    if( vphi.M() < 1.05 ):
        compute_if_absent(histos,'xq2_'+cut_name,builders['q2x']).Fill(xb,q2)
        compute_if_absent(histos,'txb_'+cut_name,builders['tx']).Fill(t,xb)
        compute_if_absent(histos,'q2t_'+cut_name,builders['kin2d']).Fill(q2,t)
        compute_if_absent(histos,'xb_el_theta_maxphi_limit_'+cut_name,builders['xeltheta']).Fill(xb, np.degrees(el.Theta()))
        compute_if_absent(histos,'xb_pr_theta_maxphi_limit_'+cut_name,builders['xprtheta']).Fill(xb, np.degrees(pro.Theta()))
        compute_if_absent(histos,'xb_kp_theta_maxphi_limit_'+cut_name,builders['xkptheta']).Fill(xb, np.degrees(kp.Theta()))
        compute_if_absent(histos,'xb_km_theta_maxphi_limit_'+cut_name,builders['xkmtheta']).Fill(xb, np.degrees(km.Theta()))

        
    compute_if_absent(histos,'cos_kp_'+cut_name,builders['cos']).Fill(cos_theta_cm_kp)    
    compute_if_absent(histos,'helicity_'+cut_name,builders['hel']).Fill(hel)
    



cut_names=['epkpkmmm2_cut', 'cplpro_cut', 'cplkm_cut', 'cplkp_cut', 'pt_cut', 'epkpXmm2_cut', 'ekpkmXmm2_cut', 'epkmXmm2_cut',
           'el_p_cut','pr_p_cut','kp_p_cut','km_p_cut','kp_theta_cut','km_theta_cut']

#no need to apply epX or kpkm cuts, both are equivalent but also restrict range of interest
excl_cut_names = ('missingE', 'epkpXmm2', 'ekpkmXmm2', 'epkmXmm2')
kin_cut_names = ('pr_p','kp_p','km_m')
    

#fx_cuts_names = ('pr_mntm_min','pr_mntm_max',

ccounter =0
for ev in ff.clas12:
    ccounter+=1
    #if( ccounter == 100 ): break

    el = TLorentzVector()
    el.SetPxPyPzE(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector()
    pro.SetPxPyPzE(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector()
    kp.SetPxPyPzE(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector()
    km.SetPxPyPzE(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)

    pr_status = ev.pr_status
    kp_status = ev.kp_status
    km_status = ev.km_status
    
    prb = ev.prb
    kpb = ev.kpb
    kmb = ev.kmb

    all_final_state_status_FD=False
    if( (pr_status >= 2000 and pr_status < 4000 and kp_status >= 2000 and kp_status<4000 and km_status >=2000 and km_status < 4000) ): 
        all_final_state_status_FD=True
        
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
    q = beam-el  
    q2 = -q.M2()  
    xb = q2/(2*target.M()*q.E())
    t = 2*target.M()*(pro.E()-target.M()) 
    nu = q2/(2*target.M()*xb)
    y = nu/beam.E() 

    v3l = beam.Vect().Cross(el.Vect())
    v3h = pro.Vect().Cross((beam-el).Vect()) 
    trento = np.degrees(TMath.ACos(v3l.Dot(v3h)/(v3l.Mag()*v3h.Mag())))
    if(v3l.Dot(pro.Vect())<0): trento=-trento
    
    boost_phicm = vphi.BoostVector()
    boost_phicm = -boost_phicm
    boost_kp = TLorentzVector(kp)
    boost_km = TLorentzVector(km)
    boost_kp.Boost(boost_phicm)
    boost_km.Boost(boost_phicm)
    cos_theta_cm_kp = TMath.Cos(boost_kp.Vect().Theta())
    # fill kinematic histograms here
    if( all_final_state_status_FD ):
                
        epkpXmm2 = epkpX.M2()
        epkmXmm2 = epkmX.M2()
        ekpkmXmm2 = ekpkmX.M2()
        epkpkmXmm2 = epkpkmX.M2()
                
        #exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6,  cplpro<30, cplkm<20, cplkp<20, pt<0.5, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
        ##                          el.P() > 1.5, pro.P() < 4, kp.P() < 2.5, km.P() < 2.5, np.degrees(kp.Theta())<35, np.degrees(km.Theta())<40)
                
                
        #OG CUT
        #exclusivity_cut_results = (np.absolute(epkpkmX.E()) < 0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
        '''
        #round 1 cuts
        epkpkmxe_min = 0.08 - 0.073*sig
        epkpkmxe_max = 0.08 + 0.073*sig

        epkpX_min = 0.246 - 0.065*sig
        epkpX_max = 0.246 + 0.065*sig

        epkmX_min = 0.257 - 0.0389*sig
        epkmX_max = 0.257 + 0.0389*sig

        ekpkmX_min = 0.908 - 0.0863*sig
        ekpkmX_max = 0.908 + 0.0863*sig
        '''

        #round 2 after fitting dist. with background - extracted fit parameters and reran code then fit again for final set
        
        exclusivity_cut_results = (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min , epkpXmm2 < epkpX_max and epkpXmm2 > epkpX_min, ekpkmXmm2 < ekpkmX_max and ekpkmXmm2 > ekpkmX_min, epkmXmm2 < epkmX_max and epkmXmm2 > epkmX_min)#, (epX.M2() > 0.8) and (epX.M2()<1.2), (vphi.M() < 1.1) )
        
        kin_cut_results = (el.P() > 0.5, pro.P() < 4, kp.P() < 3.5, km.P() < 3.5)
        fx_cut_results = ( ekpkmXmm2 > 0.8 and ekpkmXmm2 < 1.15, epkpXmm2 > 0.08 and epkpXmm2 < 0.48, epkmXmm2 > 0.08 and epkmXmm2 < 0.48, cplpro < 6.0, cplkm < 6.0, cplkm < 6.0 , pt  < 0.120, epkpkmX.E() > -0.175 and epkpkmX.E() < 0.320, epkpkmXmm2 < 0.0075  )

        excl_cut_results = (  excl_cuts.passEPKPKMX(epkpkmX),
                              excl_cuts.passEPKPX(epkpX),
                              excl_cuts.passEKPKMX(ekpkmX),
                              excl_cuts.passEPKMX(epkmX)
                          )
                              
        exclusivity_cut_results2 = exclusivity_cut_results
            
        # (epkpkmX.E() < 0.4 and epkpkmX.E() >-0.2 , epkpXmm2 < 0.5 and epkpXmm2 > 0.1, ekpkmXmm2 < 1.2 and ekpkmXmm2 > 0.65, epkmXmm2 < 0.5 and epkmXmm2 > 0.05 )#, (epX.M2() > 0.8) and (epX.M2()<1.2))
        
        remove_lambda = lmda.M() < 1.49 or lmda.M() > 1.587
        remove_lambda2 = lmda.M() < 1.75 or lmda.M() > 1.89
        select_phi = vphi.M() < max_phi_mass_signal_limit and vphi.M() > min_phi_mass_signal_limit 

        coplanarity_cut_results = (cplpro<2.5, cplkm<2.5, cplkp<2.5)

        compute_if_absent(histos,'pro_beta_vs_p_raw',builders['betap']).Fill(pro.P(),prb)
        compute_if_absent(histos,'kp_beta_vs_p_raw',builders['betap']).Fill(kp.P(),kpb)
        compute_if_absent(histos,'km_beta_vs_p_raw',builders['betap']).Fill(km.P(),kmb)
        

        plotHistos(el,pro,kp,km,ev.hel,'raw')    
        if(exclusivity_cut_results[0] and exclusivity_cut_results[1] and exclusivity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_c012')
        if(exclusivity_cut_results[0] and exclusivity_cut_results[1] and exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_c013')
        if(exclusivity_cut_results[0] and exclusivity_cut_results[2] and exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_c023')
        if(exclusivity_cut_results[1] and exclusivity_cut_results[2] and exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_c123')
                

        #plot those events that fail c012,etc.
        if(not exclusivity_cut_results[0] and not exclusivity_cut_results[1] and not exclusivity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'fail_c012')    
        if(not exclusivity_cut_results[0] and not exclusivity_cut_results[1] and not exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'fail_c013')
        if(not exclusivity_cut_results[0] and not exclusivity_cut_results[2] and not exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'fail_c023')
        if(not exclusivity_cut_results[1] and not exclusivity_cut_results[2] and not exclusivity_cut_results[3]):
            plotHistos(el,pro,kp,km,ev.hel,'fail_c123')

        if( all(exclusivity_cut_results) and coplanarity_cut_results[0] and coplanarity_cut_results[1]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_excl_c01')
             
        if( all(exclusivity_cut_results) and coplanarity_cut_results[1] and coplanarity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_excl_c12')
            
        if( all(exclusivity_cut_results) and coplanarity_cut_results[0] and coplanarity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_excl_c02')

        if( coplanarity_cut_results[0] and coplanarity_cut_results[1]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_coplan_c01')
             
        if( coplanarity_cut_results[1] and coplanarity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_coplan_c12')
            
        if( coplanarity_cut_results[0] and coplanarity_cut_results[2]):
            plotHistos(el,pro,kp,km,ev.hel,'pass_coplan_c02')

        if( all(coplanarity_cut_results) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_coplan')


        coplancc_indx=0
        for coplancc in coplanarity_cut_results:
            if( coplancc ):
                plotHistos(el,pro,kp,km,ev.hel,'pass_coplancc'+str(coplancc_indx))
            coplancc_indx+=1
    
        cc_idx=0
        any_false=True
        for cc in exclusivity_cut_results :
            if(cc):
                plotHistos(el,pro,kp,km,ev.hel,excl_cut_names[cc_idx])
            cc_idx+=1
            if( not cc ):
                any_false=False


        ### study cut combos            
        if( all(exclusivity_cut_results2) and all(coplanarity_cut_results)):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_pass_all_coplanarity')

        if( all(kin_cut_results) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_kinematic_cuts')

        if( all(fx_cut_results) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_FX_cuts')

        if passEPKPKMX(epkpkmX) :
            plotHistos(el,pro,kp,km,ev.hel,'pass_epkpkmX_cut_only')
            
        if passEPKPX(epkpX) :
            plotHistos(el,pro,kp,km,ev.hel,'pass_epkpX_cut_only')

        if passEPKMX(epkmX) :
            plotHistos(el,pro,kp,km,ev.hel,'pass_epkmX_cut_only')

        if passEKPKMX(ekpkmX) :
            plotHistos(el,pro,kp,km,ev.hel,'pass_ekpkmX_cut_only')

        
        if( all(exclusivity_cut_results2) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity')
            if select_phi:
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_select_phi')
            if select_phi and remove_lambda:
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_select_phi_rmv_lambda')

            if remove_lambda :
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_phi_rmv_lambda')

            if remove_lambda and remove_lambda2:
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_rmv_both_lambdas')
                
            if select_phi and remove_lambda and pro.P() < 3.0:
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_select_phi_rmv_lambda_pro_mntm_cut')



            compute_if_absent(histos,'pro_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(pro.P(),prb)
            compute_if_absent(histos,'kp_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(kp.P(),kpb)
            compute_if_absent(histos,'km_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(km.P(),kmb)
            
            if( remove_lambda and remove_lambda2):
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_remove_lambda')
                compute_if_absent(histos,'pro_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(pro.P(),prb)
                compute_if_absent(histos,'kp_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(kp.P(),kpb)
                compute_if_absent(histos,'km_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(km.P(),kmb)

                if( vphi.M() < max_phi_mass_signal_limit and vphi.M() > min_phi_mass_signal_limit):
                    if( ev.hel > 0 ):
                        compute_if_absent(histos,'pos_asy_pass_all_exclusivity_remove_lmbd_vphim',builders['asy']).Fill(trento)
                    if( ev.hel < 0 ):            
                            compute_if_absent(histos,'neg_asy_pass_all_exclusivity_remove_lmda_vphim',builders['asy']).Fill(trento)

            if( vphi.M() < max_phi_mass_signal_limit and vphi.M() > min_phi_mass_signal_limit):
                if( ev.hel > 0 ):
                    compute_if_absent(histos,'pos_asy_pass_all_exclusivity_vphim',builders['asy']).Fill(trento)
                if( ev.hel < 0 ):            
                    compute_if_absent(histos,'neg_asy_pass_all_exclusivity_vphim',builders['asy']).Fill(trento)


            coplancc_indx=0
            for coplancc in coplanarity_cut_results:
                if(coplancc):
                    plotHistos(el,pro,kp,km,ev.hel,'pass_all_excl_pass_coplancc'+str(coplancc_indx))
                coplancc_indx+=1


        if( all(exclusivity_cut_results2) and (lmda.M() > 1.6 or lmda.M() < 1.49 ) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_pass_lmdacut')


        ### fail all cuts
        result_cut012=returnListResult(3,exclusivity_cut_results2)
        result_cut123=returnListResult(0,exclusivity_cut_results2)
        result_cut023=returnListResult(1,exclusivity_cut_results2)
        result_cut013=returnListResult(2,exclusivity_cut_results2)
        

        if( False in exclusivity_cut_results2 ):
            #plotHistos(el,pro,kp,km,ev.hel,'fail_at_least_one_exclusivity')
            if( ev.hel > 0 ):
                compute_if_absent(histos,'pos_asy_fail_all_exclusivity',builders['asy']).Fill(trento)
            if( ev.hel < 0 ):            
                compute_if_absent(histos,'neg_asy_fail_all_exclusivity',builders['asy']).Fill(trento)

        if( not result_cut012 ):
            plotHistos(el,pro,kp,km,ev.hel,'fail_at_least_c012')
        if( not result_cut123 ):
            plotHistos(el,pro,kp,km,ev.hel,'fail_at_least_c123')
        if( not result_cut023 ):
            plotHistos(el,pro,kp,km,ev.hel,'fail_at_least_c023')
        if( not result_cut013 ):
            plotHistos(el,pro,kp,km,ev.hel,'fail_at_least_c013')




fout = TFile('exclusive_phi_kinematics_final_ebc_heb_FDonly_'+tag+'_'+base,'recreate')
for _, hist in histos.items():
    hist.Write()

fout.Close()
