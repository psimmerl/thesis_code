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

phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution

sig=4.0
epkpkmxe_min = 0.08 - 0.0398*sig
epkpkmxe_max = 0.08 + 0.0398*sig

epkpX_min = 0.248 - 0.059*sig
epkpX_max = 0.248 + 0.059*sig

epkmX_min = 0.248 - 0.055*sig
epkmX_max = 0.248 + 0.055*sig

ekpkmX_min = 0.904 - 0.0719*sig
ekpkmX_max = 0.904 + 0.0719*sig

def passEvent( epkpkmX, epkpX, epkmX, ekpkmX ):
    exclusivity_cut_results = (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min , epkpX.M2() < epkpX_max and epkpX.M2() > epkpX_min, ekpkmX.M2() < ekpkmX_max and ekpkmX.M2() > ekpkmX_min, epkmX.M2() < epkmX_max and epkmX.M2() > epkmX_min)    
    return exclusivity_cut_results


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

def buildMCResMntmPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiPr(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)

def buildMCResMntmKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiKp(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)

def buildMCResMntmKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 4.0, 250, -0.15, 0.15)
def buildMCResThetaKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, 0.0, 60.0, 250, -2.0, 2.0)
def buildMCResPhiKm(title):
    print('building {}'.format(title))
    return TH2F(title,title,250, -180, 180, 250, -2.0, 2.0)

def buildQ2(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0,beam.E())
def buildx(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0.0,1.05)
def buildt(title):
    print('building {}'.format(title))
    return TH1F(title,title,100,0,5.0)


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
    'mcresthetael': buildMCResThetaEl,
    'mcresphiel': buildMCResPhiEl,

    'mcresmntmpr': buildMCResMntmPr,
    'mcresthetapr': buildMCResThetaPr,
    'mcresphipr': buildMCResPhiPr,

    'mcresmntmkp': buildMCResMntmKp,
    'mcresthetakp': buildMCResThetaKp,
    'mcresphikp': buildMCResPhiKp,

    'mcresmntmkm': buildMCResMntmKm,
    'mcresthetakm': buildMCResThetaKm,
    'mcresphikm': buildMCResPhiKm,

    'q2':buildQ2,
    't':buildt,
    'x':buildx    


}

    
for ev in ff.clas12:

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

    compute_if_absent(histos,'q2_prior_cuts',builders['q2']).Fill(q2)
    compute_if_absent(histos,'t_prior_cuts',builders['t']).Fill(t)
    compute_if_absent(histos,'xb_prior_cuts',builders['x']).Fill(xb)


    if( (pr_status >= 2000 and pr_status < 4000 and kp_status >= 2000 and kp_status<4000 and km_status >=2000 and km_status < 4000) ): 
        all_final_state_status_FD=True
            compute_if_absent(histos,'q2_prior_cuts_FD',builders['q2']).Fill(q2)
            compute_if_absent(histos,'t_prior_cuts_FD',builders['t']).Fill(t)
            compute_if_absent(histos,'xb_prior_cuts_FD',builders['x']).Fill(xb)


        if( all( passEvent(epkpkmX, epkpX, epkmX, ekpkmX) ) ):

            q = beam-el  
            q2 = -q.M2()  
            xb = q2/(2*target.M()*q.E())
            t = 2*target.M()*(pro.E()-target.M()) 

            ###############
            ##kin
            compute_if_absent(histos,'q2',builders['q2']).Fill(q2)
            compute_if_absent(histos,'t',builders['t']).Fill(t)
            compute_if_absent(histos,'xb',builders['x']).Fill(xb)



            ################
            ## epX 
            compute_if_absent(histos,'epX_vs_theta_pr',builders['mmthetaEPX']).Fill(np.degrees(pro.Theta()), epX.M2())
            compute_if_absent(histos,'epX_vs_p_pr',builders['mmpEPX']).Fill(pro.P(), epX.M2())
            
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
            delta_p_pr =  pro.P() - ekpkmX.P() 
            delta_theta_pr = np.degrees(-ekpkmX.Theta() + pro.Theta())
            
            delta_p_km = -epkpX.P() + km.P()
            delta_theta_km = np.degrees(-epkpX.Theta() + km.Theta())
            
            delta_p_kp = -epkmX.P() + kp.P()
            delta_theta_kp = np.degrees(-epkmX.Theta() + kp.Theta())
            
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


        
        
        
        
        
fout = TFile('data_res_phi_kin_ebc_heb_'+tag+'.root','recreate')
for _, hist in histos.items():
    hist.Write()

fout.Close()



    
    
