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

#######################################################
##### MACRO TO SELECT EVENTS FOR SIDEBAND SUBTRACTION 

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


phi_mass_fit=1.021 #1.020
phi_mass_resolution = 5.64e-03  
# 2.5 defined by signal to background ration of about 0.7
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution
hmass = TH1F('hmass','hmass',100,0.8,1.5)
#bin_phimass_cut = hmass.GetXaxis().FindBin(phi_mass_signal_limit)
max_bin_phimass_cut = hmass.GetXaxis().FindBin(max_phi_mass_signal_limit)
min_bin_phimass_cut = hmass.GetXaxis().FindBin(min_phi_mass_signal_limit)


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
    return TH1F(title,title,20,-180.0,180.0)
def buildHel(title):
    print('building {}'.format(title))
    return TH1F(title,title,4,-2,2)
def buildQ2x(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,1,100,0.0,beam.E())
def buildtx(title):
    print('building {}'.format(title))
    return TH2F(title,title,100,0,6,100,0.0,1.0)

def buildBetaP(title):
    print('building {}'.format(title))
    return TH2F(title,title, 200, 0.0, beam.E(), 200, 0.0, 1.25)

    
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
    'tx':buildtx,
    'betap': buildBetaP
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
    
    boost_phicm = vphi.BoostVector()
    boost_phicm = -boost_phicm
    boost_kp = TLorentzVector(kp)
    boost_km = TLorentzVector(km)
    boost_kp.Boost(boost_phicm)
    boost_km.Boost(boost_phicm)
    cos_theta_cm_kp = TMath.Cos(boost_kp.Vect().Theta())

    kp_as_pip = TLorentzVector(kp.Px(), kp.Py(), kp.Pz(), np.sqrt(kp.P()**2 + 0.1396**2))
    vphi_kp_as_pip = kp_as_pip + km
    
    compute_if_absent(histos,'vphimass_'+cut_name,builders['mass']).Fill(vphi.M())
    compute_if_absent(histos,'lambdamass_'+cut_name,builders['massL']).Fill(lmda.M())
    compute_if_absent(histos,'lambdaphimass_'+cut_name,builders['mass2D']).Fill(lmda.M(), vphi.M())

    if( hel > 0 ):
        compute_if_absent(histos,'vphimass_pos_helicity_'+cut_name,builders['mass']).Fill(vphi.M())
    if( hel < 0 ):
        compute_if_absent(histos,'vphimass_neg_helicity_'+cut_name,builders['mass']).Fill(vphi.M())
    


cut_names=['epkpkmmm2_cut', 'cplpro_cut', 'cplkm_cut', 'cplkp_cut', 'pt_cut', 'epkpXmm2_cut', 'ekpkmXmm2_cut', 'epkmXmm2_cut',
           'el_p_cut','pr_p_cut','kp_p_cut','km_p_cut','kp_theta_cut','km_theta_cut']

#no need to apply epX or kpkm cuts, both are equivalent but also restrict range of interest
excl_cut_names = ('missingE', 'epkpXmm2', 'ekpkmXmm2', 'epkmXmm2')
    
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
                
        #exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6,  cplpro<30, cplkm<20, cplkp<20, pt<0.5, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
        ##                          el.P() > 1.5, pro.P() < 4, kp.P() < 2.5, km.P() < 2.5, np.degrees(kp.Theta())<35, np.degrees(km.Theta())<40)
                
                
        #OG CUT
        #exclusivity_cut_results = (np.absolute(epkpkmX.E()) < 0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
        sig=5
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
        epkpkmxe_min = 0.08 - 0.0398*sig
        epkpkmxe_max = 0.08 + 0.0398*sig

        epkpX_min = 0.248 - 0.059*sig
        epkpX_max = 0.248 + 0.059*sig

        epkmX_min = 0.248 - 0.055*sig
        epkmX_max = 0.248 + 0.055*sig

        ekpkmX_min = 0.904 - 0.0719*sig
        ekpkmX_max = 0.904 + 0.0719*sig
        
        exclusivity_cut_results = (epkpkmX.E() < epkpkmxe_max and epkpkmX.E() > epkpkmxe_min , epkpXmm2 < epkpX_max and epkpXmm2 > epkpX_min, ekpkmXmm2 < ekpkmX_max and ekpkmXmm2 > ekpkmX_min, epkmXmm2 < epkmX_max and epkmXmm2 > epkmX_min)#, (epX.M2() > 0.8) and (epX.M2()<1.2), (vphi.M() < 1.1) )
        
        exclusivity_cut_results2 = exclusivity_cut_results

        # (epkpkmX.E() < 0.4 and epkpkmX.E() >-0.2 , epkpXmm2 < 0.5 and epkpXmm2 > 0.1, ekpkmXmm2 < 1.2 and ekpkmXmm2 > 0.65, epkmXmm2 < 0.5 and epkmXmm2 > 0.05 )#, (epX.M2() > 0.8) and (epX.M2()<1.2))
        
        remove_lambda = lmda.M() < 1.49 or lmda.M() > 1.587
        remove_lambda2 = lmda.M() < 1.75 or lmda.M() > 1.89

        coplanarity_cut_results = (cplpro<2.5, cplkm<2.5, cplkp<2.5)

        plotHistos(el,pro,kp,km,ev.hel,'raw')    
           
       
        ### study cut combos            
        if( all(exclusivity_cut_results2) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity')

            compute_if_absent(histos,'pro_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(pro.P(),prb)
            compute_if_absent(histos,'kp_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(kp.P(),kpb)
            compute_if_absent(histos,'km_beta_vs_p_pass_all_exclusivity',builders['betap']).Fill(km.P(),kmb)
            
            if( remove_lambda and remove_lambda2 ):
                plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_remove_lambda')
                compute_if_absent(histos,'pro_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(pro.P(),prb)
                compute_if_absent(histos,'kp_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(kp.P(),kpb)
                compute_if_absent(histos,'km_beta_vs_p_pass_all_exclusivity_remove_lambda',builders['betap']).Fill(km.P(),kmb)

                if( vphi.M() < 1.05 ):
                    if( ev.hel > 0 ):
                        compute_if_absent(histos,'pos_asy_pass_all_exclusivity_remove_lmbd_vphim',builders['asy']).Fill(trento)
                    if( ev.hel < 0 ):            
                            compute_if_absent(histos,'neg_asy_pass_all_exclusivity_remove_lmda_vphim',builders['asy']).Fill(trento)

            if( vphi.M() < 1.05 ):
                if( ev.hel > 0 ):
                    compute_if_absent(histos,'pos_asy_pass_all_exclusivity_vphim',builders['asy']).Fill(trento)
                if( ev.hel < 0 ):            
                    compute_if_absent(histos,'neg_asy_pass_all_exclusivity_vphim',builders['asy']).Fill(trento)

            if( exclusivity_cut_results2 and  remove_lambda and remove_lambda2 ):
        
                bin_phimass = hmass.GetXaxis().FindBin(vphi.M())
                max_bin_min_phi = hmass.GetXaxis().FindBin(max_phi_mass_signal_limit)
                min_bin_min_phi = hmass.GetXaxis().FindBin(min_phi_mass_signal_limit)
                min_bin_mass = max_bin_min_phi
                for ii in range(1,25):

                    max_bin_mass = min_bin_mass + 3#13
                    max_phi_mass=0
                    #print(' min %f , max %f ' % (min_mass, max_mass))


                    if( bin_phimass >= min_bin_mass and bin_phimass < max_bin_mass ):
                        compute_if_absent(histos,'vphimass_mass_binrange_'+str(ii),builders['mass']).Fill(vphi.M())

                        if( ev.hel > 0 ):
                            #print('here')
                            compute_if_absent(histos,'pos_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii),builders['asy']).Fill(trento)
                            #compute_if_absent(histos,'pos_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii),builders['asy']).Fill(trento)
                        
                        if( ev.hel < 0 ):            
                            compute_if_absent(histos,'neg_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii),builders['asy']).Fill(trento)
                            #compute_if_absent(histos,'neg_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii),builders['asy']).Fill(trento)
                    min_bin_mass = max_bin_mass

        
            
                

            


        if( all(exclusivity_cut_results2) and (lmda.M() > 1.6 or lmda.M() < 1.49 ) ):
            plotHistos(el,pro,kp,km,ev.hel,'pass_all_exclusivity_pass_lmdacut')


        ### fail all cuts
        result_cut012=returnListResult(3,exclusivity_cut_results2)
        result_cut123=returnListResult(0,exclusivity_cut_results2)
        result_cut023=returnListResult(1,exclusivity_cut_results2)
        result_cut013=returnListResult(2,exclusivity_cut_results2)
        
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
    
    
