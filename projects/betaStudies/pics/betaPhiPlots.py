from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile
from ROOT import TGraphErrors
from ROOT import TLorentzVector
from ROOT import gROOT, gBenchmark, gStyle, gPad
from collections import OrderedDict
from array import array
import numpy as np
import math
import sys, os

beam = TLorentzVector(0,0,10.604,10.604)
target = TLorentzVector(0,0,0.0,0.938272)

ff = TFile(sys.argv[1])
datatype = sys.argv[2]

base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )


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
    return TH2F(title,title, 200, -180.0, 180.0, 200, 0.0, 60.0)
def buildBetaP(title):
    print('building {}'.format(title))
    return TH2F(title,title, 200, 0.0, beam.E(), 200, 0.0, 1.25)

builders={
    'ptheta':buildPTheta,
    'phitheta':buildPhiTheta,
    'betap': buildBetaP
}

for ev in ff.clas12:
    
    el_indx = ev.el_indx
    pr_indx = ev.pr_indx
    kp_indx = ev.kp_indx
    km_indx = ev.km_indx

    pr_stat = ev.pr_status
    kp_stat = ev.kp_status
    km_stat = ev.km_status

    glob_event = ev.global_event
    comb_indx  = ev.comb_index
    
    el = TLorentzVector()
    el.SetPxPyPzE(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector()
    pro.SetPxPyPzE(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector()
    kp.SetPxPyPzE(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector()
    km.SetPxPyPzE(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)
     
    prb = ev.prb
    kpb = ev.kpb
    kmb = ev.kmb

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
    
    lmda = pro+kp 
    vphi = kp+km 

    epkpXmm2 = epkpX.M2()
    epkpkmXmm2 = epkpkmX.M2()
    epkmXmm2 = epkmX.M2()
    ekpkmXmm2 = ekpkmX.M2()
    
    #print(' global event %d - combination %d' % (glob_event, comb_indx) )
    
    compute_if_absent(histos,'el_ptheta',builders['ptheta']).Fill( el.P(),  np.degrees(el.Theta()) )
    compute_if_absent(histos,'pr_ptheta',builders['ptheta']).Fill( pro.P(),  np.degrees(pro.Theta()) )
    compute_if_absent(histos,'kp_ptheta',builders['ptheta']).Fill( kp.P(),  np.degrees(kp.Theta()) )
    compute_if_absent(histos,'km_ptheta',builders['ptheta']).Fill( km.P(),  np.degrees(km.Theta()) )
    
    compute_if_absent(histos,'el_phitheta',builders['phitheta']).Fill( np.degrees(el.Phi()),  np.degrees(el.Theta()) )
    compute_if_absent(histos,'pr_phitheta',builders['phitheta']).Fill(np.degrees(pro.Phi()),  np.degrees(pro.Theta()) )
    compute_if_absent(histos,'kp_phitheta',builders['phitheta']).Fill( np.degrees(kp.Phi()),  np.degrees(kp.Theta()) )
    compute_if_absent(histos,'km_phitheta',builders['phitheta']).Fill( np.degrees(km.Phi()),  np.degrees(km.Theta()) )

    exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6,  cplpro<30, cplkm<20, cplkp<20, pt<0.5, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)

    if( pr_stat >= 2000 and pr_stat < 4000 ): 
            compute_if_absent(histos,'pr_betap',builders['betap']).Fill(pro.P(), prb)
    

    if( kp_stat >= 2000 and kp_stat < 4000 ):
            compute_if_absent(histos,'kp_betap',builders['betap']).Fill(kp.P(), kpb)
    

    if( km_stat >= 2000 and km_stat < 4000 ): 
            compute_if_absent(histos,'km_betap',builders['betap']).Fill(km.P(), kmb)
    

    if( all(exclusivity_cut_results) ):

        compute_if_absent(histos,'el_ptheta_final',builders['ptheta']).Fill( el.P(),  np.degrees(el.Theta()) )
        compute_if_absent(histos,'pr_ptheta_final',builders['ptheta']).Fill( pro.P(),  np.degrees(pro.Theta()) )
        compute_if_absent(histos,'kp_ptheta_final',builders['ptheta']).Fill( kp.P(),  np.degrees(kp.Theta()) )
        compute_if_absent(histos,'km_ptheta_final',builders['ptheta']).Fill( km.P(),  np.degrees(km.Theta()) )

        compute_if_absent(histos,'el_phitheta_final',builders['phitheta']).Fill( np.degrees(el.Phi()),  np.degrees(el.Theta()) )
        compute_if_absent(histos,'pr_phitheta_final',builders['phitheta']).Fill(np.degrees(pro.Phi()),  np.degrees(pro.Theta()) )
        compute_if_absent(histos,'kp_phitheta_final',builders['phitheta']).Fill( np.degrees(kp.Phi()),  np.degrees(kp.Theta()) )
        compute_if_absent(histos,'km_phitheta_final',builders['phitheta']).Fill( np.degrees(km.Phi()),  np.degrees(km.Theta()) )
        
        
        if( pr_stat >= 2000 and pr_stat < 4000 ): 
            compute_if_absent(histos,'pr_betap_final',builders['betap']).Fill(pro.P(), prb)
        
            
        if( kp_stat >= 2000 and kp_stat < 4000 ):
            compute_if_absent(histos,'kp_betap_final',builders['betap']).Fill(kp.P(), kpb)
            
            
        if( km_stat >= 2000 and km_stat < 4000 ): 
            compute_if_absent(histos,'km_betap_final',builders['betap']).Fill(km.P(), kmb)
                
                
            
    

    
fout = TFile('phi_beta_study'+base,'recreate')
for _, hist in histos.items():
    hist.Write()

fout.Close()
