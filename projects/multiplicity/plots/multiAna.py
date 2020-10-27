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

pr_multi = {}
kp_multi = {}
km_multi = {}

pr_multi_fd={}
pr_multi_cd={}

kp_multi_fd={}
kp_multi_cd={}

km_multi_fd={}
km_multi_cd={}



def buildMultiplicity(title):
    print('building {}'.format(title))
    return TH1F(title,title, 19, 1.0, 20.0)

builders={
    'multiplicity' : buildMultiplicity
}

event_comb = []
event_comb_cd = []
event_comb_fd = []
event_comb_final = []
event_comb_final_cd = []
event_comb_final_fd = []


cc=0
prev_event=0
for ev in ff.clas12:
    cc+=1
    if(cc >= 10000 ):break 

    el_indx = ev.el_indx
    pr_indx = ev.pr_indx
    kp_indx = ev.kp_indx
    km_indx = ev.km_indx

    pr_stat = ev.pr_status
    kp_stat = ev.kp_status
    km_stat = ev.km_status

    glob_event = ev.global_event
    comb_indx  = ev.comb_index

    if( pr_stat >= 2000 and pr_stat < 4000 ):
        event_comb_fd.append(glob_event)
    else:
        event_comb_cd.append(glob_event)

    el = TLorentzVector()
    el.SetPxPyPzE(ev.el_px,ev.el_py,ev.el_pz,ev.el_e)
    pro = TLorentzVector()
    pro.SetPxPyPzE(ev.pr_px,ev.pr_py,ev.pr_pz,ev.pr_e)
    kp = TLorentzVector()
    kp.SetPxPyPzE(ev.kp_px,ev.kp_py,ev.kp_pz,ev.kp_e)
    km = TLorentzVector()
    km.SetPxPyPzE(ev.km_px,ev.km_py,ev.km_pz,ev.km_e)

    #print('proton %f and kaonP %f and kaonM %f  ' % (pro.P(), kp.P() , km.P()))
    #print(' current event %d and previous event %d' % (glob_event,prev_event))
        
    if(  glob_event > prev_event):
        prev_event=glob_event
        pr_multi[glob_event]=set()
        kp_multi[glob_event]=set()
        km_multi[glob_event]=set()
    
        pr_multi_fd[glob_event]=set()
        kp_multi_fd[glob_event]=set()
        km_multi_fd[glob_event]=set()
    
        pr_multi_cd[glob_event]=set()
        kp_multi_cd[glob_event]=set()
        km_multi_cd[glob_event]=set()
    if( glob_event == prev_event ):
        #print('fill set here %d ' % (glob_event))
        pr_multi[glob_event].add(pro.P())
        kp_multi[glob_event].add(kp.P())
        km_multi[glob_event].add(km.P())

        #fill pr
        if( pr_stat >= 2000 and pr_stat < 4000 ):
            pr_multi_fd[glob_event].add(pro.P())
        else:
            pr_multi_cd[glob_event].add(pro.P())
        #fill kp
        if( kp_stat >= 2000 and kp_stat < 4000):
            kp_multi_fd[glob_event].add(kp.P())
        else:
            kp_multi_cd[glob_event].add(kp.P())        

        #fill km
        if( km_stat >= 2000 and km_stat < 4000 ):
            km_multi_fd[glob_event].add(km.P())
        else:
            km_multi_cd[glob_event].add(km.P())
    
        #print(' after current event %d and previous event %d' % (glob_event,prev_event))
    
    

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
    event_comb.append(glob_event)

    #exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6,  cplpro<30, cplkm<20, cplkp<20, pt<0.5, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)


    exclusivity_cut_results = (np.absolute(epkpkmX.E()) < 0.4 and np.absolute(epkpkmX.E()) >-0.2 , epkpXmm2 < 0.5 and epkpXmm2 > 0.1, ekpkmXmm2 < 1.2 and ekpkmXmm2 > 0.65, epkmXmm2 < 0.5 and epkmXmm2 > 0.05)



    if( all(exclusivity_cut_results) ):
        event_comb_final.append(glob_event)

        if( pr_stat >= 2000 and pr_stat < 4000 ):
            event_comb_final_fd.append(glob_event)
        else:
            event_comb_final_cd.append(glob_event)
        


proton_counts = [ len(value) for key,value in pr_multi.items() ]
proton_counts_fd = [ len(value) for key,value in pr_multi_fd.items() ]
proton_counts_cd = [ len(value) for key,value in pr_multi_cd.items() ]

kp_counts = [ len(value) for key,value in kp_multi.items() ]
kp_counts_fd = [ len(value) for key,value in kp_multi_fd.items() ]
kp_counts_cd = [ len(value) for key,value in kp_multi_cd.items() ]

km_counts = [ len(value) for key,value in km_multi.items() ]
km_counts_fd = [ len(value) for key,value in km_multi_fd.items() ]
km_counts_cd = [ len(value) for key,value in km_multi_cd.items() ]



for ii in proton_counts:
    compute_if_absent(histos,'proton_counts_all_events',builders['multiplicity']).Fill(ii) 
for ii in proton_counts_fd:
    compute_if_absent(histos,'proton_counts_all_events_FD',builders['multiplicity']).Fill(ii) 
for ii in proton_counts_cd:
    compute_if_absent(histos,'proton_counts_all_events_CD',builders['multiplicity']).Fill(ii) 

for ii in kp_counts:
    compute_if_absent(histos,'kp_counts_all_events',builders['multiplicity']).Fill(ii) 
for ii in kp_counts_fd:
    compute_if_absent(histos,'kp_counts_all_events_FD',builders['multiplicity']).Fill(ii) 
for ii in kp_counts_cd:
    compute_if_absent(histos,'kp_counts_all_events_CD',builders['multiplicity']).Fill(ii) 

for ii in km_counts:
    compute_if_absent(histos,'km_counts_all_events',builders['multiplicity']).Fill(ii) 
for ii in km_counts_fd:
    compute_if_absent(histos,'km_counts_all_events_FD',builders['multiplicity']).Fill(ii) 
for ii in km_counts_cd:
    compute_if_absent(histos,'km_counts_all_events_CD',builders['multiplicity']).Fill(ii) 


m_event_comb= dict((x,event_comb.count(x)) for x in set(event_comb))
m_event_comb_cd= dict((x,event_comb_cd.count(x)) for x in set(event_comb_cd))
m_event_comb_fd= dict((x,event_comb_fd.count(x)) for x in set(event_comb_fd))

m_event_comb_final_temp= dict((x,event_comb_final.count(x)) for x in set(event_comb_final))
m_event_comb_final_cd_temp= dict((x,event_comb_final_cd.count(x)) for x in set(event_comb_final_cd))
m_event_comb_final_fd_temp= dict((x,event_comb_final_fd.count(x)) for x in set(event_comb_final_fd))

#get the number of combinations for an event that passed by finding it in the uncut dictionary
m_event_comb_final = { key:value for key,value in m_event_comb.items() if key in m_event_comb_final_temp } 
m_event_comb_final_cd = { key:value for key,value in m_event_comb.items() if key in m_event_comb_final_cd_temp } 
m_event_comb_final_fd = { key:value for key,value in m_event_comb.items() if key in m_event_comb_final_fd_temp } 

event_values = list((m_event_comb.values()))
event_values_fd = list((m_event_comb_fd.values()))
event_values_cd = list((m_event_comb_cd.values()))

event_values_final = list((m_event_comb_final.values()))
event_values_final_fd = list((m_event_comb_final_fd.values()))
event_values_final_cd = list((m_event_comb_final_cd.values()))


for ii in event_values:
    compute_if_absent(histos,'all_combination_all_events',builders['multiplicity']).Fill(ii)

for ii in event_values_cd:
    compute_if_absent(histos,'CD_combination_all_events',builders['multiplicity']).Fill(ii)

for ii in event_values_fd:
    compute_if_absent(histos,'FD_combination_all_events',builders['multiplicity']).Fill(ii)

for ii in event_values_final:
    compute_if_absent(histos,'all_combination_final_events',builders['multiplicity']).Fill(ii)

for ii in event_values_final_cd:
    compute_if_absent(histos,'CD_combination_final_events',builders['multiplicity']).Fill(ii)

for ii in event_values_final_fd:
    compute_if_absent(histos,'FD_combination_final_events',builders['multiplicity']).Fill(ii)



fout = TFile('phi_beta_study'+base,'recreate')
for _, hist in histos.items():
    hist.Write()

fout.Close()
ff.Close()
