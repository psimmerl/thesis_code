from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from collections import OrderedDict

import numpy as np
import os, sys

from array import array
import math

def readInBinningInfo(fname):
    fin = open(fname,'r')
    temp = []
    for ll in fin:
        print ll
        temp.append(float(ll))
    return temp    

def readInOverlapBinningInfo(fname):
    ovrlp_list = []
    fin = open(fname,'r')
    for ll in fin:
        temp_ll = ll.split(' ')
        print(temp_ll)
        ovrlp_list.append([float(temp_ll[0]), float(temp_ll[1])])
    return ovrlp_list

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
def buildMass(title):
    print('building {}'.format(title))
    return TH1F(title,title,150,0.8,1.5)
def buildMassSmallBins(title):
    print('building {}'.format(title))
    return TH1F(title,title,300,0.8,1.5)
def buildAsy(title):
    print('building {}'.format(title))
    return TH1F(title,title,10,-180.0,180.0)
def buildAsyMod2(title):
    print('building {}'.format(title))
    return TH1F(title,title,12,-180.0,180.0)
def buildAsyMod3(title):
    bin_low_trento=[-180, -144, -108, -72, -36, 36, 72, 108, 144, 180]
    print('building {}'.format(title))
    return TH1F(title,title,len(bin_low_trento)-1,array('d',bin_low_trento))
def buildAsyMod4(title):
    bin_low_trento=[-180, -150, -120, -90, -60, -30, 30, 60, 90, 120, 150, 180]
    print('building {}'.format(title))
    return TH1F(title,title,len(bin_low_trento)-1,array('d',bin_low_trento))
def buildQ2(title):
    print('building {}'.format(title))
    return TH1F(title,title,100, 0.0, 12)
def buildXb(title):
    print('building {}'.format(title))
    return TH1F(title,title,100, 0.0, 1.15)
def buildT(title):
    print('building {}'.format(title))
    return TH1F(title,title,100, 0.0, 6)
def buildW(title):
    print('building {}'.format(title))
    return TH1F(title,title,100, 1.5, 4)
def buildXbT(title):
    print('building {}'.format(title))
    return TH2F(title,title,200, 0.0, 6.0, 200, 0.0, 1.15)
def buildWT(title):
    print('building {}'.format(title))
    return TH2F(title,title,200, 0.0, 4, 200, 0.0, 6.0)

builders={
    'mass': buildMass,
    'mass_small_bins': buildMassSmallBins,
    'asy': buildAsy, # 10 trento bins -> nominal bins for thesis m2
    'asy_mod2': buildAsyMod2, ## 12 trento bins
    'asy_mod3': buildAsyMod3, ## same as 'asy' but with two middle bins merged
    'asy_mod4': buildAsyMod4, ## same as 'asy' but with two middle bins merged
    'q2': buildQ2,
    'xb': buildXb,
    't': buildT,
    'W': buildW,
    'xbt' : buildXbT,
    'wt': buildWT
}


# new mass bins, larger sizing
field_setting='inbNoutb'
bin_var = '11'
mass_bins = readInBinningInfo('bin_mass_rangesVar'+bin_var+'_'+field_setting+'.txt')

overlap_mass_bins = readInOverlapBinningInfo('bin_mass_slidingVar4_'+field_setting+'.txt')

q2_bins=readInBinningInfo('bin_limits_q2_'+field_setting+'.txt')
xb_bins=readInBinningInfo('bin_limits_xb_'+field_setting+'.txt')
t_bins=readInBinningInfo('bin_limits_t_'+field_setting+'.txt')
w_bins=readInBinningInfo('bin_limits_w_'+field_setting+'.txt')

q2_centers = readInBinningInfo('bin_centers_q2_'+field_setting+'.txt')
xb_centers = readInBinningInfo('bin_centers_xb_'+field_setting+'.txt')
t_centers = readInBinningInfo('bin_centers_t_'+field_setting+'.txt')
w_centers = readInBinningInfo('bin_centers_w_'+field_setting+'.txt')


### define the binning scheme for kpkm mass for bins in trento
h_phi_binning = TH1F("h_phi_binning","h_phi_binning",10, -180, 180) # nominal 
bin_low_trento_c3=[-180, -144, -108, -72, -36, 36, 72, 108, 144, 180]
h_phi_binning_config3 = TH1F("h_phi_binning_config3","h_phi_binning_config3", len(bin_low_trento_c3)-1, array('d',bin_low_trento_c3))

h_phi_binning_config2 = TH1F("h_phi_binning_config2","h_phi_binning_config2", 12, -180, 180) # nominal -> systematic check 
bin_low_trento_c4=[-180, -150, -120, -90, -60, -30, 30, 60, 90, 120, 150, 180]
h_phi_binning_config4 = TH1F("h_phi_binning_config4","h_phi_binning_config4", len(bin_low_trento_c4)-1, array('d',bin_low_trento_c4))

equal_q2_bins = []
equal_xb_bins = []
equal_t_bins = []
equal_w_bins = []
## after splitting xb into equal sized bins, look at q2
equal_q2_xb_bins=[]
# list for masses
im_kpkm_bins = []

ff = TFile(sys.argv[1])
ccounter=0

min_mass_kpkm = mass_bins[1]
max_mass_kpkm = mass_bins[2]

print(' min {} - max {} mass ranges'.format(min_mass_kpkm, max_mass_kpkm))

for ev in ff.clas12:
    ccounter+=1
        
    kpkm_mass = ev.phimass
    hel = ev.helicity
    trento  = ev.trento
    q2 = ev.q2
    t = ev.t
    xb = ev.xb
    w =  ev.w

    hel_stat = 0
    if hel<0: hel_stat='neg' 
    elif hel>0: hel_stat='pos'
    
    if hel == 0 : continue
    
    mass_bin = -1
    for bb in range(1,len(mass_bins)):
        if kpkm_mass > mass_bins[bb-1] and kpkm_mass < mass_bins[bb]:
            mass_bin=bb

    sliding_mass_bin = []
    for bb in range(0,len(overlap_mass_bins)):
        if kpkm_mass > overlap_mass_bins[bb][0] and kpkm_mass < overlap_mass_bins[bb][1]:
            sliding_mass_bin.append(bb)

    q2_bin = -1
    for bb in range(1, len(q2_bins)):
        if q2 > q2_bins[bb-1] and q2 < q2_bins[bb]:
            q2_bin=bb

    xb_bin = -1
    for bb in range(1, len(xb_bins)):
        if xb > xb_bins[bb-1] and xb <xb_bins[bb]:
            xb_bin=bb
    
    t_bin = -1
    for bb in range(1, len(t_bins)):
        if t > t_bins[bb-1] and t < t_bins[bb]:
            t_bin=bb

    w_bin = -1
    for bb in range(1, len(w_bins)):
        if w > w_bins[bb-1] and w < w_bins[bb]:
            w_bin=bb
        
  
    # plot entire mass spectrum
    compute_if_absent(histos,'mass',builders['mass']).Fill( kpkm_mass )
    compute_if_absent(histos,'mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
    # plot mass spectrum for each mass bin for plotting purpose
    if mass_bin >= 0: compute_if_absent(histos,'mass_bin'+str(mass_bin),builders['mass']).Fill( kpkm_mass )
    if mass_bin >= 0: compute_if_absent(histos,'mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
    # plot mass spectrum but for the sliding mass bins
    for sb in sliding_mass_bin:
        if sb >= 0: compute_if_absent(histos,'sliding_mass_bin'+str(sb),builders['mass']).Fill( kpkm_mass )
        if sb >= 0: compute_if_absent(histos,'sliding_mass_small_bin'+str(sb),builders['mass_small_bins']).Fill( kpkm_mass )    

    # plot all trento phi spectrum
    compute_if_absent(histos,'trento',builders['asy']).Fill( trento )

    # plot trento in each mass bin for each helicity state
    compute_if_absent(histos,'trento_'+hel_stat,builders['asy']).Fill( trento )
    if( mass_bin >= 0 ): 
        compute_if_absent(histos,'trento_'+hel_stat+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )           ## 10 trento bins 
        compute_if_absent(histos,'mod2_trento_'+hel_stat+'_mass_bin'+str(mass_bin),builders['asy_mod2']).Fill( trento ) ## 12 trento bins   
        compute_if_absent(histos,'mod3_trento_'+hel_stat+'_mass_bin'+str(mass_bin),builders['asy_mod3']).Fill( trento ) ## merged middle bins 8 trento b
        compute_if_absent(histos,'mod4_trento_'+hel_stat+'_mass_bin'+str(mass_bin),builders['asy_mod4']).Fill( trento ) ## merged middle bins 10 trento b
    # plot trento in each of the sliding mass bins for each helicity state
    for sb in sliding_mass_bin:
        if( sb >= 0 ): compute_if_absent(histos,'trento_'+hel_stat+'_sliding_mass_bin'+str(sb),builders['asy']).Fill( trento )

    # plot -t raw, then in bins of q2 and xb 
    compute_if_absent(histos,'t',builders['t']).Fill( t )    
    compute_if_absent(histos,'t_'+hel_stat,builders['t']).Fill( t )    
    if( kpkm_mass > min_mass_kpkm and kpkm_mass < max_mass_kpkm and q2_bin >= 0 and xb_bin >= 0 ):
        compute_if_absent(histos,'t_q2b'+str(q2_bin)+'_xbb'+str(xb_bin),builders['t']).Fill( t )    

    # plot the W spectrum here, so we can look into binning it
    # also check xb vs t for different Q2 vales
    if kpkm_mass > min_mass_kpkm and kpkm_mass < max_mass_kpkm:
        # plot the final q2,-t, xb, w dist here as well. These are used 
        # in the final thesis work
        compute_if_absent(histos,'q2_final',builders['q2']).Fill( q2 )    
        compute_if_absent(histos,'t_final',builders['t']).Fill( t )    
        compute_if_absent(histos,'xb_final',builders['xb']).Fill( xb )    
        compute_if_absent(histos,'w_final',builders['W']).Fill( w )    

        compute_if_absent(histos,'wt_final',builders['wt']).Fill( w, t )    
        if q2 < 2.5:
            compute_if_absent(histos,'xbt_lowQ2_final',builders['xbt']).Fill( t, xb )    
            compute_if_absent(histos,'wt_lowQ2_final',builders['wt']).Fill( w, t )    
        else:
            compute_if_absent(histos,'xbt_highQ2_final',builders['xbt']).Fill( t, xb )
            compute_if_absent(histos,'wt_highQ2_final',builders['wt']).Fill( w, t )        
        if w < 2.1:
            compute_if_absent(histos,'xbt_wbin1',builders['xbt']).Fill( t, xb )    
        elif w > 2.1 and w < 2.5:
            compute_if_absent(histos,'xbt_wbin2',builders['xbt']).Fill( t, xb )    
        elif w > 2.5 and w < 3.1:
            compute_if_absent(histos,'xbt_wbin3',builders['xbt']).Fill( t, xb ) 
            
    
    
    # get the trento phi distribution per mass bins in q2, -t, and xb
    if( t_bin >= 0 and mass_bin >= 0 ): 
        compute_if_absent(histos,'t_bin'+str(t_bin)+'_mass_bin'+str(mass_bin),builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos,'t_bin'+str(t_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'trento_'+hel_stat+'_t_bin'+str(t_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )
    if( q2_bin >= 0 and mass_bin >= 0 ): 
        compute_if_absent(histos,'q2_bin'+str(q2_bin)+'_mass_bin'+str(mass_bin),builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos,'q2_bin'+str(q2_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'trento_'+hel_stat+'_q2_bin'+str(q2_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )
    if( xb_bin >= 0 and mass_bin >= 0 ): 
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_mass_bin'+str(mass_bin),builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'trento_'+hel_stat+'_xb_bin'+str(xb_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )
    if( w_bin >= 0 and mass_bin >= 0 ): 
        compute_if_absent(histos,'w_bin'+str(w_bin)+'_mass_bin'+str(mass_bin),builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos,'w_bin'+str(w_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'trento_'+hel_stat+'_w_bin'+str(w_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )        
        if q2<2.75:
            compute_if_absent(histos,'w_bin'+str(w_bin)+'_mass_bin'+str(mass_bin)+'_lowq2',builders['mass']).Fill( kpkm_mass )
            compute_if_absent(histos,'w_bin'+str(w_bin)+'_mass_small_bin'+str(mass_bin)+'_lowq2',builders['mass_small_bins']).Fill( kpkm_mass )
            compute_if_absent(histos,'trento_'+hel_stat+'_w_bin'+str(w_bin)+'_mass_bin'+str(mass_bin)+'_lowq2',builders['asy']).Fill( trento )

        
    if( kpkm_mass > min_mass_kpkm and kpkm_mass < max_mass_kpkm and q2_bin >= 0 ) : compute_if_absent(histos,'t_q2b'+str(q2_bin),builders['t']).Fill( t )    # integrated over xb
    if( kpkm_mass > min_mass_kpkm and kpkm_mass < max_mass_kpkm and xb_bin >= 0 ) : compute_if_absent(histos,'t_xbb'+str(xb_bin),builders['t']).Fill( t )    # integrated over q2

    ## we want to do limited multiD analysis 
    # do it for low q2 and low xb bin, and high q2 and high xb bin
    if( q2_bin == 2 and xb_bin == 2 ):
        compute_if_absent(histos,'trento_'+hel_stat+'_xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
    if( q2_bin == 1 and xb_bin == 1 ):
        compute_if_absent(histos,'trento_'+hel_stat+'_xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_bin'+str(mass_bin),builders['asy']).Fill( trento )
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_small_bin'+str(mass_bin),builders['mass_small_bins']).Fill( kpkm_mass )
        compute_if_absent(histos,'xb_bin'+str(xb_bin)+'_q2_bin'+str(q2_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
        

    ## get events within signal
    ## original binning was 1.01 to 1.02866 
    ## poster presentation at 1.01 to 1.03
    ## binvar8 was 1.0126 to 1.0268
    if( kpkm_mass > min_mass_kpkm and kpkm_mass < max_mass_kpkm ): #1.05 ): changed to be consistent with binning choice
        equal_q2_bins.append(q2)
        equal_xb_bins.append(xb)
        equal_t_bins.append(t)
        equal_w_bins.append(w)
    if kpkm_mass > max_mass_kpkm:
        im_kpkm_bins.append(kpkm_mass)

    ### Histograms and procedure for method 2 -FX method where the phi mass in fitted in bins of phi_trento
    ## this method is aimed at extracting bsa while subtracting background in signal region
    phi_trento_bin = h_phi_binning.FindBin(trento)
    phi_trento_bin_config2 = h_phi_binning_config2.FindBin(trento)
    phi_trento_bin_config3 = h_phi_binning_config3.FindBin(trento)
    phi_trento_bin_config4 = h_phi_binning_config4.FindBin(trento)

    compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin)+'_'+hel_stat, builders['mass']).Fill( kpkm_mass )
    #check what happens with 12 bins
    compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin_config2)+'_config2_'+hel_stat, builders['mass']).Fill( kpkm_mass )
    ## check phi mass binned using trento config 3, 10 bins but middle two merged so 9 bins
    compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin_config3)+'_config3_'+hel_stat, builders['mass']).Fill( kpkm_mass )
    ## check phi mass binned using trento config 4, 12 bins but middle two merged, so 11 bins
    compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin_config4)+'_config4_'+hel_stat, builders['mass']).Fill( kpkm_mass )

    ## similarly get the same distribution for the x, q2, -t bins integrated over one another
    if( xb_bin>= 0 ):
        compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin)+'_'+hel_stat+'_xb_bin'+str(xb_bin), builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos, 'xb_bin_'+str(xb_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
    if( q2_bin>= 0 ):
        compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin)+'_'+hel_stat+'_q2_bin'+str(q2_bin), builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos, 'q2_bin_'+str(q2_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
    if( t_bin>= 0 ):
        compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin)+'_'+hel_stat+'_t_bin'+str(t_bin), builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos, 't_bin_'+str(t_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
    if( w_bin >= 0):
        compute_if_absent(histos, 'trento_bin'+str(phi_trento_bin)+'_'+hel_stat+'_w_bin'+str(w_bin), builders['mass']).Fill( kpkm_mass )
        compute_if_absent(histos, 'w_bin_'+str(w_bin)+'_mass_small_bin',builders['mass_small_bins']).Fill( kpkm_mass )
        if q2 < 2.75:
            compute_if_absent(histos, 'w_bin_'+str(w_bin)+'_mass_small_bin_lowq2',builders['mass_small_bins']).Fill( kpkm_mass )


fout = TFile('mon_asy_pidtype4a_rmvres12_binningV'+bin_var+'_cutvar0_-'+field_setting+'_v2.root','recreate')
for _, hist in histos.items():
    hist.Write()



ff.Close()
fout.Close()


print( len(equal_q2_bins))
### get binning description 
def getEqualBinInfo( data_list, n_bins, title ):
    print('Producing Binning Information for {}'.format(title))
    data_list.sort()
    final_bins = np.array_split(np.array(data_list), n_bins)
    print('--> Number of entries in first bin {}'.format(len(final_bins[0])))
    print('--> Number of entries in second bin {}'.format(len(final_bins[1])))
    cc=0
    code_out = []
    bin_center = []
    for bb in final_bins:
        #print(bb)
        arithemetic_bin_avg = sum(bb)/len(bb)
        # dp is datapoint
        first_dp = bb[0]
        last_dp = bb[-1]
        print('BIN {}: '.format(cc))
        print('---> Arithmetic mean: {}'.format(arithemetic_bin_avg))
        print('---> Data Range: {0} - {1} '.format(first_dp, last_dp))
        code_out.append(first_dp)
        code_out.append(last_dp)
        bin_center.append(arithemetic_bin_avg)

        cc+=1
    bin_limits=[]
    print(code_out)
    bin_limits.append(code_out[0])
    for ii in range(0, len(code_out)):
        if ii%2 != 0:
            print(' bin {0} edge is {1}'.format(ii, code_out[ii]))
            bin_limits.append(code_out[ii])
            
    print('Python Code Limits: {0}_bins={1}'.format(title,bin_limits))
    print('Python Code Bin Center: {0}_bin_centers={1}'.format(title,bin_center))

    print(' Saving to Text Files ' )
    f_bin_lim_out = open('bin_limits_'+title+'.txt','w')
    f_bin_cntr_out = open('bin_centers_'+title+'.txt','w')
    for ll in bin_limits:
        f_bin_lim_out.write('{:.4f} \n'.format(ll))#str(ll)+' \n')
    for ll in bin_center:
        f_bin_cntr_out.write('{:.4f} \n'.format(ll))#str(ll)+' \n')

    f_bin_lim_out.close()
    f_bin_cntr_out.close()
             

### equal
### get bin info here, then go back and rerun the analyis with the new binning from these functions 
#getEqualBinInfo(equal_q2_bins, 2, 'q2_'+field_setting+'_pidtype3' )
#getEqualBinInfo(equal_xb_bins, 2, 'xb_'+field_setting+'_pidtype3' )
#getEqualBinInfo(equal_t_bins,  2, 't_'+field_setting+'_pidtype3' )
#getEqualBinInfo(equal_w_bins,  2, 'w_'+field_setting+'_pidtype3' )
### check binning of the mass dist

print('-------> split the mass spectrum into sublists equal in size to the signal region')
n_ev_in_phi=len(equal_q2_bins)
print('--> number of events in phi peak to split for chunks {}'.format(n_ev_in_phi))
print(len(im_kpkm_bins))
im_kpkm_bins.sort()
im_kpkm_chunks = [im_kpkm_bins[x:x+n_ev_in_phi] for x in range(0, len(im_kpkm_bins), n_ev_in_phi)]
print('--> first chunk after phi peak ranges from')
print(len(im_kpkm_chunks[0]))
print(' min kpkm mass {}'.format(im_kpkm_chunks[0][0]))
print(' max kpkm mass {}'.format(im_kpkm_chunks[0][-1]))


'''
print('------> split the mass spectrum of equal number of entries over sliding window')
# like [1, 2, 3, 4, 5, 6] -> [1,2,3] [2, 3, 4] [3, 4, 5] [4, 5, 6]
print(' first value in kpkm mass list {}'.format(im_kpkm_bins[0]))
mm_slider=[]
mm=0
special_binning = open('bin_mass_slidingVar4_'+field_setting+'.txt','w')
special_binning.write('{} {} \n'.format(0.985, 1.0107))
special_binning.write('{} {} \n'.format(1.0107, 1.0287))
while mm < len(im_kpkm_bins):
    print(' mm start {} end {}'.format(mm, mm + n_ev_in_phi))
    temp_mm_slider = im_kpkm_bins[mm : mm + n_ev_in_phi]
    print(' window range {} - {}'.format(temp_mm_slider[0], temp_mm_slider[-1]))
    mm_slider.append([temp_mm_slider[0],temp_mm_slider[-1]])
    mm += n_ev_in_phi/2
    print(' advance mm to {}'.format(mm))
    special_binning.write('{0:.4f} {1:.4f} \n'.format(temp_mm_slider[0], temp_mm_slider[-1]))
special_binning.close()
'''
    
