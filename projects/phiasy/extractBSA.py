from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import kRed
from ROOT import Fit
from array import array
import math
import sys


def getBSA(h_pos_asym,h_neg_asym):
    tren_phi = array('d')
    asy = array('d')
    tren_phi_er = array('d')
    asy_er = array('d')

    for i in range(1, h_pos_asym.GetNbinsX()+1):
        print(" bin number %d " % (i) )
        
        P = h_pos_asym.GetBinContent(i)
        N = h_neg_asym.GetBinContent(i)
        bc = h_pos_asym.GetXaxis().GetBinCenter(i)
        
        nom=P-N
        den=P+N
        y=0
        y_err=0
        print(" number of positive values %d, number of negative values %d" % (P, N))        
        if(den==0):
           y=0
           y_err=0
        else:
            y=nom/den
            y_err = math.sqrt( (1 - y**2)/(P+N))

        print(" bin center %d, asy %d, bin center error %d, asy error %d " % (bc, y, 0, y_err))

        tren_phi.append(bc)
        asy.append(y*(1/0.85355))
        tren_phi_er.append(0)
        asy_er.append(y_err)
    g_asy = TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
    return g_asy

def writeBSAToTxt(h_pos_asym,h_neg_asym, f_out):
    for i in range(1, h_pos_asym.GetNbinsX()+1):
        print(" bin number %d " % (i) )
        
        P = h_pos_asym.GetBinContent(i)
        N = h_neg_asym.GetBinContent(i)
        bc = h_pos_asym.GetXaxis().GetBinCenter(i)
        
        nom=P-N
        den=P+N
        y=0
        y_err=0
        print(" number of positive values %d, number of negative values %d" % (P, N))        
        if(den==0):
           y=0
           y_err=0
        else:
            y=nom/den
            y_err = math.sqrt( (1 - y**2)/(P+N))

        print(" bin center %d, asy %d, bin center error %d, asy error %d " % (bc, y, 0, y_err))

        f_out.write(str(bc) + " " + str(y*(1/0.85355)) + " 0 " + str(y_err) + "\n")
    

fin = TFile(sys.argv[1])
signal_sigma = 2.50

f_out = open("bsa_binned_out.txt","w+")

hhs = {}
for kk in fin.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


###########################################################3
## get the measured BSA


h_pos_asym = fin.Get("pos_asy_pass_all_exclusivity_remove_lmbd_vphim")#pos_asy_pass_all_exclusivity_vphim")
h_neg_asym = fin.Get("neg_asy_pass_all_exclusivity_remove_lmda_vphim")#neg_asy_pass_all_exclusivity_vphim")
h_pos_asym.RebinX(2)
h_neg_asym.RebinX(2)

# define the fit for bsa
fasy = TF1("fasy","[0]*sin(x*3.1415926/180.0)",-180,180)
fasy.SetParameter(0,1)

c1a = TCanvas('c1a','c1a',1000,1000)
g_asy2 = getBSA(h_pos_asym,h_neg_asym)
writeBSAToTxt(h_pos_asym,h_neg_asym, f_out)

g_asy2.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy2.Draw("AP")
g_asy2.GetHistogram().SetMaximum(1.0)
g_asy2.GetHistogram().SetMinimum(-1.0)
g_asy2.Fit("fasy")
g_asy2.Draw("AP")
fasy.Draw("same")
c1a.SaveAs("phi_asy_final_"+str(signal_sigma)+"sigma.pdf")


measured_bsa_fit  = fasy.GetParameter(0)
measured_bsa_fit_err = fasy.GetParError(0)
print(' -- the measured bsa is %f +/- %f ' % (measured_bsa_fit, measured_bsa_fit_err) )



#########################################################################
## get bsa due to background 


background_sigma=2*signal_sigma
phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
min_phi_mass_cut = phi_mass_fit - signal_sigma*phi_mass_resolution

phi_mass_signal_limit = phi_mass_fit + signal_sigma*phi_mass_resolution

max_signal = phi_mass_fit + signal_sigma*phi_mass_resolution
min_signal = phi_mass_fit - signal_sigma*phi_mass_resolution


min_mass = phi_mass_signal_limit


asy_fit_par = array('d')
asy_fit_err = array('d')
asy_massrange = array('d')
asy_massrange_error = array('d')


for ii in range(1,20):

    fasy = TF1("fasy","[0]*sin(x*3.1415926/180.0)",-180,180)
    fasy.SetParameter(0,1)


    #pos_asy_pass_all_exclusivity_vphim_massbinrange
    h_neg_binrange = hhs['neg_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii)]
    h_neg_binrange.RebinX(2)
    h_pos_binrange = hhs['pos_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii)]
    h_pos_binrange.RebinX(2)

    g_asy_binrange=getBSA(h_pos_binrange,h_neg_binrange)


    ##corresponding mass side band
    ## Phi mass region used for asymmetry measurement
    phim = hhs['vphimass_pass_all_exclusivity']
    c2 = TCanvas("c2"+str(ii),"c2"+str(ii),900,1800)
    c2.Divide(1,2)
    c2.cd(1)
    phim.SetTitle("Selected Invariant Mass of K^{+}K^{-} - Sideband;IM K^{+}K^{-} (GeV);counts")
    phim.GetXaxis().CenterTitle()
    phim.GetYaxis().CenterTitle()
    phim.Draw("PE")
    c2.Update()
    p1=c2.GetPad(1)
    ymax=gPad.GetUymax()

    
    max_mass = phi_mass_signal_limit + (background_sigma*phi_mass_resolution)*ii
    center_mass_value = (min_mass + max_mass)/2.0

    print(' min mass range %f , max mass range %f ' % (min_mass,max_mass))

    phisel = TBox(min_mass,0,max_mass,ymax)

    min_mass=max_mass

    phisel.SetFillColorAlpha(kRed,0.15)
    phisel.Draw('same')
    c2.cd(2)
    g_asy_binrange.SetTitle('SideBand Asy range '+str(ii))
    g_asy_binrange.Fit('fasy')
    g_asy_binrange.Draw('AP')
    g_asy_binrange.GetHistogram().SetMaximum(1.0)
    g_asy_binrange.GetHistogram().SetMinimum(-1.0)
    g_asy_binrange.Draw('AP')

    fasy.Draw('same')

    asy_fit_par.append(fasy.GetParameter(0))
    asy_fit_err.append(fasy.GetParError(0))
    asy_massrange.append(center_mass_value)
    asy_massrange_error.append(0.0)

    #if( ii % 2 == 0 ):
     #   myc+=1
      #  c2a.cd(myc)

      #c2.SaveAs("phiPeakCutForAsyMassRegion"+str(ii)+"_"+str(signal_sigma)+"sigma.pdf")




#fit_bg_bsa = TF1('fit_bg_bsa','[0]*x + [1]',phi_mass_signal_limit, 1.47)
fit_bg_bsa = TF1('fit_bg_bsa','[0]',phi_mass_signal_limit, 1.47)


c1 = TCanvas('c1','c1',1000,1000)
c1.cd(1)
g_asy_over_massrange = TGraphErrors(len(asy_fit_par),asy_massrange,asy_fit_par,asy_massrange_error,asy_fit_err)
g_asy_over_massrange.SetTitle('Extracted Asymetry For Different Mass Ranges')
g_asy_over_massrange.SetMarkerStyle(8)
g_asy_over_massrange.Draw('AP')
g_asy_over_massrange.GetHistogram().SetMaximum(0.250)
g_asy_over_massrange.GetHistogram().SetMinimum(-0.250)
g_asy_over_massrange.Draw('AP')
g_asy_over_massrange.Fit('fit_bg_bsa')
fit_bg_bsa.Draw('same')
#c1.SaveAs('graph_sidebank_extracted_asy_fit_binrange_'+str(signal_sigma)+'sigma.pdf')

extrapolated_bsa = fit_bg_bsa.GetParameter(0)#*phi_mass_fit + fit_bg_bsa.GetParameter(1)
#extrapolated_slope_err = fit_bg_bsa.GetParError(0)
extrapolated_b_err = fit_bg_bsa.GetParError(0)
extrapolated_bsa_err = extrapolated_b_err# math.sqrt( (extrapolated_slope_err**2 )*(phi_mass_fit**2) + extrapolated_b_err**2 )

print(' ---extrapolating the BSA from background---' )
print(' extrapolated background contribution is %f ' % ( extrapolated_bsa) )


phim = hhs['vphimass_pass_all_exclusivity']
c2 = TCanvas("c2","c2",900,500)
c2.cd(1)
phim.SetTitle("Selected Invariant Mass of K^{+}K^{-} for BSA;IM K^{+}K^{-} (GeV);counts")
phim.GetXaxis().CenterTitle()
phim.GetYaxis().CenterTitle()
phim.Draw("PE")
c2.Update()
p1=c2.GetPad(1)
ymax=gPad.GetUymax()
phisel = TBox(min_signal,0,phi_mass_signal_limit,ymax)
phisel.SetFillColorAlpha(kRed,0.15)
phisel.Draw('same')
#c2.SaveAs("phiPeakCutForAsy_"+str(signal_sigma)+"sigma.pdf")


#####################################################################
# calculate the BSA from measured data

# determine number of events from bg and signal
phim = fin.Get('vphimass_pass_all_exclusivity')
gaus = "[0]*TMath::Gaus(x,[1],[2])"
hvyside = "(x>[3])*2/TMath::Pi()*TMath::ATan( (x-[3])/[4])"
bck="[5]*TMath::Power((x - 1.9),[6])*TMath::Exp([7]*(x-1.9))"

form = gaus + " + " + hvyside + " * " + bck
bckform=hvyside + "*" + bck

#set fit parameters
fitparm = [300, 1.02, 0.015, 0.987, 0.0045, 50, 2,2]

phif = TF1("fitPhi",form, 0.98, 1.2)

for pp in range(len(fitparm)):
    print("par number %d, value of %f" %(pp, fitparm[pp]) )
    phif.SetParameter(pp,fitparm[pp])
    #phif.SetParLimits(pp, 0.1*fitparm2[pp], 2.5*fitparm2[pp])

phif.SetNpx(10000)    
phif.FixParameter(3,fitparm[3])
phim.Fit("fitPhi","BRN")

# get the number of events from background and signal from the fit results 
# define region from 0.98 to 1.036 for integral
bin_width = phim.GetXaxis().GetBinWidth(1)
fsignal=TF1('fsignal',form, 0.98, max_signal)
fsignal.SetParameter(0,phif.GetParameter(0))
fsignal.SetParameter(1,phif.GetParameter(1))
fsignal.SetParameter(2,phif.GetParameter(2))
fsignal.SetParameter(3,phif.GetParameter(3))
fsignal.SetParameter(4,phif.GetParameter(4))
fsignal.SetParameter(5,phif.GetParameter(5))
fsignal.SetParameter(6,phif.GetParameter(6))
fsignal.SetParameter(7,phif.GetParameter(7))

fbck=TF1('fbck2',bckform, 0.98, max_signal)
fbck.SetParameter(3,phif.GetParameter(3))
fbck.SetParameter(4,phif.GetParameter(4))
fbck.SetParameter(5,phif.GetParameter(5))
fbck.SetParameter(6,phif.GetParameter(6))
fbck.SetParameter(7,phif.GetParameter(7))


vari = phif.GetParameter(2)
mean = phif.GetParameter(1)

n_fit = fsignal.Integral(min_signal, max_signal)/bin_width;
n_bck = fbck.Integral(min_signal, max_signal)/bin_width;
print(' min range %f to max range %f ' %( min_signal, max_signal) )

n_signal = n_fit - n_bck
print('signal %f ' % (n_fit-n_bck))
print('bck %f ' % (n_bck))
ratio = n_signal/n_bck
print('S/B RATIO %f ' % (ratio) )


#determine the true bsa
measured_bsa_w = measured_bsa_fit * ( n_signal + n_bck )/n_signal
bck_bsa_w = extrapolated_bsa * (n_bck/n_signal)
true_bsa = measured_bsa_fit * ( n_signal + n_bck )/n_signal - extrapolated_bsa * (n_bck/n_signal)
print(' measured bsa weighted %f and background bsa weighted %f ' % (measured_bsa_w , bck_bsa_w  ) )
print(' true bsa %f ' % (true_bsa) )


#determine the error in the bsa
bsa_m_err = measured_bsa_fit_err
bsa_e_err = extrapolated_bsa_err

print(' bsa measured error %f ' % (bsa_m_err))
print(' bsa extrapolated error %f ' % (bsa_e_err))

sig_err = math.sqrt(n_signal)
bck_err = math.sqrt(n_bck)

total_err_c1 = (bsa_m_err**2) * (( n_signal + n_bck )/n_signal)**2
total_err_c2 = (bsa_e_err**2) * (-n_bck/n_signal)**2
total_err_c3 = ((bck_err**2) * (((measured_bsa_fit - extrapolated_bsa)/n_signal)**2))
total_err_c4 = (sig_err**2) *  ((( n_bck/(n_signal**2)) *  (extrapolated_bsa - measured_bsa_fit))**2)

print(' total errr = sqrt( %f + %f + %f + %f )' % (total_err_c1, total_err_c2, total_err_c3, total_err_c4) )

total_err = math.sqrt((bsa_m_err**2) * (( n_signal + n_bck )/n_signal)**2 + (bsa_e_err**2 * (-n_bck/n_signal)**2) + ((bck_err**2) * ((measured_bsa_fit - extrapolated_bsa)/n_signal)**2) + ((sig_err**2) * ( (n_bck/(n_signal**2))*( extrapolated_bsa - measured_bsa_fit))**2) )
print(' true bsa %f ' % (true_bsa) )
print(' total err %f ' % (total_err) )

#calculate the BSA error

