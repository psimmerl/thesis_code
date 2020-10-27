from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TMath, TFormula
from ROOT import TLegend, TLatex
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad, kBlack, kRed, kBlue, kMagenta, kPink
from array import array
import numpy as np
import math
import sys


gStyle.SetOptStat(0000)

fin = TFile(sys.argv[1])

out_name = sys.argv[2]
#pass_all_exclusivity_rmv_both_lambdas
phim = fin.Get('mass_small_bin') #im_kpkm_pass_all_pass_additional')#vphimass_pass_all_exclusivity')#_rmv_both_lambdas')

gaus = "[0]*TMath::Gaus(x,[1],[2])"
hvyside = "(x>[4])*2/TMath::Pi()*TMath::ATan( (x-[4])/[5])"
hvyside2 = "(x>[3])*2/TMath::Pi()*TMath::ATan( (x-[3])/[4])"
bck = "[3]*(-x*[6] + [7])" #TMath::Exp(-TMath::Power(TMath::Abs([6]*(x-[4])),[7]))"
bck = "TMath::Exp(-TMath::Power(TMath::Abs([6]*(x-[4])),[7]))"
form = gaus +"+"+ hvyside + "*" + bck


#A*(IM-IM_th)^B * exp[-C*(IM-IM_th)]
bck2="[5]*TMath::Power((x - 0.9),[6])*TMath::Exp([7]*(x-0.9))"
form2 = gaus + " + " + hvyside2 + " * " + bck2
bckform=hvyside2 + "*" + bck2
#fitparm2 = [200, 1.019, 0.2, 100, 0.987, 2, 2]
fitparm2 = [850, 1.019, 0.015, 0.987, 0.945, 50, 2,2] #[350, 1.02, 0.015, 0.987, 0.0045, 50, 2,2]

#print(" my function to fit phi is " + form)
phif = TF1("fitPhi",form2, 0.98, 1.3) # was 1.2
print(phif)
fitparm = [300.2, 1.02, 0.022, 50, 0.987, 0.0045, 2.19, 3]

for pp in range(len(fitparm2)):
    print("par number %d, value of %f" %(pp, fitparm2[pp]) )
    phif.SetParameter(pp,fitparm2[pp])
    #phif.SetParLimits(pp, 0.1*fitparm2[pp], 2.5*fitparm2[pp])

phif.SetNpx(10000)    
phif.FixParameter(3,fitparm2[3])
phif.FixParameter(4,fitparm2[4])
phim.Fit("fitPhi","BRN")

bck = TF1("bck2",form2,0.98,1.3)
bck.SetNpx(10000)
for pp in range(len(fitparm2)):
    print("par number %d, value of %f" %(pp, fitparm2[pp]) )
    bck.SetParameter(pp,phif.GetParameter(pp))

bck.FixParameter(0,0)

# get the number of events from background and signal from the fit results 
# define region from 0.98 to 1.036 for integral
bin_width = phim.GetXaxis().GetBinWidth(1)
print(' bin width to normalize {}'.format(bin_width))
fsignal=TF1('fsignal',form2, 0.98, 1.2036)
#fsignal.SetNpx(10000)
fsignal.SetParameter(0,phif.GetParameter(0))
fsignal.SetParameter(1,phif.GetParameter(1))
fsignal.SetParameter(2,phif.GetParameter(2))
fsignal.SetParameter(3,phif.GetParameter(3))
fsignal.SetParameter(4,phif.GetParameter(4))
fsignal.SetParameter(5,phif.GetParameter(5))
fsignal.SetParameter(6,phif.GetParameter(6))
fsignal.SetParameter(7,phif.GetParameter(7))

fbck=TF1('fbck2',bckform, 0.98, 1.2036)
fbck.SetParameter(3,phif.GetParameter(3))
fbck.SetParameter(4,phif.GetParameter(4))
fbck.SetParameter(5,phif.GetParameter(5))
fbck.SetParameter(6,phif.GetParameter(6))
fbck.SetParameter(7,phif.GetParameter(7))

lineshape  = "[3]*TMath::Gaus(x,[4],[5]) + (x>[6])*2/TMath::Pi()*TMath::ATan( (x-[6])/[7]) + [8]*TMath::Power((x - 0.9),[9])*TMath::Exp([10]*(x-0.9))"
f_fom = TF1("fom","gaus(0)/TMath::Sqrt("+lineshape+")",1.019, 1.04)
f_fom.SetParameter(0,phif.GetParameter(0))
f_fom.SetParameter(1,phif.GetParameter(1))
f_fom.SetParameter(2,phif.GetParameter(2))
f_fom.SetParameter(3,phif.GetParameter(0))
f_fom.SetParameter(4,phif.GetParameter(1))
f_fom.SetParameter(5,phif.GetParameter(2))
f_fom.SetParameter(6,phif.GetParameter(3))
f_fom.SetParameter(7,phif.GetParameter(4))
f_fom.SetParameter(8,phif.GetParameter(5))
f_fom.SetParameter(9,phif.GetParameter(6))
f_fom.SetParameter(10,phif.GetParameter(7))

c0a = TCanvas("c0a","c0a",900,900)
c0a.Divide(1,3)
c0a.cd(1)

sig_gaus = TF1('sig_gaus','gaus(0)',0.98,1.05)
sig_gaus.SetParameter(0,phif.GetParameter(0))
sig_gaus.SetParameter(1,phif.GetParameter(1))
sig_gaus.SetParameter(2,phif.GetParameter(2))
sig_gaus.SetNpx(10000)

fsignal.SetNpx(10000)
fsignal.Draw()
sig_gaus.SetLineColor(kBlack)
sig_gaus.Draw("same")

c0a.cd(2)
f_fom.Draw()
c0a.cd(3)

mean =  phif.GetParameter(1)
sigma = phif.GetParameter(2)

x_center = array('d')                                                                                                                                                                                  
y_val = array('d')                                                                                                                                                                                       
x_err = array('d')                                                                                                                                                                                       
y_err = array('d') 
xaxis_range = np.arange(0, 7*sigma, 0.0001)
for ii in range(1, len(xaxis_range)-1):
    delta_x = xaxis_range[ii]
    xx_max = mean + delta_x
    xx_min = mean - delta_x
    print(' integral range: {} - {}'.format(xx_min,xx_max))

    #xx_from_peak =  xaxis_range[ii] - sig_gaus.GetParameter(1)
    xx_from_peak =  xaxis_range[ii]
    int_lineshape = fsignal.Integral(xx_min, xx_max)
    int_gaus = sig_gaus.Integral(xx_min, xx_max)
    fom = int_gaus/np.sqrt(int_lineshape)
    x_center.append(xx_from_peak/sigma)
    print('center {0} : int sig {1} , eval bck {2} , fom {3}'.format(xx_from_peak/sigma, int_gaus, int_lineshape,fom)) 
    y_val.append(fom)
    x_err.append(0)
    y_err.append(0)

g_fom = TGraphErrors(len(x_center),x_center,y_val, x_err, y_err)
g_fom.SetTitle("")
g_fom.SetMarkerStyle(20)
g_fom.SetMarkerColor(kBlack)
g_fom.SetMarkerSize(1)
g_fom.GetXaxis().SetLimits(0, 5)
g_fom.GetHistogram().SetMinimum(0)
g_fom.Draw("AL")



c0a.SaveAs("fom.pdf")



'''

### get signal and associated errors here

temp_min = phif.GetParameter(1) - phif.GetParameter(2)
temp_max = phif.GetParameter(1) + 2.5*phif.GetParameter(2)
print(' min range {} max range {}'.format(temp_min,temp_max))
phi_n_signal = fsignal.Integral(temp_min, temp_max)/bin_width;
phi_n_bck = fbck.Integral(temp_min, temp_max)/bin_width;

print(' integral of sig {}'.format(fsignal.Integral(temp_min, temp_max)))
print(' integral of bck {}'.format(fbck.Integral(temp_min, temp_max)))

print(' bin width: {} '.format(bin_width) )
print(' integrated signal value /binwidth: {}'.format(phi_n_signal))
print(' integrated bck value/binwidth: {} '.format(phi_n_bck))


## get the phi signal with N - dN AND S - dS and max with (+) signs
fit_sigma_err = phif.GetParError(2)
fit_norm_err = phif.GetParError(0)
fit_sigma = phif.GetParameter(2)
fit_norm  = phif.GetParameter(0)

max_fsignal = TF1('max_fsignal',form2, 0.98, 1.2036)
max_fsignal.SetParameter(0,phif.GetParameter(0)+phif.GetParError(0)) # N to change
max_fsignal.SetParameter(1,phif.GetParameter(1))
max_fsignal.SetParameter(2,phif.GetParameter(2)+phif.GetParError(2)) # S to change
max_fsignal.SetParameter(3,phif.GetParameter(3))
max_fsignal.SetParameter(4,phif.GetParameter(4))
max_fsignal.SetParameter(5,phif.GetParameter(5))
max_fsignal.SetParameter(6,phif.GetParameter(6))
max_fsignal.SetParameter(7,phif.GetParameter(7))

min_fsignal = TF1('min_fsignal',form2, 0.98, 1.2036)
min_fsignal.SetParameter(0,phif.GetParameter(0)-phif.GetParError(0)) # N to change
min_fsignal.SetParameter(1,phif.GetParameter(1))
min_fsignal.SetParameter(2,phif.GetParameter(2)-phif.GetParError(2)) # S to change
min_fsignal.SetParameter(3,phif.GetParameter(3))
min_fsignal.SetParameter(4,phif.GetParameter(4))
min_fsignal.SetParameter(5,phif.GetParameter(5))
min_fsignal.SetParameter(6,phif.GetParameter(6))
min_fsignal.SetParameter(7,phif.GetParameter(7))

## get integral of max and min functions
max_phi_n_signal = max_fsignal.Integral(temp_min, temp_max)/bin_width;
min_phi_n_signal = min_fsignal.Integral(temp_min, temp_max)/bin_width;

print(' ---- > getting errors for signal here ' )
print(' max_phi_n_signal {}'.format(max_phi_n_signal))
print(' min_phi_n_signal {}'.format(min_phi_n_signal))
print(' background {}'.format(phi_n_bck))
print(' max sig - bck: {}'.format(max_phi_n_signal-phi_n_bck))
print(' min sig - bck: {}'.format(min_phi_n_signal-phi_n_bck))


sig_gaus = TF1('sig_gaus','gaus(0)',0.98,1.05)
sig_gaus.SetParameter(0,phif.GetParameter(0))
sig_gaus.SetParameter(1,phif.GetParameter(1))
sig_gaus.SetParameter(2,phif.GetParameter(2))
sig_gaus.SetLineColor(kMagenta)
c1 = TCanvas("c1","c1",900,500)
c1.cd(1)
phim.SetTitle("")#Invariant Mass of Charged Kaons from epK^{+}K^{-};IM K^{+}K^{-} (GeV);counts")
#phim.GetXaxis().CenterTitle()
#phim.GetYaxis().CenterTitle()
phim.GetXaxis().SetRangeUser(0.95, 1.140)
phim.Draw("PE")
bck.SetLineColor(4)

phif.SetLineWidth(3)
bck.SetLineStyle(2)
phif.Draw("same")
bck.Draw("same")
sig_gaus.Draw('same')
label=TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.05)
label.SetTextColor(1)
label.DrawLatex(0.1, 0.925,'Invariant Mass of Charged Kaons from epK^{+}K^{-}')
label.DrawLatex(0.5, 0.015,'Invariant Mass K^{+}K^{-} (GeV)')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'Counts')
label.SetTextAngle(0)
d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_norm*fit_sigma_err ,2)  ) * phi_n_signal # use specific range, FX integrals -inf to inf.
label.DrawLatex(0.52, 0.82, 'Number of Phi Events: {0} #pm {1}'.format(int(phi_n_signal - phi_n_bck), int(d_nsignal) ) )
label.DrawLatex(0.52, 0.76, '#mu = {0:6.5} GeV #pm {1:2.2f} MeV'.format(phif.GetParameter(1), phif.GetParError(1)*1000.0))
label.DrawLatex(0.52, 0.7,  '#sigma = {0:2.2f} #pm {1:2.2f} MeV'.format(phif.GetParameter(2)*1000.0, phif.GetParError(2)*1000.0))
label.DrawLatex(0.52, 0.64, 'FWHM = {0:2.2f} MeV'.format(2*TMath.Sqrt(2.0*TMath.Log(2.0))*phif.GetParameter(2)*1000.0))
c1.SaveAs("phiPeak_{}.png".format(out_name))
c1.SaveAs("phiPeak_{}.pdf".format(out_name))


############################################
############################################
c2 = TCanvas("c2","c1",900,500)
c2.cd(1)

phim.Draw("PE")
phim.SetMaximum(1650)
min_fsignal.SetLineColor(kBlack)
max_fsignal.SetLineColor(kBlack)
min_fsignal.SetLineWidth(1)
max_fsignal.SetLineWidth(1)
min_fsignal.SetLineStyle(3)
max_fsignal.SetLineStyle(4)

phif.Draw("same")
bck.Draw("same")
sig_gaus.Draw('same')
min_fsignal.Draw("same")
max_fsignal.Draw("same")

label.DrawLatex(0.1, 0.925,'Invariant Mass of Charged Kaons from epK^{+}K^{-}')
label.DrawLatex(0.5, 0.015,'Invariant Mass K^{+}K^{-} (GeV)')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'Counts')
label.SetTextAngle(0)

label.DrawLatex(0.53, 0.82, 'Number of Phi Events: {}'.format(int(phi_n_signal - phi_n_bck)))
label.DrawLatex(0.53, 0.76, 'Max Phi Events-Bck: {}'.format(int(max_phi_n_signal - phi_n_bck)))
label.DrawLatex(0.53, 0.7,  'Min of Phi Events-Bck: {}'.format(int(min_phi_n_signal - phi_n_bck)))
label.DrawLatex(0.53, 0.64, '#Delta(Max-Nom): {}'.format(int(max_phi_n_signal - phi_n_signal)))
label.DrawLatex(0.53, 0.58, '#Delta(Nom-Min): {}'.format(int(phi_n_signal - min_phi_n_signal)))

label.DrawLatex(0.53, 0.5,  'Raw Nom. Phi Events: {}'.format(int(phi_n_signal)))
label.DrawLatex(0.53, 0.44, 'Raw Bck Events: {}'.format(int(phi_n_bck)))

#label.DrawLatex(0.56, 0.76, '#mu = {0:6.5} #pm {1:6.5} GeV'.format(phif.GetParameter(1), phif.GetParError(1)))
#label.DrawLatex(0.56, 0.7,  '#sigma = {0:6.6f} #pm {1:6.6f} GeV'.format(phif.GetParameter(2), phif.GetParError(2)))
#label.DrawLatex(0.56, 0.64, 'FWHM = {0:6.6f} GeV'.format(2*TMath.Sqrt(2.0*TMath.Log(2.0))*phif.GetParameter(2)))
c2.SaveAs("phiPeak_detailed_{}.png".format(out_name))
c2.SaveAs("phiPeak_detailed_{}.pdf".format(out_name))

##############################################################
## plot figure of merit to determine where to define the signal bins
#sig_gaus.Eval( n_sigma )

#fit_signal_to_background = TF1("fit_signal_to_background","sig_gaus/fbck2",sig_gaus.GetParameter(1), sig_gaus.GetParameter(1) + 7*sig_gaus.GetParameter(2))
c3 = TCanvas("C3","C3",600,1200)
c3.Divide(1,2)
c3.cd(1)

x_center_sb = array('d')
y_val_sb = array('d')
x_err_sb = array('d')
y_err_sb = array('d')

for xx in np.arange(sig_gaus.GetParameter(1) - 0.1*sig_gaus.GetParameter(2), sig_gaus.GetParameter(1) +3*sig_gaus.GetParameter(2), 0.0001):
    xx_from_peak = xx - sig_gaus.GetParameter(1)
    eval_signal = sig_gaus.Eval(xx)
    eval_bck = bck.Eval(xx)
    sb_ratio = eval_signal/( eval_bck )
    
    x_center_sb.append(xx_from_peak/sig_gaus.GetParameter(2))
    #print('center {0} : eval sig {1} , eval bck {2} , fom {3}'.format(xx_from_peak/sig_gaus.GetParameter(2), eval_signal, eval_bck, fom)) 
    y_val_sb.append(sb_ratio)
    x_err_sb.append(0)
    y_err_sb.append(0)

g_sb = TGraphErrors(len(x_center_sb),x_center_sb,y_val_sb, x_err_sb, y_err_sb)
g_sb.SetTitle("")
g_sb.SetMarkerStyle(20)
g_sb.SetMarkerColor(kBlack)
g_sb.SetMarkerSize(1)
g_sb.GetXaxis().SetLimits(0, 5)
g_sb.GetHistogram().SetMinimum(0)
g_sb.Draw("AP")

label.DrawLatex(0.1, 0.925,'Signal/Background from {0:.4f}-{1:.4f}'.format(sig_gaus.GetParameter(1) - 0.1*sig_gaus.GetParameter(2), sig_gaus.GetParameter(1) +3*sig_gaus.GetParameter(2)))
label.DrawLatex(0.5, 0.025,'N_{#sigma} from #phi peak')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'Signal/Background')
label.SetTextAngle(0)


c3.cd(2)
x_center = array('d')
y_val = array('d')
x_err = array('d')
y_err = array('d')

xaxis_range = np.arange(sig_gaus.GetParameter(1) - 0.1*sig_gaus.GetParameter(2), sig_gaus.GetParameter(1) +3*sig_gaus.GetParameter(2), 0.0001)
for ii in range(0, len(xaxis_range)-1):
    xx_min = xaxis_range[ii]
    xx_max = xaxis_range[ii+1]

    #xx_from_peak =  xaxis_range[ii] - sig_gaus.GetParameter(1)
    xx_from_peak =  (xx_min+xx_max)/2.0 - sig_gaus.GetParameter(1)
    print('{} - {} '.format(xx_min, xx_max))
    int_signal = sig_gaus.Integral(xx_min, xx_max)/bin_width
    int_lineshape = fsignal.Integral(xx_min, xx_max)/bin_width
    
    fom = int_signal/np.sqrt( int_lineshape )
    x_center.append(xx_from_peak/sig_gaus.GetParameter(2))
    #print('center {0} : eval sig {1} , eval bck {2} , fom {3}'.format(xx_from_peak/sig_gaus.GetParameter(2), eval_signal, eval_bck, fom)) 
    print('center {0} : integral sig {1} , integral bck {2} , fom {3}'.format(xx_from_peak/sig_gaus.GetParameter(2), int_signal, int_lineshape, fom)) 
    y_val.append(fom)
    x_err.append(0)
    y_err.append(0)

g_fom = TGraphErrors(len(x_center),x_center,y_val, x_err, y_err)
g_fom.SetTitle("")
g_fom.SetMarkerStyle(20)
g_fom.SetMarkerColor(kBlack)
g_fom.SetMarkerSize(1)
g_fom.GetXaxis().SetLimits(0, 5)
g_fom.GetHistogram().SetMinimum(0)
g_fom.Draw("AP")

label.DrawLatex(0.1, 0.925,'FOM from {0:.4f}-{1:.4f}'.format(sig_gaus.GetParameter(1) - 0.1*sig_gaus.GetParameter(2), sig_gaus.GetParameter(1) +3*sig_gaus.GetParameter(2)))
label.DrawLatex(0.5, 0.025,'N_{#sigma} from #phi peak')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'FOM')
label.SetTextAngle(0)


c3.SaveAs("phiPeak_fom_{}.pdf".format(out_name))

'''
'''
phiraw = fin.Get('vphimass_raw')
phipall= fin.Get('vphimass_pass_all_exclusivity')
phipalllmbd = fin.Get('vphimass_pass_all_exclusivity_pass_lmdacut')


vari = phif.GetParameter(2)
mean = phif.GetParameter(1)

print( ' using mean of %f and sigma of %f ' % (mean, vari))
for ii in range(1,5):
    print(' --- checking S/B ration at level of %d ---' % (ii))
    max_sig_range = mean + ii*vari
    min_sig_range = mean - ii*vari

    n_signal = fsignal.Integral(min_sig_range, max_sig_range)/bin_width;
    n_bck = fbck.Integral(min_sig_range, max_sig_range)/bin_width;
    print(' min range %f to max range %f ' %( min_sig_range, max_sig_range) )

    print('signal %f ' % (n_signal-n_bck))
    print('bck %f ' % (n_bck))
    ratio = (n_signal-n_bck)/n_bck
    print('S/B RATIO %f ' % (ratio) )




############################################
############################################
#compare phi mass 
c2 = TCanvas('c2','c2',900,900)
c2.cd(1)
phiraw.SetTitle('Invariant Mass of Charged Kaons from ep->epK^{+}K^{-}; I.M. K^{+}K^{-} (GeV);counts')
phiraw.SetLineColor(kBlack)
phipall.SetLineColor(kRed)
#phipalllmbd.SetLineColor(kRed)
phiraw.Draw()
phipall.Draw('same')
l2 = TLegend(0.6,0.5, 0.9, 0.9)
l2.AddEntry(phiraw,'All Events')
l2.AddEntry(phipall,'Final Events')
l2.Draw('same')
#phipalllmbd.Draw('same')
c2.SaveAs('hist_phi_mass_peak_comparison_{}.png'.format(out_name))

'''
