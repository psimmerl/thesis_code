from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TMath, TFormula
from ROOT import TLegend, TLatex, TLine
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad, kBlack, kRed, kBlue, kMagenta, kPink

import argparse

from array import array
import numpy as np
import math
import sys


gStyle.SetOptStat(0000)

ap = argparse.ArgumentParser()
ap.add_argument(
    '-i',
    '--input_file',
    required=True
)
ap.add_argument(
    '-o',
    '--output_prefix',
    required=True
)
args = ap.parse_args()

input_rootfile = args.input_file 
out_name = args.output_prefix#sys.argv[2]
fin = TFile(input_rootfile)#sys.argv[1])

field_setting=""
if 'inbNout' in out_name:
    field_setting="inbNoutb"
elif 'inb' in out_name:
    field_setting="inb"
elif 'outb' in out_name:
    field_setting="outb"
else:
    field_setting=""


phim = fin.Get('mass_small_bin')
# use below histo for binnings in xb and q2 and -t
#phim = fin.Get('xb_bin2_q2_bin2_mass_small_bin')
#xb_bin_'+str(xb_bin)+'_mass_small_bin

gaus = "[0]*TMath::Gaus(x,[1],[2])"
bck = "[3]*(x>2*0.4937)*TMath::Power(x-2*0.4937, [4])*TMath::Exp(-[5]*(x - 2*0.4937) )"

form = gaus + " + " + bck
bckform = bck

# in and outbending separately
#fitparm3 = [ 1000 , 1.01937 , 0.00457582 , 2452.97 , 0.339580 , 0.1 ]

#rmv res1
#fitparm3 = [ 1000 , 1.01937 , 0.00427582 , 2452.97 , 0.339580 , 0.1 ]

#rmv res2
#fitparm3 = [ 1000 , 1.01984 , 0.004505 , 2452.97 , 0.339580 , 0.1 ]

#rmv res 1 and 2
fitparm3 = [ 350 , 1.0200 , 0.004207, 2452.97 , 0.6753 , 22 ]


print(" my function to fit phi is " + form)
max_fit_range = 1.065
phif = TF1("fitPhi",form, 2*0.4937, max_fit_range)

for pp in range(0,len(fitparm3)):
    print("par number %d, value of %f" %(pp, fitparm3[pp]) )
    phif.SetParameter(pp,fitparm3[pp])

phif.SetNpx(10000)    
phim.Fit("fitPhi","BRN")


fbck = TF1("bck",form, 2*0.4937, max_fit_range)
for pp in range(3,len(fitparm3)):
    print("par number %d, value of %f" %(pp, phif.GetParameter(pp) ))
    fbck.SetParameter(pp, phif.GetParameter(pp))
    print(' fbck par {} is {}'.format(pp, fbck.GetParameter(pp)))
fbck.SetNpx(10000)
fbck.FixParameter(0,0)


# get the number of events from background and signal from the fit results 
# define region from 0.98 to 1.036 for integral
bin_width = phim.GetXaxis().GetBinWidth(1)
print(' bin width to normalize {}'.format(bin_width))
fsignal=TF1('fsignal',form, 0.98, max_fit_range)
#fsignal.SetNpx(10000)
for ii in range(0, len(fitparm3)):
    fsignal.SetParameter(ii,phif.GetParameter(ii))

### get signal and associated errors here

temp_min = phif.GetParameter(1) - 3.5*phif.GetParameter(2)
temp_max = phif.GetParameter(1) + 3.5*phif.GetParameter(2)
print(' min range {} max range {}'.format(temp_min,temp_max))

phi_n_signal = fsignal.Integral(temp_min, temp_max)/bin_width
phi_n_bck = fbck.Integral(temp_min, temp_max)/bin_width

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

max_fsignal = TF1('max_fsignal',form, 0.98, 1.2036)
for ii in range(0, phif.GetNpar()):
    max_fsignal.SetParameter(ii, phif.GetParameter(ii))
max_fsignal.SetParameter(0,phif.GetParameter(0)+phif.GetParError(0)) # N to change
max_fsignal.SetParameter(2,phif.GetParameter(2)+phif.GetParError(2)) # S to change

min_fsignal = TF1('min_fsignal',form, 0.98, 1.2036)
for ii in range(0, phif.GetNpar()):
    min_fsignal.SetParameter(ii, phif.GetParameter(ii))
min_fsignal.SetParameter(0,phif.GetParameter(0)-phif.GetParError(0)) # N to change
min_fsignal.SetParameter(2,phif.GetParameter(2)-phif.GetParError(2)) # S to change

## get integral of max and min functions
max_phi_n_signal = max_fsignal.Integral(temp_min, temp_max)/bin_width
min_phi_n_signal = min_fsignal.Integral(temp_min, temp_max)/bin_width

print(' ---- > getting errors for signal here ' )
print(' max_phi_n_signal {}'.format(max_phi_n_signal))
print(' min_phi_n_signal {}'.format(min_phi_n_signal))
print(' background {}'.format(phi_n_bck))
print(' max sig - bck: {}'.format(max_phi_n_signal-phi_n_bck))
print(' min sig - bck: {}'.format(min_phi_n_signal-phi_n_bck))

sig_gaus = TF1('sig_gaus','gaus(0)',0.98,max_fit_range)
sig_gaus.SetParameter(0,phif.GetParameter(0))
sig_gaus.SetParameter(1,phif.GetParameter(1))
sig_gaus.SetParameter(2,phif.GetParameter(2))
sig_gaus.SetLineColor(kMagenta)
c1 = TCanvas("c1","c1",900,500)
c1.cd(1)
phim.SetTitle("")
phim.GetXaxis().SetRangeUser(0.95, 1.140)
phim.Draw("PE")
phif.Draw("same")
fbck.SetLineColor(4)
print(fbck)
phif.SetLineWidth(3)
fbck.SetLineStyle(2)
phif.Draw("same")
fbck.Draw("same")
sig_gaus.Draw('same')

# labels 
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

label.SetTextColorAlpha(1, 0.376)
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.08)
label.DrawLatex(0.125, 0.847, field_setting)
label.SetTextSize(0.05)
label.SetTextColorAlpha(1, 1)

normC = np.sqrt( 2.0 * TMath.Pi())/bin_width
d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_norm*fit_sigma_err ,2)  ) * normC # use specific range, FX integrals -inf to inf.
Nsignal = phif.GetParameter(0)*phif.GetParameter(2)*(TMath.Sqrt(2*TMath.Pi())/bin_width)
label.DrawLatex(0.52, 0.82, 'Number of Phi Events: {0} #pm {1}'.format(int(Nsignal), int(d_nsignal) ) )
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
fbck.Draw("same")
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
label.DrawLatex(0.53, 0.52, '#Delta Sig: {0:.4f}'.format(int(d_nsignal)))

label.DrawLatex(0.53, 0.5,  'Raw Nom. Phi Events: {}'.format(int(phi_n_signal)))
label.DrawLatex(0.53, 0.44, 'Raw Bck Events: {}'.format(int(phi_n_bck)))
label.DrawLatex(0.53, 0.39, 'fit N: {}'.format(int(phif.GetParameter(0))))


label.SetTextColorAlpha(1, 0.376)
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.08)
label.DrawLatex(0.125, 0.847, field_setting)
label.SetTextSize(0.05)
label.SetTextColorAlpha(1, 1)
c2.SaveAs("phiPeak_detailed_{}.png".format(out_name))
c2.SaveAs("phiPeak_detailed_{}.pdf".format(out_name))

#######################################################################################################################################
## plot figure of merit to determine where to define the signal bins

mean = sig_gaus.GetParameter(1)
sigma = sig_gaus.GetParameter(2)
xaxis_range = np.arange(0, 7*sigma, 0.0001)

c3 = TCanvas("C3","C3",600,1200)
c3.Divide(1,2)

c3.cd(1)

x_center = array('d')
y_val = array('d')
x_err = array('d')    
y_err = array('d') 
for ii in range(1, len(xaxis_range)-1):
    delta_x = xaxis_range[ii]
    xx_max = mean + delta_x
    xx_min = mean - delta_x
    #print(' integral range: {} - {}'.format(xx_min,xx_max))

    xx_from_peak =  xaxis_range[ii]
    int_lineshape = fsignal.Integral(xx_min, xx_max)
    int_gaus = sig_gaus.Integral(xx_min, xx_max)
    fom = int_gaus/np.sqrt(int_lineshape)
    x_center.append(xx_from_peak/sigma)
    y_val.append(fom)
    x_err.append(0)
    y_err.append(0)

g_fom = TGraphErrors(len(x_center),x_center,y_val, x_err, y_err)
g_fom.SetTitle("")
g_fom.SetMarkerStyle(20)
g_fom.SetMarkerColor(kBlack)
g_fom.SetMarkerSize(1)
g_fom.GetXaxis().SetLimits(0, 7.5)
g_fom.GetHistogram().SetMinimum(0)
g_fom.SetLineWidth(3)
g_fom.Draw("AL")

label.DrawLatex(0.1, 0.925,'FOM from {0:.4f}-{1:.4f} GeV'.format(sig_gaus.GetParameter(1) - 7*sig_gaus.GetParameter(2), sig_gaus.GetParameter(1) + 7*sig_gaus.GetParameter(2)))
label.DrawLatex(0.5, 0.025,'#pm N_{#sigma} from #phi peak')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'FOM')
label.SetTextAngle(0)
label.SetTextColorAlpha(1, 0.376)
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.08)
label.DrawLatex(0.125, 0.847, field_setting)
label.SetTextSize(0.05)
label.SetTextColorAlpha(1, 1)

x_max_of_fom = x_center[y_val.index(max(y_val))]
max_of_fom = max(y_val)
print(' max fom {} at x {}'.format(max_of_fom, x_max_of_fom))
c3.Update()
y_axis_max = gPad.GetUymax()
l_max_of_fom = TLine(x_max_of_fom, 0.0, x_max_of_fom, y_axis_max)
l_max_of_fom.SetLineColor(kRed)
l_max_of_fom.SetLineWidth(3)
l_max_of_fom.Draw("same")
label.DrawLatex(0.33, 0.5,'Max FOM at #pm{0:.4f} #sigma'.format(x_max_of_fom))
label.DrawLatex(0.33, 0.42,'K^{+}K^{-}'+' range: {0:.4f}-{1:.4f} GeV'.format(mean-x_max_of_fom*sigma, mean+x_max_of_fom*sigma))

c3.cd(2)

x_center_sb = array('d')
y_val_sb = array('d')
x_err_sb = array('d')
y_err_sb = array('d')

for ii in range(1, len(xaxis_range)-1):
    delta_x = xaxis_range[ii]
    xx_max = mean + delta_x
    xx_min = mean - delta_x

    xx_from_peak =  xaxis_range[ii]
    int_bck = fbck.Integral(xx_min, xx_max)
    int_gaus = sig_gaus.Integral(xx_min, xx_max)
    sb_ratio = int_gaus/int_bck
    x_center_sb.append(xx_from_peak/sigma)
    #print('center {0} : eval sig {1} , eval bck {2} , fom {3}'.format(xx_from_peak/sigma, int_signal, int_bck, sb_ratio)) 
    y_val_sb.append(sb_ratio)
    x_err_sb.append(0)
    y_err_sb.append(0)

g_sb = TGraphErrors(len(x_center_sb),x_center_sb,y_val_sb, x_err_sb, y_err_sb)
g_sb.SetTitle("")
g_sb.SetMarkerStyle(20)
g_sb.SetMarkerColor(kBlack)
g_sb.SetMarkerSize(1)
g_sb.GetXaxis().SetLimits(0, 7.5)
g_sb.GetHistogram().SetMinimum(0)
g_sb.SetLineWidth(3)
g_sb.Draw("AL")

label.DrawLatex(0.1, 0.925,'Signal/Background from {0:.4f}-{1:.4f} GeV'.format(sig_gaus.GetParameter(1)-7*sig_gaus.GetParameter(2) , sig_gaus.GetParameter(1) + 7*sig_gaus.GetParameter(2)))
label.DrawLatex(0.5, 0.025,'#pm N_{#sigma} from #phi peak')
label.SetTextAngle(90)
label.DrawLatex(0.04, 0.5,'Signal/Background')
label.SetTextAngle(0)
label.SetTextColorAlpha(1, 0.376)
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.08)
label.DrawLatex(0.125, 0.847, field_setting)
label.SetTextSize(0.05)
label.SetTextColorAlpha(1, 1)

c3.Update()
y_axis_max = gPad.GetUymax()
l_max_of_fom_sb = TLine(x_max_of_fom, 0.0, x_max_of_fom, y_axis_max)
l_max_of_fom_sb.SetLineColor(kRed)
l_max_of_fom_sb.SetLineWidth(3)
l_max_of_fom_sb.Draw("same")

sb_ratio_at_max_fom = y_val_sb[y_val.index(max(y_val))]
# determine s/b based on sigma maximizing fom
int_bck = fbck.Integral(mean-x_max_of_fom*sigma, mean+x_max_of_fom*sigma)
int_gaus = sig_gaus.Integral(mean-x_max_of_fom*sigma, mean+x_max_of_fom*sigma)
sb_ratio = int_gaus/int_bck
print('--------------------> range to calculated S/B ratio {} - {}'.format(mean-x_max_of_fom*sigma, mean+x_max_of_fom*sigma))
print('--------------------> directly calculating S/B RATIO gives {}'.format(sb_ratio))
label.DrawLatex(0.35, 0.75,'S/B ratio at Max FOM {0:.4f}'.format(sb_ratio_at_max_fom))

c3.SaveAs("phiPeak_fom_{}.pdf".format(out_name))


###################################################
## calculate the FOM and S/B for the sideband here
xx_side_bands = []
xx_sigma_sb = np.arange(1, 3.5, 0.5)
print(xx_sigma_sb)
sb=0
delta_kpkm_sb = 0.0294
for ss in xx_sigma_sb:
    xx_side_bands.append( [mean - sigma*ss, mean + sigma*ss] )

    int_bck = fbck.Integral(mean - sigma*ss, mean + sigma*ss)
    int_gaus = sig_gaus.Integral(mean - sigma*ss, mean + sigma*ss)
    sb_ratio = int_gaus/int_bck
    print(' for use using SideBand subtraction band number {}'.format(sb))
    print(' sigma {}'.format(ss))
    print(' signal {}'.format(int_gaus))
    print(' bck {}'.format(int_bck))
    print(' signal to bck {}'.format(sb_ratio))
    
    print(' low mass {}'.format(mean - sigma*ss))
    print(' high mass {}'.format(mean + sigma*ss))
    print(mean + sigma*ss + delta_kpkm_sb)
    sb+=1

## special calculation of S/B done here
## this is for method 1 in the thesis
min_range_to_use = 1.0107
max_range_to_use = 1.0287
int_cal_bck = fbck.Integral(min_range_to_use, max_range_to_use)
int_cal_sig = sig_gaus.Integral(min_range_to_use, max_range_to_use)
sb_ratio_to_use = int_cal_sig/int_cal_bck
print('defined signal range from {} - {} has S/B of {}'.format(min_range_to_use, max_range_to_use,sb_ratio_to_use))



