from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import kRed
from array import array
import math
import sys

fin = TFile(sys.argv[1])

h_pos_asym = fin.Get("pos_asy_pass_all_exclusivity_vphim")
h_neg_asym = fin.Get("neg_asy_pass_all_exclusivity_vphim")
#h_pos_asym.RebinX(2)
#h_neg_asym.RebinX(2)

#remove lambda from phi sample with hard cut
h_pos_asym_nlbd = fin.Get("pos_asy_pass_all_exclusivity_vphim")
h_neg_asym_nlbd = fin.Get("neg_asy_pass_all_exclusivity_vphim")#remove_lmda_vphim")
#h_pos_asym_nlbd.RebinX(2)
#h_neg_asym_nlbd.RebinX(2)

h_pos_asym_fail = fin.Get("pos_asy_fail_all_exclusivity")
h_neg_asym_fail = fin.Get("neg_asy_fail_all_exclusivity")
#h_pos_asym_fail.RebinX(2)
#h_neg_asym_fail.RebinX(2)

h_p_copy = TH1F(h_pos_asym)
h_n_copy = TH1F(h_neg_asym)

h_p_copy_nlbd = TH1F(h_pos_asym_nlbd)
h_n_copy_nlbd = TH1F(h_neg_asym_nlbd)

h_p_copy_fail = TH1F(h_pos_asym_fail)
h_n_copy_fail = TH1F(h_neg_asym_fail)

h_p_copy2 = TH1F(h_pos_asym)
h_n_copy2 = TH1F(h_neg_asym)

h_p_copy.Add(h_n_copy,-1)
h_p_copy_nlbd.Add(h_n_copy_nlbd,-1)
h_p_copy_fail.Add(h_n_copy_fail,-1)



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

# create fit function for asy
fasy = TF1("fasy","[0]*sin(x*3.1415926/180.0)",0.0,361.0)
fasy.SetParameter(0,1)

c1 = TCanvas("c1", "canvas", 800, 800)
gStyle.SetOptStat(0000)
c1.Divide(2,2)

c1.cd(1)
h_pos_asym_nlbd.SetTitle("Counts with Pos. Helicity; #phi_{trento}; Counts")
h_pos_asym_nlbd.Draw("HIST+E")

c1.cd(2)
h_neg_asym_nlbd.SetTitle("Counts with Neg. Helicity; #phi_{trento}; Counts")
h_neg_asym_nlbd.Draw("HIST+E")

c1.cd(3)
h_p_copy_nlbd.SetTitle("Pos.-Neg. Helicity; #phi_{trento}; Counts")
h_p_copy_nlbd.Draw("HIST+E")

c1.cd(4)
g_asy = getBSA(h_pos_asym_nlbd,h_neg_asym_nlbd) #TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
g_asy.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy.Fit("fasy")
g_asy.Draw("AP")
g_asy.GetHistogram().SetMaximum(1.0)
g_asy.GetHistogram().SetMinimum(-1.0)
g_asy.Draw("AP")
fasy.Draw("same")
c1.SaveAs("phi_asy_detailed_nlbd.pdf")

c1a = TCanvas('c1a','c1a',1000,1000)
g_asy2 = getBSA(h_pos_asym_nlbd,h_neg_asym_nlbd) #TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
g_asy2.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy2.Draw("AP")
g_asy2.GetHistogram().SetMaximum(1.0)
g_asy2.GetHistogram().SetMinimum(-1.0)
g_asy2.Fit("fasy")
g_asy2.Draw("AP")
fasy.Draw("same")
c1a.SaveAs("phi_asy_final.pdf")


############################################3
## check Background asymmetry - everything that failed cuts

c2 = TCanvas("c2", "canvas2", 800, 800)
gStyle.SetOptStat(0000)
c2.Divide(2,2)

c2.cd(1)
h_pos_asym_fail.SetTitle("Counts with Pos. Helicity; #phi_{trento}; Counts")
h_pos_asym_fail.Draw()

c2.cd(2)
h_neg_asym_fail.SetTitle("Counts with Neg. Helicity; #phi_{trento}; Counts")
h_neg_asym_fail.Draw()

c2.cd(3)
h_p_copy_fail.SetTitle("Pos.-Neg. Helicity; #phi_{trento}; Counts")
h_p_copy_fail.Draw()

c2.cd(4)
g_asy_fail = getBSA(h_pos_asym_fail,h_neg_asym_fail)  #TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
g_asy_fail.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy_fail.Fit("fasy")
g_asy_fail.Draw("AP")
fasy.Draw("same")
c2.SaveAs("phi_asy_fail.pdf")

c3 = TCanvas('c3','c3',800,800)
gStyle.SetOptStat(0000)
c3.Divide(2,2)

c3.cd(1)
h_pos_asym_nlbd.SetTitle("Counts with Pos. Helicity; #phi_{trento}; Counts")
h_pos_asym_nlbd.Draw()

c3.cd(2)
h_neg_asym_nlbd.SetTitle("Counts with Neg. Helicity; #phi_{trento}; Counts")
h_neg_asym_nlbd.Draw()

c3.cd(3)
h_p_copy_nlbd.SetTitle("Pos.-Neg. Helicity; #phi_{trento}; Counts")
h_p_copy_nlbd.Draw()

c3.cd(4)
g_asy_nlbd = getBSA(h_pos_asym_nlbd,h_neg_asym_nlbd) #TGraphErrors(len(tren_phi), tren_phi, asy, tren_phi_er, asy_er)
g_asy_nlbd.SetTitle("Measured Asym; #phi_{trento}; Asym")
g_asy_nlbd.Draw("AP")
g_asy_nlbd.GetHistogram().SetMaximum(1.0)
g_asy_nlbd.GetHistogram().SetMinimum(-1.0)
g_asy_nlbd.Fit("fasy")
g_asy_nlbd.Draw("AP")
fasy.Draw("same")
c3.SaveAs("phi_asy_nlbd.pdf")

############################################3
############################################
###########################################
###########################################
phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution

## Phi mass region used for asymmetry measurement
phim =fin.Get('vphimass_pass_all_exclusivity')
c2 = TCanvas("c2","c2",900,500)
c2.cd(1)
phim.SetTitle("Selected Invariant Mass of K^{+}K^{-} for BSA;IM K^{+}K^{-} (GeV);counts")
phim.GetXaxis().CenterTitle()
phim.GetYaxis().CenterTitle()
phim.Draw("PE")
c2.Update()
p1=c2.GetPad(1)
ymax=gPad.GetUymax()
phisel = TBox(min_phi_mass_signal_limit,0,max_phi_mass_signal_limit,ymax)
phisel.SetFillColorAlpha(kRed,0.15)
phisel.Draw('same')
c2.SaveAs("phiPeakCutForAsy.pdf")
