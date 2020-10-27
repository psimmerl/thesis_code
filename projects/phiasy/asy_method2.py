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

#vphimass_neg_helicity_trentobin_9_pass_all_exclusivity


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

