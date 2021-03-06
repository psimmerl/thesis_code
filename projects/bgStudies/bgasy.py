from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad, kRed
from array import array
import math
import sys

fin = TFile(sys.argv[1])

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj



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
asy_fit_par = array('d')
asy_fit_err = array('d')
asy_massrange = array('d')
asy_massrange_error = array('d')

c2a = TCanvas('c2a','c2a',1800,1200)
c2a.Divide(3,2)
myc=1

phi_mass_fit = 1.021
phi_mass_resolution = 0.564e-03 # 4.25752e-03  
phi_mass_signal_limit = phi_mass_fit + 5*phi_mass_resolution
min_mass = phi_mass_signal_limit

for ii in range(1,23):

    fasy = TF1("fasy","[0]*sin(x*3.1415926/180.0)",0.0,360.0)
    fasy.SetParameter(0,1)


    #pos_asy_pass_all_exclusivity_vphim_massbinrange
    h_neg_binrange = hhs['neg_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii)]
    h_neg_binrange.RebinX(2)
    h_pos_binrange = hhs['pos_asy_pass_all_exclusivity_vphim_massbinrange'+str(ii)]
    h_pos_binrange.RebinX(2)

    g_asy_binrange=getBSA(h_pos_binrange,h_neg_binrange)


    ##corresponding mass side band
    ## Phi mass region used for asymmetry measurement
    phim = hhs['vphimass_pass_all_exclusivity_remove_lambda']
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


    max_mass = phi_mass_signal_limit+ (5*0.00564)*ii
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

    c2.SaveAs("phiPeakCutForAsyMassRegion"+str(ii)+".pdf")




fit_bg_bsa = TF1('fit_bg_bsa','[0]*x + [1]',phi_mass_signal_limit, 1.47)

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
c1.SaveAs('graph_sidebank_extracted_asy_fit_binrange.pdf')

par0 = fit_bg_bsa.GetParameter(0)
par1 = fit_bg_bsa.GetParameter(1)
print('%f %f ' % (par0, par1))



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
phisel = TBox(min_phi_mass_cut,0,phi_mass_signal_limit,ymax)
phisel.SetFillColorAlpha(kRed,0.15)
phisel.Draw('same')
c2.SaveAs("phiPeakCutForAsy.pdf")



    
