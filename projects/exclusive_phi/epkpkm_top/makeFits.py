from ROOT import TFile, TH1F, TH2F, TF1, TF1
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012
gStyle.SetOptStat(00000)

#define fitting function
def fitHisto(htemp):
    binmax = htemp.GetMaximumBin()
    histmax = htemp.GetMaximum()

    binlow=binmax
    binhigh=binmax

    while(htemp.GetBinContent(binhigh) >= 0.65*histmax and binhigh<=htemp.GetNbinsX()):
        binhigh+=1
    while(htemp.GetBinContent(binlow) >= 0.65*histmax and binlow>=1):
        binlow-=1

    xlow = htemp.GetBinLowEdge(binlow)
    xhigh= htemp.GetBinLowEdge(binhigh+1)

    htemp.Fit('gaus','','',xlow,xhigh)
    fit_temp = htemp.GetListOfFunctions().FindObject('gaus')
    #fit_temp.SetLineWidth(3)

    return fit_temp

def fitHistoGP(htemp,mass):
    binmax = htemp.GetMaximumBin()
    histmax = htemp.GetMaximum()

    binlow=binmax
    binhigh=binmax

    while(htemp.GetBinContent(binhigh) >= 0.65*histmax and binhigh<=htemp.GetNbinsX()):
        binhigh+=1
    while(htemp.GetBinContent(binlow) >= 0.65*histmax and binlow>=1):
        binlow-=1

    xlow = htemp.GetBinLowEdge(binlow)
    xhigh= htemp.GetBinLowEdge(binhigh+1)
    gpol_fit = TF1('gpol_fit','gaus(0)*pol2(3)',xlow,xhigh)    
    gpol_fit.SetParameter(0,histmax)
    gpol_fit.SetParameter(1,mass)
    gpol_fit.SetParameter(2,0.01)
    gpol_fit.SetParameter(3,1)
    gpol_fit.SetParameter(4,1)
    gpol_fit.SetParameter(5,1)


    htemp.Fit('gpol_fit','R','',xlow,xhigh)
    fit_temp = gpol_fit# htemp.GetListOfFunctions().FindObject('gpol_fit')#aus+pol1')
    #fit_temp.SetLineWidth(3)

    return fit_temp

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


#fit epkmX first pass c012
fit_epkmX_c012 = fitHisto(hhs['LepkmX_pass_c012'])
c0 = TCanvas('cepkmX','cepkmX',900,900)
c0.cd(1)
hhs['LepkmX_pass_c012'].Draw()
hhs['LepkmX_pass_c012'].GetXaxis().SetRangeUser(-0.5,1.5)
c0.Update()
fit_epkmX_c012.SetLineWidth(2)
fit_epkmX_c012.Draw('same')
c0.SaveAs('hist_fit_epkmX_pass_c012.png')

#fit ekpkmX first pass c013
fit_ekpkmX_c013 = fitHisto(hhs['LekpkmX_pass_c013'])
c1 = TCanvas('cekpkmX','cekpkmX',900,900)
c1.cd(1)
hhs['LekpkmX_pass_c013'].Draw()
hhs['LekpkmX_pass_c013'].GetXaxis().SetRangeUser(0,2.0)
c1.Update()
fit_ekpkmX_c013.SetLineWidth(2)
fit_ekpkmX_c013.Draw('same')
c1.SaveAs('hist_fit_ekpkmX_pass_c013.png')


#fit epkpX first pass c023
fit_epkpX_c023 = fitHisto(hhs['LepkpX_pass_c023'])
c2 = TCanvas('cepkpX','cepkpX',900,900)
c2.cd(1)
hhs['LepkpX_pass_c023'].Draw()
hhs['LepkpX_pass_c023'].GetXaxis().SetRangeUser(-0.5,1.0)
c2.Update()
fit_epkpX_c023.SetLineWidth(2)
fit_epkpX_c023.Draw('same')
c2.SaveAs('hist_fit_epkpX_pass_c023.png')

#fit epkpkmXe first pass c123
fit_epkpkmXe_c123 = fitHisto(hhs['LepkpkmXe_pass_c123'])
c3 = TCanvas('cepkpkmXe','cepkpkmXe',900,900)
c3.cd(1)
hhs['LepkpkmXe_pass_c123'].Draw()
hhs['LepkpkmXe_pass_c123'].GetXaxis().SetRangeUser(-0.6,1.2)
c3.Update()
fit_epkpkmXe_c123.SetLineWidth(2)
fit_epkpkmXe_c123.Draw('same')
c3.SaveAs('hist_fit_epkpkmXe_pass_c123.png')


###########################################################################3
#fit epkmX first pass raw

fit_epkmX_raw = fitHistoGP(hhs['LepkmX_raw'],0.207)#0.207)
c0 = TCanvas('cepkmXr','cepkmXr',900,900)
c0.cd(1)
hhs['LepkmX_raw'].SetTitle('epK^{-}X MM2 Raw Spectrum; mm2 (GeV^{2}); counts')
hhs['LepkmX_raw'].Draw()
hhs['LepkmX_raw'].GetXaxis().SetRangeUser(-0.5,1.5)
c0.Update()
fit_epkmX_raw.SetLineWidth(2)
fit_epkmX_raw.Draw('same')
c0.SaveAs('hist_fit_epkmX_raw.png')
print('epkmX raw fit results: mean %f sigma %f' % (fit_epkmX_raw.GetParameter(1), fit_epkmX_raw.GetParameter(2)))


#fit ekpkmX first pass c013
fit_ekpkmX_raw = fitHistoGP(hhs['LekpkmX_raw'],1)
c1 = TCanvas('cekpkmXr','cekpkmXr',900,900)
c1.cd(1)
hhs['LekpkmX_raw'].SetTitle('eK^{+}K^{-}X MM2 Raw Spectrum; mm2 (GeV^{2}); counts')
hhs['LekpkmX_raw'].Draw()
hhs['LekpkmX_raw'].GetXaxis().SetRangeUser(0,2.0)
c1.Update()
fit_ekpkmX_raw.SetLineWidth(2)
fit_ekpkmX_raw.Draw('same')
c1.SaveAs('hist_fit_ekpkmX_raw.png')
print('ekpkmX raw fit results: mean %f sigma %f' % (fit_ekpkmX_raw.GetParameter(1), fit_ekpkmX_raw.GetParameter(2)))


#fit epkpX first pass raw
fit_epkpX_raw = fitHistoGP(hhs['LepkpX_raw'],0.247)
c2 = TCanvas('cepkpXr','cepkpXr',900,900)
c2.cd(1)
hhs['LepkpX_raw'].SetTitle('epK^{+}X MM2 Raw Spectrum; mm2 (GeV^{2}); counts')
hhs['LepkpX_raw'].Draw()
hhs['LepkpX_raw'].GetXaxis().SetRangeUser(-0.5,1.0)
c2.Update()
fit_epkpX_raw.SetLineWidth(2)
fit_epkpX_raw.Draw('same')
c2.SaveAs('hist_fit_epkpX_raw.png')
print('epkpX raw fit results: mean %f sigma %f' % (fit_epkpX_raw.GetParameter(1), fit_epkpX_raw.GetParameter(2)))


#fit epkpkmXe first pass c123
fit_epkpkmXe_raw = fitHistoGP(hhs['LepkpkmXe_raw'],0.01)
c3 = TCanvas('cepkpkmXer','cepkpkmXer',900,900)
c3.cd(1)
hhs['LepkpkmXe_raw'].SetTitle('Missing Energy Raw epK^{+}K^{-}X; missing energy (GeV);counts')
hhs['LepkpkmXe_raw'].Draw()
hhs['LepkpkmXe_raw'].GetXaxis().SetRangeUser(-0.6,1.2)
c3.Update()
fit_epkpkmXe_raw.SetLineWidth(2)
fit_epkpkmXe_raw.Draw('same')
c3.SaveAs('hist_fit_epkpkmXe_raw.png')
print('epkpkmXe raw fit results: mean %f sigma %f' % (fit_epkpkmXe_raw.GetParameter(1), fit_epkpkmXe_raw.GetParameter(2)))

c4 = TCanvas('c4','c4',1000,1000)
c4.Divide(2,2)
c4.cd(1)
hhs['LepkmX_raw'].Draw()
c4.cd(2)
hhs['LekpkmX_raw'].Draw()
c4.cd(3)
hhs['LepkpX_raw'].Draw()
c4.cd(4)
hhs['LepkpkmXe_raw'].Draw()
c4.SaveAs('hist_all_raw_fitted_plots.png')

