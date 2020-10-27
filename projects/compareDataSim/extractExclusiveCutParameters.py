from ROOT import TFile, TH1F, TH2F, TF1, TLine, TLegend
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012


hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()    
    hhs[obj.GetName()] = obj

tag = sys.argv[2]
mon_out_file_name='extract_excl_cuts_'+tag+'.pdf'
can = TCanvas('can','can',2000,2000)
can.Print('{}['.format(mon_out_file_name))
can.Divide(2,2)

can.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs['epkpkmXe_pass_c123'].Draw()
mean1 = hhs['epkpkmXe_pass_c123'].GetMean()
sig1 =  hhs['epkpkmXe_pass_c123'].GetRMS()
fit1 = TF1('fit1','gaus', mean1-1.5*sig1, mean1+1.5*sig1)
fit1.SetParameter(1,mean1)
fit1.SetParameter(2,sig1)
hhs['epkpkmXe_pass_c123'].Fit(fit1,'R','',mean1-1.5*sig1, mean1+1.5*sig1)
fit1.Draw('same')

can.cd(2)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_pass_c013'].Draw()
mean2 = hhs['ekpkmX_pass_c013'].GetMean()
sig2 =  hhs['ekpkmX_pass_c013'].GetRMS()
fit2 = TF1('fit2','gaus', mean2-2*sig2, mean2+2*sig2)
hhs['ekpkmX_pass_c013'].Fit(fit2,'R','',mean2-2*sig2, mean2+2*sig2)
fit2.Draw('same')


can.cd(3)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_pass_c023'].Draw()
mean3 = hhs['epkpX_pass_c023'].GetMean()
sig3 =  hhs['epkpX_pass_c023'].GetRMS()
fit3 = TF1('fit3','gaus', mean3-2*sig3, mean3+2*sig3)
hhs['epkpX_pass_c023'].Fit(fit3,'R','',mean3-2*sig3, mean3+2*sig3)
fit3.Draw('same')


can.cd(4)
hhs['epkmX_pass_c012'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_pass_c012'].Draw()
mean4 = hhs['epkmX_pass_c012'].GetMean()
sig4 =  hhs['epkmX_pass_c012'].GetRMS()
fit4 = TF1('fit4','gaus', mean4-2*sig4, mean4 + 2*sig4)
hhs['epkmX_pass_c012'].Fit(fit4,'R','',mean4-2*sig4, mean4+2*sig4)
fit4.Draw('same')

can.Print(mon_out_file_name)
can.Print('{}]'.format(mon_out_file_name))



