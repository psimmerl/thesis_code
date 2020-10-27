from ROOT import TFile, TH1F, TH2F, TF1, TLine, TLegend
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

gStyle.SetOptStat(0000)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()    
    hhs[obj.GetName()] = obj

hhs_mc = {}
ffmc = TFile(sys.argv[2])
for kk in ffmc.GetListOfKeys():
    obj = kk.ReadObj()
    hhs_mc[obj.GetName()] = obj

tag=sys.argv[3]
mon_out_file_name='monitor_data_sim_final_excl_cuts'+tag+'.pdf'
can = TCanvas('can','can',2000,1000)
can.Print('{}['.format(mon_out_file_name))
can.Divide(4,2)


can.cd(1)
hhs['epkpkmXe_raw'].SetTitle('epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].SetLineColor(kBlack)
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs['epkpkmXe_fail_at_least_c123'].SetLineColor(kBlue)
hhs['epkpkmXe_raw'].Draw()
hhs['epkpkmXe_pass_c123'].Draw('same')

## test number of bins 
nbinsa = hhs['epkpkmXe_raw'].GetXaxis().GetNbins()
nbinsb = hhs['epkpkmXe_pass_c123'].GetXaxis().GetNbins()
nbinsc = hhs['epkpkmXe_fail_at_least_one_exclusivity'].GetXaxis().GetNbins()
print(' number of bins raw %d, number of bins passsing %d, number of bins fail %d ' %(nbinsa, nbinsb, nbinsc) )

hhs['epkpkmXe_fail_at_least_c123'].Draw('same')
can.cd(2)
hhs['ekpkmX_raw'].SetTitle('eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].SetLineColor(kBlack)
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs['ekpkmX_fail_at_least_c013'].SetLineColor(kBlue)
hhs['ekpkmX_raw'].Draw()
hhs['ekpkmX_pass_c013'].Draw('same')
hhs['ekpkmX_fail_at_least_c013'].Draw('same')
hhs['ekpkmX_raw'].SetMinimum(0)
can.cd(3)
hhs['epkpX_raw'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].SetLineColor(kBlack)
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs['epkpX_fail_at_least_c023'].SetLineColor(kBlue)
hhs['epkpX_raw'].Draw()
hhs['epkpX_pass_c023'].Draw('same')
hhs['epkpX_fail_at_least_c023'].Draw('same')
hhs['epkpX_fail_at_least_c023'].SetMinimum(0)
can.cd(4)
hhs['epkmX_raw'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].SetLineColor(kBlack)
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs['epkmX_fail_at_least_c012'].SetLineColor(kBlue)
hhs['epkmX_raw'].Draw()
hhs['epkmX_pass_c012'].Draw('same')
hhs['epkmX_fail_at_least_c012'].Draw('same')
hhs['epkmX_fail_at_least_c012'].SetMinimum(0)

can.cd(5)
hhs_mc['epkpkmXe_raw'].SetTitle('SIM. epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs_mc['epkpkmXe_raw'].SetLineColor(kBlack)
hhs_mc['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs_mc['epkpkmXe_fail_at_least_c123'].SetLineColor(kBlue)
hhs_mc['epkpkmXe_raw'].SetAxisRange(-0.2, 0.5,"X");
hhs_mc['epkpkmXe_raw'].Draw()
hhs_mc['epkpkmXe_pass_c123'].Draw('same')

## test number of bins 
nbinsa = hhs['epkpkmXe_raw'].GetXaxis().GetNbins()
nbinsb = hhs['epkpkmXe_pass_c123'].GetXaxis().GetNbins()
nbinsc = hhs['epkpkmXe_fail_at_least_one_exclusivity'].GetXaxis().GetNbins()
print(' number of bins raw %d, number of bins passsing %d, number of bins fail %d ' %(nbinsa, nbinsb, nbinsc) )

hhs_mc['epkpkmXe_fail_at_least_c123'].Draw('same')
can.cd(6)
hhs_mc['ekpkmX_raw'].SetTitle('SIM. eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs_mc['ekpkmX_raw'].SetLineColor(kBlack)
hhs_mc['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs_mc['ekpkmX_fail_at_least_c013'].SetLineColor(kBlue)
hhs_mc['ekpkmX_raw'].SetAxisRange(0.5, 1.2,"X");
hhs_mc['ekpkmX_raw'].Draw()
hhs_mc['ekpkmX_pass_c013'].Draw('same')
hhs_mc['ekpkmX_fail_at_least_c013'].Draw('same')
hhs_mc['ekpkmX_raw'].SetMinimum(0)
can.cd(7)
hhs_mc['epkpX_raw'].SetTitle('SIM. epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs_mc['epkpX_raw'].SetLineColor(kBlack)
hhs_mc['epkpX_pass_c023'].SetLineColor(kRed)
hhs_mc['epkpX_fail_at_least_c023'].SetLineColor(kBlue)
hhs_mc['epkpX_raw'].SetAxisRange(-0.25, 0.5,"X");
hhs_mc['epkpX_raw'].Draw()
hhs_mc['epkpX_pass_c023'].Draw('same')
hhs_mc['epkpX_fail_at_least_c023'].Draw('same')
hhs_mc['epkpX_fail_at_least_c023'].SetMinimum(0)
can.cd(8)
hhs_mc['epkmX_raw'].SetTitle('SIM. epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs_mc['epkmX_raw'].SetLineColor(kBlack)
hhs_mc['epkmX_pass_c012'].SetLineColor(kRed)
hhs_mc['epkmX_fail_at_least_c012'].SetLineColor(kBlue)
hhs_mc['epkmX_raw'].SetAxisRange(-0.25, 0.5,"X");
hhs_mc['epkmX_raw'].Draw()
hhs_mc['epkmX_pass_c012'].Draw('same')
hhs_mc['epkmX_fail_at_least_c012'].Draw('same')
hhs_mc['epkmX_fail_at_least_c012'].SetMinimum(0)

can.Print(mon_out_file_name)
can.Clear()
can.SetCanvasSize(1200,1200)
can.Divide(2,2)

can.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs_mc['epkpkmXe_pass_c123'].SetLineColor(kBlue)
hhs_mc['epkpkmXe_pass_c123'].Scale(hhs['epkpkmXe_pass_c123'].GetMaximum()/hhs_mc['epkpkmXe_pass_c123'].GetMaximum())
hhs['epkpkmXe_pass_c123'].SetAxisRange(-0.5, 1.0,"X");
hhs['epkpkmXe_pass_c123'].Draw()
hhs_mc['epkpkmXe_pass_c123'].Draw('same+hist')
l1=TLegend(0.8,0.8,0.9,0.9)
l1.AddEntry(hhs['epkpkmXe_pass_c123'],'DATA')
l1.AddEntry(hhs_mc['epkpkmXe_pass_c123'],'SIM')
l1.Draw('same')

can.cd(2)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs_mc['ekpkmX_pass_c013'].SetLineColor(kBlue)
hhs_mc['ekpkmX_pass_c013'].Scale(hhs['ekpkmX_pass_c023'].GetMaximum()/hhs_mc['ekpkmX_pass_c023'].GetMaximum())
hhs['ekpkmX_pass_c013'].SetAxisRange(0.5, 1.5,"X");
hhs['ekpkmX_pass_c013'].Draw()
hhs_mc['ekpkmX_pass_c013'].Draw('same+hist')
l2=TLegend(0.1,0.8,0.2,0.9)
l2.AddEntry(hhs['ekpkmX_pass_c013'],'DATA')
l2.AddEntry(hhs_mc['ekpkmX_pass_c013'],'SIM')
l2.Draw('same')

can.cd(3)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs_mc['epkpX_pass_c023'].SetLineColor(kBlue)
hhs_mc['epkpX_pass_c023'].Scale(hhs['epkpX_pass_c023'].GetMaximum()/hhs_mc['epkpX_pass_c023'].GetMaximum())
hhs['epkpX_pass_c023'].SetAxisRange(-0.5, 1.0,"X");
hhs['epkpX_pass_c023'].Draw()
hhs_mc['epkpX_pass_c023'].Draw('same+hist')
l3=TLegend(0.8,0.8,0.9,0.9)
l3.AddEntry(hhs['epkpX_pass_c023'],'DATA')
l3.AddEntry(hhs_mc['epkpX_pass_c023'],'SIM')
l3.Draw('same')

can.cd(4)
hhs['epkmX_pass_c012'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs_mc['epkmX_pass_c012'].SetLineColor(kBlue)
hhs_mc['epkmX_pass_c012'].Scale(hhs['epkmX_pass_c012'].GetMaximum()/hhs_mc['epkmX_pass_c012'].GetMaximum())
hhs['epkmX_pass_c012'].SetAxisRange(-0.5, 1.0,"X");
hhs['epkmX_pass_c012'].Draw()
hhs_mc['epkmX_pass_c012'].Draw('same+hist')
l4=TLegend(0.8,0.8,0.9,0.9)
l4.AddEntry(hhs['epkmX_pass_c012'],'DATA')
l4.AddEntry(hhs_mc['epkmX_pass_c012'],'SIM')
l4.Draw('same')

can.Print(mon_out_file_name)
can.Clear()


can.Print('{}]'.format(mon_out_file_name))
