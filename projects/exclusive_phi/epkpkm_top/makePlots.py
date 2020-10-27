from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

#gStyle.SetOptStat(1111)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


c0 = TCanvas("c0","c0",1200,450)
c0.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X Energy With All Cuts;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].SetLineColor(kBlack)
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs['epkpkmXe_raw'].Draw()
hhs['epkpkmXe_pass_c123'].Draw('same')
c0.SaveAs('hist_epkpkmXe_cut_plots.png')


c1 = TCanvas("c1","c1",1200,450)
c1.cd(1)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].SetLineColor(kBlack)
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs['ekpkmX_raw'].Draw()
hhs['ekpkmX_pass_c013'].Draw('same')
hhs['ekpkmX_raw'].SetMinimum(0)
c1.SaveAs('hist_ekpkmXe_cut_plots.png')

c2 = TCanvas("c2","c2",1200,450)
c2.cd(1)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].SetLineColor(kBlack)
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs['epkpX_raw'].Draw()
hhs['epkpX_pass_c023'].Draw('same')
hhs['epkpX_pass_c023'].SetMinimum(0)
c2.SaveAs('hist_epkpX_cut_plots.png')


c3 = TCanvas("c3","c3",1200,450)
c3.cd(1)
hhs['epkmX_pass_c012'].SetTitle('epK^{-}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].SetLineColor(kBlack)
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs['epkmX_raw'].Draw()
hhs['epkmX_pass_c012'].Draw('same')
hhs['epkmX_pass_c012'].SetMinimum(0)
c3.SaveAs('hist_epkmX_cut_plots.png')


c3a=TCanvas('c3a','c3a',1000,1000)
c3a.Divide(2,2)
c3a.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].SetLineColor(kBlack)
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs['epkpkmXe_raw'].Draw()
hhs['epkpkmXe_pass_c123'].Draw('same')
c3a.cd(2)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].SetLineColor(kBlack)
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs['ekpkmX_raw'].Draw()
hhs['ekpkmX_pass_c013'].Draw('same')
hhs['ekpkmX_raw'].SetMinimum(0)
c3a.cd(3)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].SetLineColor(kBlack)
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs['epkpX_raw'].Draw()
hhs['epkpX_pass_c023'].Draw('same')
hhs['epkpX_pass_c023'].SetMinimum(0)
c3a.cd(4)
hhs['epkmX_raw'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].SetLineColor(kBlack)
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs['epkmX_raw'].Draw()
hhs['epkmX_pass_c012'].Draw('same')
hhs['epkmX_pass_c012'].SetMinimum(0)
c3a.SaveAs('hist_all_before_after.png')

sig=4
epkpkmxe_min = 0.08 - 0.0398*sig
epkpkmxe_max = 0.08 + 0.0398*sig

epkpX_min = 0.248 - 0.059*sig
epkpX_max = 0.248 + 0.059*sig

epkmX_min = 0.248 - 0.055*sig
epkmX_max = 0.248 + 0.055*sig

ekpkmX_min = 0.904 - 0.0719*sig
ekpkmX_max = 0.904 + 0.0719*sig






c3b=TCanvas('c3b','c3b',1000,1000)
c3b.Divide(2,2)
c3b.cd(1)
hhs['epkpkmXe_raw'].SetTitle('epK^{+}K^{-}X E Raw with cut lines;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].Draw() 
c3b.Update()
p1=c3b.GetPad(1)
ymax=gPad.GetUymax()
l1min = TLine(epkpkmxe_min,0.0,epkpkmxe_min,gPad.GetUymax())
l1max = TLine(epkpkmxe_max,0.0,epkpkmxe_max,gPad.GetUymax())
l1min.SetLineColor(kRed)
l1max.SetLineColor(kRed)
l1min.Draw('same')
l1max.Draw('same')

c3b.cd(2)
hhs['ekpkmX_raw'].SetTitle('eK^{+}K^{-}X Raw with cut lines;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].Draw()
c3b.Update()
p1=c3b.GetPad(1)
ymax=gPad.GetUymax()
l4min = TLine(ekpkmX_min,0.0,ekpkmX_min,ymax)
l4max = TLine(ekpkmX_max,0.0,ekpkmX_max,ymax)
l4min.SetLineColor(kRed)
l4max.SetLineColor(kRed)
l4min.Draw('same')
l4max.Draw('same')
c3b.cd(3)
hhs['epkpX_raw'].SetTitle('epK^{+}X Raw with cut lines;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].Draw()
c3b.Update()
p1=c3b.GetPad(1)
ymax=gPad.GetUymax()
l2min = TLine(epkpX_min,0.0,epkpX_min,ymax)
l2max = TLine(epkpX_max,0.0,epkpX_max,ymax)
l2min.SetLineColor(kRed)
l2max.SetLineColor(kRed)
l2min.Draw('same')
l2max.Draw('same')
c3b.cd(4)
hhs['epkmX_raw'].SetTitle('epK^{-}X Raw with cut lines;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].Draw()
c3b.Update()
p1=c3b.GetPad(1)
ymax=gPad.GetUymax()
l3min = TLine(epkmX_min,0.0,epkmX_min,ymax)
l3max = TLine(epkmX_max,0.0,epkmX_max,ymax)
l3min.SetLineColor(kRed)
l3max.SetLineColor(kRed)
l3min.Draw('same')
l3max.Draw('same')
c3b.SaveAs('hist_all_raw_with_cutlines.png')






c4 = TCanvas("c4","c4",900,900)
gStyle.SetOptStat(000000)
c4.cd(1)
hhs['vphimass_raw'].SetTitle("Invariant Mass of K^{+}K^{-}; mass (GeV); counts")
hhs['vphimass_raw'].SetLineColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetLineColor(kRed)
hhs['vphimass_raw'].Draw()
hhs['vphimass_pass_all_exclusivity'].Draw('same')
c4.SaveAs('hist_kpkm_cut_plots.png')

n_total_entries=hhs['epkmX_raw'].GetEntries()
print(n_total_entries)

hist_names = ('epkpkmXe','ekpkmX','epkpX','epkmX','vphimass')
excl_cut_names = ('missingE', 'epkpXmm2', 'ekpkmXmm2', 'epkmXmm2')
excl_cut_rates_all = []

for cn in excl_cut_names:
    #get blanked number of pass rates 
    n_pass = hhs['epkpX_'+cn].GetEntries()
    excl_cut_rates_all.append(n_pass)



d_excl_cut_rates_narrow={}
#get number in relevant region
for hn in hist_names:
    print(' looking at histogram %s ' % (hn))
    excl_cut_rates_narrow=[]
    for cn in excl_cut_names:
        n_pass_narrow=hhs[hn+'_'+cn].Integral(4,190,'')
        n_tot_narrow=hhs[hn+'_raw'].Integral(4,190,'')    
        print('for cut %s, %d pass out of %d ' % (cn,n_pass_narrow, n_tot_narrow))
        excl_cut_rates_narrow.append(n_pass_narrow)
    d_excl_cut_rates_narrow[hn]=excl_cut_rates_narrow



    


h_pass_rate = TH1F('h_pass_rate','h_pass_rate',len(excl_cut_names) ,0,len(excl_cut_names)-1)
bb=1
for p in excl_cut_rates_all:
    h_pass_rate.SetBinContent(bb,p)
    bb+=1

cp = TCanvas('cp','cp',900,900)
h_pass_rate.Draw()
cp.SaveAs('hist_pass_rates.png')

#plot rates for each histogrm in narrow region
for key,value in d_excl_cut_rates_narrow.items():

    h_pass_rate_narrow = TH1F('h_pass_rate_narrow_for_'+key,'h_pass_rate_narrow_for_'+key,len(excl_cut_names) ,0,len(excl_cut_names)-1)
    bb=1
    for p in value:
        h_pass_rate_narrow.SetBinContent(bb,p)
        bb+=1


    cpn = TCanvas('cpn'+key,'cpn'+key,900,900)
    h_pass_rate_narrow.Draw()
    cpn.SaveAs('hist_pass_rates_narrow_for'+key+'.png')


#ploet signal to background for each cut that is applied
phi_sb, phi_cut = array( 'd' ), array( 'd' ) 
cn_indx=0
for cn in excl_cut_names:
    nphi_signal = hhs['vphimass_'+cn].GetMaximum()#Integral(21,25,'')
    nphi_approx_bckgrnd = hhs['vphimass_'+cn].Integral(26,31,'')/6
    signal_to_bck = nphi_signal/nphi_approx_bckgrnd
    phi_sb.append(signal_to_bck)
    phi_cut.append(cn_indx)
    cn_indx+=1
    

g_phi_sb = TGraph(len(phi_cut),phi_cut,phi_sb)
g_phi_sb.SetTitle('S/B for K^{+}K^{-} per cut')
csb = TCanvas('csb','csb',900,900)
csb.cd(1)
g_phi_sb.SetMarkerColor( 4 )
g_phi_sb.SetMarkerStyle( 21 )
g_phi_sb.Draw('AP')
csb.SaveAs('vphimass_signal_to_bck.png')










