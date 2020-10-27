from ROOT import TH1F, TFile, TCanvas
from ROOT import gStyle, gPad
import matplotlib.pyplot as plt
import sys


fin = TFile(sys.argv[1],'open')


#select events with at least one proton + at least one kaonPlus + at least one kaonMinus          
h_multi_int = fin.Get('h_multiplicity_integrated')
h_multi_pr_int_allkaons= fin.Get('h_multiplicity_pr_integrated_allkaons')
h_multi_kp_int_allprkm = fin.Get('h_multiplicity_kp_integrated_allprkm')
h_multi_km_int_allprkp = fin.Get('h_multiplicity_km_integrated_allprkp')

#select events with only proton + at least one kaonPlus + at least one kaonMinus    
h_multi_singlepr_pr_int_allk = fin.Get('h_multiplicity_singlepr_pr_integrated_allkaons')
h_multi_singlepr_kp_int_allprkm = fin.Get('h_multiplicity_singlepr_kp_integrated_allprkm')
h_multi_singlepr_km_int_allprkp = fin.Get('h_multiplicity_singlepr_km_integrated_allprkp')

#select events with at least one proton + at least one kaonPlus    
h_mult_pr_topprkp = fin.Get('h_multiplicity_pr_integrated_topprkp')
h_mult_kp_topprkp = fin.Get('h_multiplicity_kp_integrated_topprkp')

#select events with at least one proton + at least one kaonMin
h_mult_pr_topprkm = fin.Get('h_multiplicity_pr_integrated_topprkp')
h_mult_km_topprkm = fin.Get('h_multiplicity_km_integrated_topprkm')

#select events with at least one kaonPlus + at least one kaonMinus       
h_mult_kp_topkpkm = fin.Get('h_multiplicity_kp_integrated_topkpkm')
h_mult_km_topkpkm = fin.Get('h_multiplicity_km_integrated_topkpkm')

gStyle.SetOptStat('e')

c1 = TCanvas("c1","c1",900,900) 
c1.cd(1)
gPad.SetLogy()
h_multi_pr_int_allkaons.SetTitle("Number of Protons for AL 1PR&1KP&1KM;Number of protons per event; Counts")
h_multi_pr_int_allkaons.Draw()
c1.SaveAs("number_of_protons_for_al_1pr1kp1km.png")

def autopct_generator(limit):
    def inner_autopct(pct):
        return ('%.2f' % pct) if pct > limit else ''
    return inner_autopct

def makePieChart(htemp , name, particle):
    pr_int_all_kaons_lab = []
    n_pr_int_all_kaons = []
    fig1, ax1 = plt.subplots()
    for bb in range(1,htemp.GetXaxis().GetNbins()):
        nev = htemp.GetBinContent(bb)
        fract = nev/htemp.GetEntries() *100
        bbname = str(bb)+particle
        if( fract < 1 ): 
            bbname = ''
        pr_int_all_kaons_lab.append(bbname)
        n_pr_int_all_kaons.append(nev)

    ax1.pie(n_pr_int_all_kaons,labels=pr_int_all_kaons_lab, startangle = 140,autopct=autopct_generator(1))
    ax1.axis('equal')
    fig1.savefig('pie'+name)

print("Making proton pie chart now")
makePieChart(h_multi_pr_int_allkaons,"number_of_protons_for_al_1pr1kp1km","Pr")

c2 = TCanvas("c2","c2",900,900)
c2.cd(1)
gPad.SetLogy()
h_multi_kp_int_allprkm.SetTitle("Number of K^{+} for AL 1PR&1KP&1KM;Number of K^{+} per event; Counts")
h_multi_kp_int_allprkm.Draw()
c2.SaveAs("number_of_kp_for_al_1pr1kp1km.png")
makePieChart(h_multi_kp_int_allprkm,"number_of_kp_for_al_1pr1kp1km","Kp")


c3 = TCanvas("c3","c3",900,900)
c3.cd(1)
gPad.SetLogy()
h_multi_km_int_allprkp.SetTitle("Number of K^{-} for AL 1PR&1KP&1KM;Number of K^{-} in event; Counts")
h_multi_km_int_allprkp.Draw()
c3.SaveAs("number_of_km_for_al_1pr1kp1km.png")
makePieChart(h_multi_km_int_allprkp,"number_of_km_for_al_1pr1kp1km","Km")

c4 = TCanvas("c4","c4",900,900)
c4.cd(1)
gPad.SetLogy()
h_multi_singlepr_kp_int_allprkm.SetTitle("Number of K^{+} for only 1PR& AL 1KP&1KM;Number of K^{+} per event;Counts")
h_multi_singlepr_kp_int_allprkm.Draw()
c4.SaveAs("number_of_kp_for_only1pr_al1kp1km.png")
makePieChart(h_multi_singlepr_kp_int_allprkm,"number_of_kp_for_only1pr_al1kp1km","Kp")


c5 = TCanvas("c5","c5",900,900)
c5.cd(1)
gPad.SetLogy()
h_multi_singlepr_km_int_allprkp.SetTitle("Number of K^{-} for only 1PR& AL 1KP&1KM;Number of K^{-} per event;Counts")
h_multi_singlepr_km_int_allprkp.Draw()
c5.SaveAs("number_of_km_for_only1pr_al1kp1km.png")
makePieChart(h_multi_singlepr_km_int_allprkp,"number_of_km_for_only1pr_al1kp1km","Km")

c6 = TCanvas("c6","c6",900,900)
c6.cd(1)
gPad.SetLogy()
h_mult_pr_topprkp.SetTitle("Number of Protons for AL 1PR&1KP;Number of Protons per event;Counts")
h_mult_pr_topprkp.Draw()
c6.SaveAs("number_of_pr_al1pr1kp.png")
makePieChart(h_mult_pr_topprkp,"number_of_pr_al1pr1kp","Pr")

c7 = TCanvas("c7","c7",900,900)
c7.cd(1)
gPad.SetLogy()
h_mult_kp_topprkp.SetTitle("Number of K^{+} for AL 1PR&1KP;Number of K^{+} per event;Counts")
h_mult_kp_topprkp.Draw()
c7.SaveAs("number_of_kp_al1pr1kp.png")
makePieChart(h_mult_kp_topprkp,"number_of_kp_al1pr1kp","Kp")

c8 = TCanvas("c8","c8",900,900)
c8.cd(1)
gPad.SetLogy()
h_mult_pr_topprkm.SetTitle("Number of Protons for AL 1PR&1KM;Number of Protons per event;Counts")
h_mult_pr_topprkm.Draw()
c8.SaveAs("number_of_pr_al1pr1km.png")
makePieChart(h_mult_pr_topprkm,"number_of_pr_al1pr1km","Pr")

c9 = TCanvas("c9","c9",900,900)
c9.cd(1)
gPad.SetLogy()
h_mult_km_topprkm.SetTitle("Number of K^{-} for AL 1PR&1KM;Number of K^{-} per event;Counts")
h_mult_km_topprkm.Draw()
c9.SaveAs("number_of_km_al1pr1km.png")
makePieChart(h_mult_km_topprkm,"number_of_km_al1pr1km","Km")

c10 = TCanvas("c10","c10",900,900)
c10.cd(1)
gPad.SetLogy()
h_mult_kp_topkpkm.SetTitle("Number of K^{+} for AL 1KP&1KM;Number of K^{+} per event; Counts")
h_mult_kp_topkpkm.Draw()
c10.SaveAs("number_of_kp_al1kp1km.png")
makePieChart(h_mult_kp_topkpkm,"number_of_kp_al1kp1km","Kp")

c11 = TCanvas("c11","c11",900,900)
c11.cd(1)
gPad.SetLogy()
h_mult_km_topkpkm.SetTitle("Number of K^{-} for AL 1KP&1KM;Number of K^{-} per event; Counts")
h_mult_km_topkpkm.Draw()
c11.SaveAs("number_of_km_al1kp1km.png")
makePieChart(h_mult_km_topkpkm,"number_of_km_al1kp1km","Km")











