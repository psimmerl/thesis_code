from ROOT import TH1F, TFile, TCanvas
from ROOT import gStyle, gPad
import matplotlib.pyplot as plt
import sys


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
        fract = nev/htemp.GetEntries() * 100
        bbname = str(bb)+particle+"("+str(nev)+")"
        if( fract < 1 ): 
            bbname = ''
        pr_int_all_kaons_lab.append(bbname)
        n_pr_int_all_kaons.append(nev)

    ax1.pie(n_pr_int_all_kaons,labels=pr_int_all_kaons_lab, startangle = 140,autopct=autopct_generator(1))
    ax1.axis('equal')
    fig1.savefig('pie'+name)


fin2 = TFile(sys.argv[1],'open')

h_combsize = fin2.Get('combsize')

makePieChart(h_combsize,"multiplicity_of_combinations.png","comb size")




