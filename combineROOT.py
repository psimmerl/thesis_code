import argparse
from ROOT import TFile, TH1F, TH2F, TF1, TLine, TLatex, TLegend, TBox
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue, kMagenta, kGreen
from ROOT import TGraph
import sys,os, math
from array import array

ap = argparse.ArgumentParser()#Paul addition
ap.add_argument(
    '-iinb',
    '--input_inb_file',
    required=True
)
ap.add_argument(
    '-ioutb',
    '--input_outb_file',
    required=True
)
ap.add_argument(
    '-o',
    '--inbNoutb_file',
    required=True
)
args = ap.parse_args()

finb = TFile(args.input_inb_file)
foutb = TFile(args.input_outb_file)

hhs = {}
for kk in finb.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj

for kk in foutb.GetListOfKeys():
    obj = kk.ReadObj()
    if obj.GetName() in hhs.keys():
        hhs[obj.GetName()].Add(obj)
    else:
        hhs[obj.GetName()] = obj


finbNoutb = TFile(args.inbNoutb_file, 'recreate')
for _, hist in hhs.items():
    hist.Write()

finb.Close()
foutb.Close()
finbNoutb.Close()
