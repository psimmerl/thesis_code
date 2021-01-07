#!/bin/bash

./run-groovy monitorPhi.groovy inb /volatile/clas12/kenjo/bclary/phi_skim/inb/*.hipo
./run-groovy monitorPhi.groovy outb /volatile/clas12/kenjo/bclary/phi_skim/outb/epKpKm_RGA_OUT.hipo
./run-groovy monitorPhi.groovy inbNoutb /volatile/clas12/kenjo/bclary/phi_skim/inb/*.hipo /volatile/clas12/kenjo/bclary/phi_skim/outb/epKpKm_RGA_OUT.hipo


#
inb=$PWD/monitor_phi_skim8_005032_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
outb=$PWD/monitor_phi_epKpKm_RGA_OUT_pidtype1_cutvar1_rmvres12-outb_16core_incl.root
inbNoutb=$PWD/monitor_phi_epKpKm_RGA_OUT_pidtype1_cutvar1_rmvres12-inbNoutb_16core_incl.root

inbphi=$PWD/phi_mass_results_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
outbphi=$PWD/phi_mass_results_pidtype1_cutvar1_rmvres12-outb_16core_incl.root
inbNoutbphi=$PWD/phi_mass_results_pidtype1_cutvar1_rmvres12-inbNoutb_16core_incl.root

#where did the v2 files coe from in mon-compare-fields
cd projects/compareFields
python2.7 mon-compare-fields.py -iinb $inb -ioutb $outb -o mon-compare-fields
python2.7 mon-cf-phi-event.py -iinb $inb -ioutb $outb -o mon-cf-phi-event
python2.7 mon-cf-phi-selection.py -iinb $inb -ioutb $outb -o mon-cf-phi-selection

cd ../exclusive_phi/epkpkm_top
python2.7 makeFinalPlots.py -i $inb -o pidtype1_inb_V2
python2.7 makeFinalPlots.py -i $outb -o pidtype1_outb_V2

python2.7 monitor-phi.py -i $inb -o inbending
python2.7 monitor-phi.py -i $outb -o outbending

python2.7 mon-phi-selection.py -i $inb -o mon-phi-selection-inb
python2.7 monitor-phi-selection.py -i $outb -o monitor-phi-selection-outb

python2.7 mon-phi-event.py -i $inb -o mon-phi-event-inb
python2.7 monitor-phi-event.py -i $outb -o monitor-phi-event-outb

#not sure where the bin_**.txt files came from
cd ../../phiasy
monAsyinb=$PWD/mon_asy_pidtype4a_rmvres12_binningV11_cutvar0_-inb_v2.root
monAsyoutb=$PWD/mon_asy_pidtype4a_rmvres12_binningV11_cutvar0_-outb_v2.root
monAsyinbNoutb=$PWD/mon_asy_pidtype4a_rmvres12_binningV11_cutvar0_-inbNoutb_v2.root
python2.7 monitorAsy.py -i $inbphi -o inb
python2.7 monitorAsy.py -i $outbphi -o outb
python2.7 monitorAsy.py -i $inbNoutbphi -o inbNoutb
##Supposed to run twice
python2.7 monitorAsy.py -i $inbphi -o inb
python2.7 monitorAsy.py -i $outbphi -o outb
python2.7 monitorAsy.py -i $inbNoutbphi -o inbNoutb

python2.7 mon-phi.py -i $monAsyinb -o mon-phi-inb
python2.7 mon-phi.py -i $monAsyoutb -o mon-phi-outb
python2.7 mon-phi.py -i $monAsyinbNoutb -o mon-phi-inbNoutb

python2.7 showBckAsy.py -o showBckAsy-inb
python2.7 showBckAsy.py -o showBckAsy-outb
python2.7 showBckAsy.py -o showBckAsy-inbNoutb

python2.7 mon-phi-method2.py -i $monAsyinb -o mon-phi-method2-inb
python2.7 mon-phi-method2.py -i $monAsyoutb -o mon-phi-method2-outb
python2.7 mon-phi-method2.py -i $monAsyinbNoutb -o mon-phi-method2-inbNoutb

cd ../exclusive_phi/epkpkm_top
#phiPeak_{}.png phiPeak_{}.pdf phiPeak_detailed_{}.png phiPeak_detailed_{}.pdf phiPeak_fom_{}.pdf
python2.7 phimassV2.py -i $monAsyinb -o inb
python2.7 phimassV2.py -i $monAsyoutb -o outb
python2.7 phimassV2.py -i $monAsyinbNoutb -o inbNoutb

#
python2.7 mon-phimass.py -i $monAsyinb -o mon-phimass-inb
python2.7 mon-phimass.py -i $monAsyoutb -o mon-phimass-outb
python2.7 mon-phimass.py -i $monAsyinbNoutb -o mon-phimass-inbNoutb


#monMC requires running simulation! Don't feel like setting up rn
