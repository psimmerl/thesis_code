#!/bin/bash

#####################################################################
## Right now only set up for some weird blend of inb and inbNoutb? ##
#####################################################################

# #exclusive_phi/epkpkm_top
#X monitor-phi.py -i monitorphi*.root -o outname
#X mon-phi-selection.py -i monitorphi*.root -o outname
#X mon-phi-event.py -i monitorphi*.root -o outname
#X phimassV2.py -i inputfrommonitorAsy.py
#X mon-phimass.py -i inputfrommonitorAsy.py
#X makeFinalPlots.py -i input.root

# #monMC
# mon-data-mc.py
# mon-mc.py
# mon-sim-phi.py 

# #compare fields
#X mon-compare-fields.py -iinb monitorphiinbdata.root -ioutb outbending.root-o outname
#X mon-cf-phi-event.py -iinb monitorphiinbdata.root -ioutb outbending.root-o outname
#X mon-cf-phi-selection.py -iinb monitorphiinbdata.root -ioutb outbending.root-o outname 

# #phiasy
#X monitorAsy.py phimassresults*.root
#X mon-phi.py -i outputfrommonitorAsy.root -o outname
#X showBckAsy.py -o outputname
#X mon-phi-method2.py -i outputfrommonitorAsy.root -o outname

#Generate root files LATER already done!

#field type is not dynamically typed!
#check field type
#./run-groovy monitorPhi.groovy /volatile/clas12/kenjo/bclary/phi_skim/inb/*.hipo
#check field type
#./run-groovy monitorPhi.groovy /volatile/clas12/kenjo/bclary/phi_skim/outb/epKpKm_RGA_OUT.hipo

#
inb=$PWD/monitor_phi_skim8_005032_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
outb=$PWD/monitor_phi_epKpKm_RGA_OUT_pidtype1_cutvar1_rmvres12-outb_16core_incl.root

inbphi=$PWD/phi_mass_results_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
outbphi=$PWD/phi_mass_results_pidtype1_cutvar1_rmvres12-outb_16core_incl.root

#where did the v2 files coe from in mon-compare0fields
cd projects/compareFields
#python2.7 mon-compare-fields.py -iinb $inb -ioutb $outb -o mon-compare-fields
#python2.7 mon-cf-phi-event.py -iinb $inb -ioutb $outb -o mon-cf-phi-event
#python2.7 mon-cf-phi-selection.py -iinb $inb -ioutb $outb -o mon-cf-phi-selection

#Need to fix for inb and outb
cd ../exclusive_phi/epkpkm_top
#python2.7 makeFinalPlots.py $inb
# python2.7 makeFinalPlots.py $outb

#python2.7 monitor-phi.py $inb inbending
#python2.7 monitor-phi.py $outb outbending

#Note field setting is field_setting='inb' Need to change
python2.7 mon-phi-selection.py -i $inb -o mon-phi-selection-inb
#python2.7 monitor-phi-selection.py -i $outb -o monitor-phi-selection-outb

#Note field setting is field_setting='inbNoutb' Need to change
python2.7 mon-phi-event.py -i $inb -o mon-phi-event-inb
# python2.7 monitor-phi-event.py -i $outb -o monitor-phi-event-outb

#field_setting='inbNoutb' need to change?
#inbnoutb need to fix for two
#not sure where the bin_**.txt files came from
cd ../../phiasy
#monAsyinb=
#monAsyoutb=
monAsy=$PWD/mon_asy_pidtype4a_rmvres12_binningV11_cutvar0_-inbNoutb_v2.root
python2.7 monitorAsy.py $inbphi
#python2.7 monitorAsy.py $outbphi

#field inbNoutb
python2.7 mon-phi.py -i $monAsy -o mon-phi-inbNoutb
#inbNoutb
python2.7 showBckAsy.py -o showBckAsy-inbNoutb
#inbNoutb
python2.7 mon-phi-method2.py -i $monAsy -o mon-phi-method2-inbNoutb

cd ../exclusive_phi/epkpkm_top
#phiPeak_{}.png phiPeak_{}.pdf phiPeak_detailed_{}.png phiPeak_detailed_{}.pdf phiPeak_fom_{}.pdf
python2.7 phimassV2.py $monAsy inbNoutb

#
python2.7 mon-phimass.py -i $monAsy -o mon-phimass-inbNoutb


#monMC requires running simulation! Don't feel like setting up rn
