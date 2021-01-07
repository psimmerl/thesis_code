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
inb=monitor_phi_skim8_005032_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
outb=monitor_phi_epKpKm_RGA_OUT_pidtype1_cutvar1_rmvres12-outb_16core_incl.root

#phi_mass_results_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
#phi_mass_results_pidtype1_cutvar1_rmvres12-outb_16core_incl.root
#monAsyinb=
#monAsyoutb=
monAsy=mon_asy_pidtype4a_rmvres12_binningV11_cutvar0_-inbNoutb_v2.root


python2.7 projects/compareFields/mon-compare-fields.py -iinb $inb -ioutb $outb -o mon-compare-fields
python2.7 projects/compareFields/mon-cf-phi-event.py -iinb $inb -ioutb $outb -o mon-cf-phi-event
python2.7 projects/compareFields/mon-cf-phi-selection.py -iinb $inb -ioutb $outb -o mon-cf-phi-selection

#Need to fix for inb and outb
python2.7 projects/exclusive_phi/epkpkm_top/makeFinalPlots.py $inb
# python2.7 projects/exclusive_phi/epkpkm_top/makeFinalPlots.py $outb



python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi.py $inb inbending
python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi.py $inb outbending

#Note field setting is field_setting='inb' Need to change
python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi-selection.py -i $inb monitor-phi-selection-inb
#python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi-selection.py -i $inb monitor-phi-selection-outb

#Note field setting is field_setting='inbNoutb' Need to change
python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi-event.py -i $inb monitor-phi-event-inb
# python2.7 projects/exclusive_phi/epkpkm_top/monitor-phi-event.py -i $inb monitor-phi-event-outb

#field_setting='inbNoutb' need to change?
#inbnoutb need to fix for two
python2.7 projects/phiasy/monitorAsy.py phi_mass_results_pidtype1_cutvar1_rmvres12-inb_16core_incl.root
#python2.7 projects/phiasy/monitorAsy.py phi_mass_results_pidtype1_cutvar1_rmvres12-outb_16core_incl.root


#field inbNoutb
python2.7 projects/phiasy/mon-phi.py -i $monAsy -o mon-phi-inbNoutb
#inbNoutb
python2.7 projects/phiasy/showBckAsy.py -o showBckAsy-inbNoutb
#inbNoutb
python2.7 projects/phiasy/mon-phi-method2.py -i $monAsy -o mon-phi-method2-inbNoutb

#phiPeak_{}.png phiPeak_{}.pdf phiPeak_detailed_{}.png phiPeak_detailed_{}.pdf phiPeak_fom_{}.pdf
python2.7 projects/exclusive_phi/epkpkm_top/phimassV2.py $monAsy inbNoutb

#
python2.7 projects/exclusive_phi/epkpkm_top/mon-phimass.py -i $monAsy -o mon-phimass-inbNoutb


#monMC requires running simulation! Don't feel like setting up rn
