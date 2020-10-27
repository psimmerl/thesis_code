package pid.electron

import org.jlab.detector.base.DetectorType
import event.Event

class ElectronFromEvent {

    def ebeam = 10.6
    def ebPID = 11
    def min_nphe = 2
    def nphe_loose = -1
    def nphe_med = 0
    def nphe_tight = 1
    
    //vertex wide enough for all sectors
    /// def min_vz = -12
    /// def max_vz = 9

    // technically the pcal energy
    def min_ecal_inner_dep = 0.07
    def ecal_inner_dep_loose = -0.01
    def ecal_inner_dep_med = 0.0
    def ecal_inner_dep_tight = 0.02

    // field setting to toggle kind of cut for electrons
    def mag_field = null

    //calorimeter fiducial cut
    //def min_u = 0
    //def max_u = 420
    //def min_v = 0
    //def max_v = 420
    //def min_w = 0
    //def max_w = 420

    
    //def s_min_u_outb = [33, 26, 34, 22, 27, 25]
    //def s_max_u_outb = [408, 408, 408, 408, 408, 408]
    
    //def s_min_v_outb = [16.0, 10.5, 17.0, 14.25, 18.0, 11.0]
    //def s_max_v_outb = [400, 400, 400, 400, 400, 400]
    
    //def s_min_w_outb = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    //def s_max_w_outb = [400, 400, 400, 400, 400, 400]

    def min_u_tight_inb = [4.5,  4.5,  4.5,  4.5,  4.5,  4.5 ]
    def min_u_med_inb   = [9.0,  9.0,  9.0,  9.0,  9.0,  9.0 ]
    def min_u_loose_inb = [13.5, 13.5, 13.5, 13.5, 13.5, 13.5]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_inb = [398, 398, 398, 398, 398, 398]
    def max_u_med_inb   = [408, 408, 408, 408, 408, 408]
    def max_u_loose_inb = [420, 420, 420, 420, 420, 420] 
    
    def min_u_cuts_inb = [min_u_loose_inb, min_u_med_inb, min_u_tight_inb]
    def max_u_cuts_inb = [max_u_loose_inb, max_u_med_inb, max_u_tight_inb]

    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_inb = [19,19,19,19,19,19] //[4.5,  4.5,  4.5,  4.5,  4.5,  4.5]
    def min_v_med_inb   = [14,14,14,14,14,14] //[9.0,  9.0,  9.0,  9.0,  9.0,  9.0]
    def min_v_loose_inb = [9,9,9,9,9,9]       //13.5, 13.5, 13.5, 13.5, 13.5, 13.5]
    
    // the maximum of v is never reached
    def max_v_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_v_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_inb = [min_v_loose_inb, min_v_med_inb, min_v_tight_inb]
    def max_v_cuts_inb = [max_v_loose_inb, max_v_med_inb, max_v_tight_inb]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    // flip the loose and medium because of weird combos in rga note
    def min_w_tight_inb = [19,19,19,19,19,19]                 // [4.5,  4.5,  4.5,  4.5,  4.5,  4.5]
    def min_w_med_inb   = [14,14,14,14,14,14]                 //[9.0,  9.0,  9.0,  9.0,  9.0,  9.0]
    def min_w_loose_inb = [9,9,9,9,9,9]                       // [13.5, 13.5, 13.5, 13.5, 13.5, 13.5]
    
    // the maximum of w is never reached
    def max_w_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_w_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_w_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_w_cuts_inb = [min_w_loose_inb, min_w_med_inb, min_w_tight_inb]
    def max_w_cuts_inb = [max_w_loose_inb, max_w_med_inb, max_w_tight_inb]
   
    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_u_tight_out = [42, 32, 38, 27.5, 32, 29]
    def min_u_med_out   = [33, 26, 34, 22,   27, 25]
    def min_u_loose_out = [28, 22, 30, 18,   22, 22]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_out = [398, 398, 398, 398, 398, 398]
    def max_u_med_out   = [408, 408, 408, 408, 408, 408] 
    def max_u_loose_out = [420, 420, 420, 420, 420, 420] 

    def min_u_cuts_outb = [min_u_loose_out, min_u_med_out, min_u_tight_out]
    def max_u_cuts_outb = [max_u_loose_out, max_u_med_out, max_u_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_out = [18.0, 12.0, 19.5,  15.5,  20.0, 13.0]
    def min_v_med_out   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
    def min_v_loose_out = [10.25, 8.0, 12.75, 12.5,  13.25, 9.0]
    
    // the maximum of v is never reached
    def max_v_tight_out = [400, 400, 400, 400, 400, 400]
    def max_v_med_out   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_out = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_outb = [min_v_loose_out, min_v_med_out, min_v_tight_out]
    def max_v_cuts_outb = [max_v_loose_out, max_v_med_out, max_v_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_w_tight_out = [14.0, 18.7, 18.7,  12.0, 16.0, 13.0]
    def min_w_med_out   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    def min_w_loose_out = [7.25, 11.0, 13.0,  5.5,  10.0, 6.0]
    
    // the maximum of w is never reached
    def max_w_tight_out = [400, 400, 400, 400, 400, 400]
    def max_w_med_out   = [400, 400, 400, 400, 400, 400]
    def max_w_loose_out = [400, 400, 400, 400, 400, 400]

    def min_w_cuts_outb = [min_w_loose_out, min_w_med_out, min_w_tight_out]
    def max_w_cuts_outb = [max_w_loose_out, max_w_med_out, max_w_tight_out]

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    def sigma_range = 3.5// to match rga note 08/2020 before it was 4.0
    def p0mean_inb = [0.116405, 0.118783, 0.115558, 0.112439, 0.120496, 0.117107]
    def p1mean_inb = [-0.0155433, 0.109417, -0.0088641, -0.0704795, 0.100674, -0.0164561]
    def p2mean_inb = [0.00722919, 0.0064115, 0.00773618, 0.00769267, 0.00393716, 0.00702487]
    def p3mean_inb = [-0.000859442, -0.000704544, -0.000853684, -0.00073044, -0.000403603, -0.000825457]
    def p0sigma_inb = [0.0164086, 0.00448347, 0.00770184, 0.00412631, 0.00687372, 0.00747944]
    def p1sigma_inb = [0.00771924, 0.020591, 0.0164028, 0.0217636, 0.0183004, 0.0184456]
    def p2sigma_inb = [-0.00125094, 0.000362296, -0.000381947, 0.000232047, -0.000439471, -0.000560982]
    def p3sigma_inb = [8.43022e-05, -2.94663e-05, 2.47408e-05, -5.70412e-06, 5.14683e-05, 5.80948e-05]

    /*
    def p0mean_outb = [0.105467, 0.115261, 0.127793, 0.113359, 0.112263, 0.113507]
    def p1mean_outb = [-0.135178, 0.135808, 0.903412, 0.598274, -0.0466815, 0.0550123]
    def p2mean_outb = [0.00996842, 0.00672508, -0.00035721, 0.00470925, 0.00588451, 0.00923385]
    def p3mean_outb = [-0.000754536, -0.000515365, -0.000108273, -0.000447278, -0.000358148, -0.00074643]
    def p0sigma_outb = [0.00683747, 0.0065199, 0.00297734, 0.00759701, 0.0093309, 0.00591988]
    def p1sigma_outb = [0.0180228, 0.0183979, 0.0250332, 0.0155001, 0.0137594, 0.0215643]
     */
    def p0mean_outb = [0.11946, 0.120319, 0.121442, 0.117984, 0.126656, 0.120343]
    def p1mean_outb = [-0.00388941, 0.566384, 0.346713, -0.100043, 0.766969, 0.394986]
    def p2mean_outb = [0.00693105, 0.00702178, 0.0040979, 0.00381667, 0.000243406, 0.00563105]
    def p3mean_outb = [-0.000721964, -0.000688366, -0.000459883, -0.000359391, -0.000113102, -0.000620295]
    def p0sigma_outb = [0.00831965, 0.0151432, 0.0122673, 0.0114592, 0.014943, 0.00711617]
    def p1sigma_outb = [0.0156242, 0.00691968, 0.0106717, 0.0116251, 0.00837599, 0.0158003]
    def p2sigma_outb = [-0.000538767, -0.000810607, -0.00120787, -0.000734793, -0.00160916, -2.0166e-05]
    def p3sigma_outb = [3.91481e-05, 2.00367e-05, 7.92897e-05, 2.99526e-05, 0.00011959, -2.08114e-05]

    def p0mean = null
    def p1mean = null
    def p2mean = null
    def p3mean = null
    def p0sigma = null
    def p1sigma = null
    def p2sigma = null
    def p3sigma = null


    //dcr1,2,3 fiducial cut
    //require functional form
    def sect_angle_coverage = 60
    
    

    //////////////////////////////////////////////
    // use the tracking bank to get the sector for tracks

    def dcR1_loose = -1
    def dcR1_med = 0
    def dcR1_tight = 1

    def dcR2_loose = -2
    def dcR2_med = 0
    def dcR2_tight = 2

    def dcR3_loose = -3
    def dcR3_med = 0
    def dcR3_tight = 3

    def radiusR1_inb = 29
    def radiusR1_outb = 29
       
    def heightR1_inb = [27, 22, 22, 27, 22, 22]
    def heightR1_outb = [15, 15, 15, 15, 15, 15]
    
    def radiusR2_inb = 38
    def radiusR2_outb = 39

    def heightR2_inb  = [40, 34, 34, 40, 34, 34]
    def heightR2_outb = [25, 25, 25, 25, 25, 25]

    def radiusR3_inb = 50
    def radiusR3_outb = 65
    
    def heightR3_inb = [47, 39, 39, 47, 39, 39]
    def heightR3_outb = [48, 48, 48, 48, 48, 48]

    // vertex position cuts
    def vz_tight = -1 
    def vz_med  = 0
    def vz_loose = 1

    def vz_min_sect_inb = -13
    def vz_max_sect_inb = 12 
    // old and dated below
    //def vz_min_sect_inb = [-9.0, -9.0, -9.2, -9.0, -8.8, -9.4]
    //def vz_max_sect_inb = [5.5,  6.5,  5.0,  6.4,  6.0,  6.0]
    
    // not sector dependent as of 08/2020
    def vz_min_sect_outb = -18 //[-14, -14, -14, -14, -14, -14] 
    def vz_max_sect_outb = 10 //[5, 5, 5, 5, 5, 5]

    def dc_fiducial_cutXY_min = [[[7.62814, -0.56319],[18.2833, -0.587275],[20.2027, -0.54605]],
				 [[9.20907, -0.586977],[10.493, -0.544243],[23.0759, -0.581959]],
				 [[12.5459, -0.631322],[20.5635, -0.618555],[26.3621, -0.576806]],
				 [[8.36343, -0.552394],[14.7596, -0.554798],[29.5554, -0.60545]],
				 [[16.3732, -0.663303],[10.0255, -0.533019],[31.6086, -0.617053]],
				 [[8.20222, -0.567211],[20.0181, -0.605458],[22.2098, -0.567599]]]

    def dc_fiducial_cutXY_max = [[[-7.49907,0.583375],[-18.8174,0.599219],[-23.9353,0.574699]],
				 [[-14.0547,0.631533],[-14.4223,0.597079],[-14.838,0.547436]],
				 [[-7.72508,0.578501],[-18.7928,0.56725],[-29.9003,0.612354]],
				 [[-6.12844,0.566777],[-13.6772,0.573262],[-26.1895,0.591816]],
				 [[-20.0718,0.670941],[-9.4775,0.511748],[-28.0869,0.590488]],
				 [[-9.52924,0.591687],[-17.8564,0.596417],[-23.5661,0.576317]]]

    def dc_fiducial_cutThetaPhi_min = [
	[[[37.289, -27.5201,1.12866, -0.00526111],[45.3103, -33.5226,1.72923, -0.0114495],[61.5709, -47.6158,3.4295, -0.0316429]],
	 [[36.6259, -27.4064,1.16617, -0.00604629],[50.3751, -37.5848,2.19621, -0.0169241],[35.1563, -26.514,1.09795, -0.00545864]],
	 [[27.2367, -20.3068,0.517752, -0.000335432],[39.0489, -28.6903,1.24306, -0.0065226],[41.0208, -30.0339,1.30776, -0.00626721]],
	 [[29.261, -21.7041,0.613556, -0.000774652],[39.5304, -29.1388,1.34116, -0.00823818],[44.5313, -33.4056,1.77581, -0.0123965]],
	 [[36.5659, -25.119,0.714074, -2.65397e-11],[31.6524, -22.6934,0.613977, -5.46634e-10],[34.7312, -24.9901,0.749061, -1.22922e-09]],
	 [[33.154, -23.8803,0.685794, -1.13236e-10],[42.6731, -31.0799,1.40425, -0.00730816],[46.4732, -35.6988,2.10144, -0.0164771]]]
    ]

    def dc_fiducial_cutThetaPhi_max = [
	[[[-35.1716, 25.102, -0.750281, 5.34679e-05],[-39.1633, 28.5551, -1.13429, 0.00419047],[-33.7705, 24.8068, -0.811239, 0.00138345]],
	 [[-36.2389, 26.7979, -1.08147, 0.0050898],[-43.643, 31.6783, -1.49203, 0.00872922],[-54.4042, 40.6516, -2.52393, 0.0205649]],
	 [[-38.3238, 26.1667, -0.777077, 0.000264835],[-34.2011, 24.2843, -0.696392, 3.75866e-12],[-36.4636, 25.8712, -0.786592, 2.24421e-10]],
	 [[-31.8019, 23.154, -0.653992, 2.69968e-05],[-34.6637, 24.6043, -0.714901, 2.02675e-10],[-36.7209, 26.2469, -0.828638, 0.000340435]],
	 [[-33.4016, 24.6901, -0.779889, 0.000430557],[-35.4583, 24.7491, -0.707953, 2.18559e-10],[-37.7335, 28.1547, -1.1986, 0.00582395]],
	 [[-34.7808, 24.6988, -0.719936, 5.73299e-10],[-54.5797, 40.9138, -2.57493, 0.0213354],[-38.4972, 28.3142, -1.21741, 0.00640373]]]
    ]


    // updated sampluing fraction cuts
    def ecal_e_sampl_mu = [
        [  0.2531 ,  0.2550 ,  0.2514 ,  0.2494 ,  0.2528 ,  0.2521 ] , 
        [ -0.6502 , -0.7472 , -0.7674 , -0.4913 , -0.3988 , -0.703  ] , 
        [  4.939  ,  5.350  ,  5.102  ,  6.440  ,  6.149  ,  4.957  ]]
    
    def ecal_e_sampl_sigm = [
        [  2.726e-3 ,  4.157e-3 ,  5.222e-3 ,  5.398e-3 ,  8.453e-3 ,  6.533e-3 ],
        [  1.062    ,  0.859    ,  0.5564   ,  0.6576   ,  0.3242   ,  0.4423   ], 
        [ -4.089    , -3.318    , -2.078    , -2.565    , -0.8223   , -1.274    ]]
    

    def s_min_u = null
    def s_max_u = null
    def s_min_v = null
    def s_max_v = null
    def s_min_w = null
    def s_max_w = null
      
    def heightR1=null
    def heightR2=null
    def heightR3=null

    def radiusR1=null
    def radiusR2=null
    def radiusR3=null

    def min_vz=null
    def max_vz=null
    def p_min=null

    def el_cut_strictness_lvl=null
    def ecal_cut_lvl = null
    def nphe_cut_lvl = null
    def vz_cut_lvl = null
    def min_u_cut_lvl = null
    def min_v_cut_lvl = null
    def min_w_cut_lvl = null
    def max_u_cut_lvl = null
    def max_v_cut_lvl = null
    def max_w_cut_lvl = null
    def dcr1_cut_lvl = null
    def dcr2_cut_lvl = null
    def dcr3_cut_lvl = null

    def anti_pion_threshold=null
            
    void setElectronCutStrictness(el_cut_strictness){
	el_cut_strictness_lvl=el_cut_strictness
	ecal_cut_lvl=el_cut_strictness["ecal_cut_lvl"]
	nphe_cut_lvl=el_cut_strictness["nphe_cut_lvl"]
	vz_cut_lvl=el_cut_strictness["vz_cut_lvl"]
	min_u_cut_lvl=el_cut_strictness["min_u_cut_lvl"]
	min_v_cut_lvl=el_cut_strictness["min_v_cut_lvl"]
	min_w_cut_lvl=el_cut_strictness["min_w_cut_lvl"]
	max_u_cut_lvl=el_cut_strictness["max_u_cut_lvl"]
	max_v_cut_lvl=el_cut_strictness["max_v_cut_lvl"]
	max_w_cut_lvl=el_cut_strictness["max_w_cut_lvl"]
	dcr1_cut_lvl=el_cut_strictness["dcr1_cut_lvl"]
	dcr2_cut_lvl=el_cut_strictness["dcr2_cut_lvl"]
	dcr3_cut_lvl=el_cut_strictness["dcr3_cut_lvl"]
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> el_cut_strictness " + el_cut_strictness)
	println("[ElectronFromEvent::setElectronCutStrictness] -> ecal_cut_lvl " + ecal_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> nphe_cut_lvl " + nphe_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> vz_cut_lvl " + vz_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_u_cut_lvl " + min_u_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_v_cut_lvl " + min_v_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_w_cut_lvl " + min_w_cut_lvl)
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_u_cut_lvl " + max_u_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_v_cut_lvl " + max_v_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_w_cut_lvl " + max_w_cut_lvl)
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr1_cut_lvl " + dcr1_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr2_cut_lvl " + dcr2_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr3_cut_lvl " + dcr3_cut_lvl)
	
	if( ecal_cut_lvl == 0 ) { min_ecal_inner_dep += ecal_inner_dep_loose }
	else if( ecal_cut_lvl == 1 ) { min_ecal_inner_dep += ecal_inner_dep_med }
	else if( ecal_cut_lvl == 2 ) { min_ecal_inner_dep += ecal_inner_dep_tight }

	if( nphe_cut_lvl == 0 ){ min_nphe += nphe_loose }
	else if( nphe_cut_lvl == 1 ){ min_nphe += nphe_med }
	else if( nphe_cut_lvl == 2 ){ min_nphe += nphe_tight }

	if( vz_cut_lvl == 0 ) { 
	    vz_min_sect_inb = vz_min_sect_inb-1.6// vz_min_sect_inb.collect{ it - 1.6 }
	    vz_max_sect_inb = vz_max_sect_inb+2.5//.collect{ it + 2.5 } 

	    vz_min_sect_outb = vz_min_sect_outb-vz_loose//.collect{ it - vz_loose }
	    vz_max_sect_outb = vz_max_sect_outb+vz_loose//.collect{ it + vz_loose } 
	}
	else if( vz_cut_lvl == 1 ) { 
	    vz_min_sect_inb = vz_min_sect_inb-0//.collect{ it - 0 }
	    vz_max_sect_inb = vz_max_sect_inb+0//.collect{ it + 0 } 

	    vz_min_sect_outb = vz_min_sect_outb-vz_med//.collect{ it - vz_med }
	    vz_max_sect_outb = vz_max_sect_outb+vz_med//.collect{ it + vz_med }  
	}
	else if( vz_cut_lvl == 2 ) { 
	    vz_min_sect_inb = vz_min_sect_inb+1//.collect{ it + 1 }
	    vz_max_sect_inb = vz_max_sect_inb-1.9//.collect{ it - 1.9 } 

	    vz_min_sect_outb = vz_min_sect_outb-vz_tight//.collect{ it - vz_tight }
	    vz_max_sect_outb = vz_max_sect_outb+vz_tight//.collect{ it + vz_tight } 
	}

	if( dcr1_cut_lvl == 0 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_loose }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_loose }
	    radiusR1_inb+=dcR1_loose
	    radiusR1_outb+=dcR1_loose
	}
	else if( dcr1_cut_lvl == 1 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_med }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_med }
	    radiusR1_inb+=dcR1_med
	    radiusR1_outb+=dcR1_med
	}
	else if( dcr1_cut_lvl == 2 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_tight }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_tight }
	    radiusR1_inb+=dcR1_tight
	    radiusR1_outb+=dcR1_tight
	}
	if( dcr2_cut_lvl == 0 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_loose }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_loose }
	    radiusR2_inb+=dcR2_loose
	    radiusR2_outb+=dcR2_loose
	}
	else if( dcr2_cut_lvl == 1 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_med }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_med }
	    radiusR2_inb+=dcR2_med
	    radiusR2_outb+=dcR2_med
	}
	else if( dcr2_cut_lvl == 2 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_tight }
	     heightR2_outb = heightR2_outb.collect{it + dcR2_tight }
	    radiusR2_inb+=dcR2_tight
	    radiusR2_outb+=dcR2_tight
	}

	if( dcr3_cut_lvl == 0 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_loose }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_loose }
	    radiusR3_inb+=dcR3_loose
	    radiusR3_outb+=dcR3_loose
	}
	else if( dcr3_cut_lvl == 1 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_med }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_med }
	    radiusR3_inb+=dcR3_med
	    radiusR3_outb+=dcR3_med
	}
	else if( dcr3_cut_lvl == 2 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_tight }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_tight }
	    radiusR3_inb+=dcR3_tight
	    radiusR3_outb+=dcR3_tight
	}	    	
	
    }

    void setElectronCutParameters(magnetic_field_config){
	println(' [ElectronFromEvent::setElectronCutParameters] -> setting electron cut parameters for field ' + magnetic_field_config)
	if( magnetic_field_config == "outbending" ){
	    mag_field = 1
	    println(' [ElectronFromEvent::setElectronCutParameters] -> setting parameters for outbending')
	    s_min_u=min_u_cuts_outb[min_u_cut_lvl] //min_u_med_out
	    s_min_v=min_v_cuts_outb[min_v_cut_lvl]
	    s_min_w=min_w_cuts_outb[min_w_cut_lvl]
	    
	    s_max_u=max_u_cuts_outb[max_u_cut_lvl]//max_u_med_out
	    s_max_v=max_v_cuts_outb[max_v_cut_lvl]//max_v_med_out
	    s_max_w=max_w_cuts_outb[max_w_cut_lvl]//max_w_med_out
	    
	    heightR1=heightR1_outb
	    heightR2=heightR2_outb
	    heightR3=heightR3_outb
	    
	    radiusR1=radiusR1_outb
	    radiusR2=radiusR2_outb
	    radiusR3=radiusR3_outb

	    min_vz=vz_min_sect_outb
	    max_vz=vz_max_sect_outb

	    p0mean  = p0mean_outb
	    p1mean  = p1mean_outb
	    p2mean  = p2mean_outb
	    p3mean  = p3mean_outb
	    p0sigma = p0sigma_outb
	    p1sigma = p1sigma_outb
	    p2sigma = p2sigma_outb
	    p3sigma = p3sigma_outb
	    p_min=0.2381 + 0.11905*ebeam
	    anti_pion_threshold = 0.2
	}
	else if( magnetic_field_config == "inbending" ){
	    mag_field = 0
	    println(' [ElectronFromEvent::setElectronCutParameters] -> setting parameters for inbending')
	    s_min_u=min_u_cuts_inb[min_u_cut_lvl] //min_u_med_out
	    s_min_v=min_v_cuts_inb[min_v_cut_lvl]
	    s_min_w=min_w_cuts_inb[min_w_cut_lvl]
	    
	    s_max_u=max_u_cuts_inb[max_u_cut_lvl]//max_u_med_out
	    s_max_v=max_v_cuts_inb[max_v_cut_lvl]//max_v_med_out
	    s_max_w=max_w_cuts_inb[max_w_cut_lvl]//max_w_med_out

	    heightR1=heightR1_inb
	    heightR2=heightR2_inb
	    heightR3=heightR3_inb

	    radiusR1=radiusR1_inb
	    radiusR2=radiusR2_inb
	    radiusR3=radiusR3_inb

	    min_vz=vz_min_sect_inb
	    max_vz=vz_max_sect_inb

	    p0mean  = p0mean_inb
	    p1mean  = p1mean_inb
	    p2mean  = p2mean_inb
	    p3mean  = p3mean_inb
	    p0sigma = p0sigma_inb
	    p1sigma = p1sigma_inb	
	    p2sigma = p2sigma_inb
	    p3sigma = p3sigma_inb
	    
	    p_min=0.2381 + 0.11905*ebeam
	    anti_pion_threshold=0.2
	}


	println('[ElectronFromEvent::setElectronCutParameters] -> electron pdg pid ' + ebPID )
	println('[ElectronFromEvent::setElectronCutParameters] -> min pcal edep ' + min_ecal_inner_dep )
	println('[ElectronFromEvent::setElectronCutParameters] -> min nphe htcc ' + min_nphe )

	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL U limits per sector : ' + s_min_u )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL U limits per sector : ' + s_max_u )
	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL V limits per sector : ' + s_min_v )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL V limits per sector : ' + s_max_v )
	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL W limits per sector : ' + s_min_w )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL W limits per sector : ' + s_max_w )

	println('[ElectronFromEvent::setElectronCutParameters] -> DC R1 height ' + heightR1 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R2 height ' + heightR2 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R3 height ' + heightR3 )

	println('[ElectronFromEvent::setElectronCutParameters] -> DC R1 radius ' + radiusR1 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R2 radius ' + radiusR2 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R3 radius ' + radiusR3 )

	println('[ElectronFromEvent::setElectronCutParameters] -> min vertexZ limits per sector ' + min_vz )
	println('[ElectronFromEvent::setElectronCutParameters] -> max vertexZ limits per sector ' + max_vz )

	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p0 mean ' + p0mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p1 mean ' + p1mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p2 mean ' + p2mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p3 mean ' + p3mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p1 sigma ' + p0sigma )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p2 sigma ' + p1sigma )

    }

    //////////////////////////////////////////////

    //def heightR2 = 35
    //def radiusR2 = 40

    //def heightR3 = 48
    //def radiusR3 = 49

    ///////////////////////////////////////////////////////////////////////////////////////
        
    def passElectronStatus = { event, index ->
        return (event.status[index] < 0)
    }

    def passElectronEBPIDCut = { event, index ->
	return (event.pid[index] == ebPID)
    }

    def passElectronChargeCut = { event, index ->
	return (event.charge[index] < 0)
    }

    def passElectronNpheCut = { event, index ->
	return (event.cherenkov_status.contains(index)) ? event.nphe[index] > min_nphe : false
    }

    // old cuts
    //def passElectronVertexCut = { event, index ->
    //if (event.pcal_status.contains(index)){	  	
    //def sec = event.pcal_sector[index]-1 // ? do we want the pcal sector or a drift chamber sector? - too bad using pcal sector
    //      return (event.vz[index].with{ it < max_vz[sec] && it > min_vz[sec] })
    //	}
    //	return false
    //}

    def passElectronVertexCut = { event, index ->
        return (event.vz[index].with{ it < max_vz && it > min_vz })		
    }


    def passElectronTrackQualityCut = { event, index ->
	if(event.cherenkov_status.contains(index) && 
	   event.pcal_status.contains(index) && 
	   event.dc1_status.contains(index) ){
	    if( event.pcal_sector[index]-1 > 0 && event.dc_sector[index]-1 > 0 && event.cherenkov_sector[index]-1 > 0){ //forces tracks to be in Forward detector
		return true	    
	    }
	}
	return false
    }

    def passElectronPCALFiducialCut = { event, index ->
	if (event.pcal_status.contains(index)){	  
	    def pcal_sector = event.pcal_sector[index]-1
	    //return (event.pcal_u[index].with{it < s_max_u[pcal_sector] && it > s_min_u[pcal_sector] } &&
	    //event.pcal_v[index].with{it < s_max_v[pcal_sector] && it > s_min_v[pcal_sector] } &&
	    //event.pcal_w[index].with{it < s_max_w[pcal_sector] && it > s_min_w[pcal_sector] })
	    return (event.pcal_v[index].with{it > s_min_v[pcal_sector] } &&
		    event.pcal_w[index].with{it > s_min_w[pcal_sector] })

	}
	return false
    }

    def passElectronSamplingFractionCut = { event, index ->
	if( event.ecal_inner_status.contains(index) ||  event.ecal_outer_status.contains(index) || event.pcal_status.contains(index) ){
	    def eidep=0
	    def eodep = 0
	    def pcaldep = 0
	    def sector = -1
	    if( event.ecal_inner_status.contains(index) ){
		eidep = event.ecal_inner_energy[index]
		//sector = event.ecal_inner_sector[index] -1
	    }
	    if( event.ecal_outer_status.contains(index) ){
		eodep = event.ecal_outer_energy[index]
		//sector = event.ecal_outer_sector[index]-1
	    }
	    if( event.pcal_status.contains(index) ){
		pcaldep = event.pcal_energy[index]
		sector = event.pcal_sector[index]-1
	    }
	    
	    def edep = eidep + eodep + pcaldep 
	    if( sector >= 0 ){
		
		def p = event.p[index]
		//def mean = p0mean[sector] * (1 + p/Math.sqrt(p*p + p1mean[sector])) + p2mean[sector]*p + p3mean[sector]*p*p
		//def sigma = p0sigma[sector] + p1sigma[sector]/Math.sqrt(p) + p2sigma[sector]*p + p3sigma[sector]*p*p

		def mean = ecal_e_sampl_mu[0][sector] + ecal_e_sampl_mu[1][sector]/1000*Math.pow(p - ecal_e_sampl_mu[2][sector],2)
		def sigma = ecal_e_sampl_sigm[0][sector] + ecal_e_sampl_sigm[1][sector]/(10*(p - ecal_e_sampl_sigm[2][sector]))

		def upper_cut = mean + sigma_range * sigma
		def lower_cut = mean - sigma_range * sigma
		
		if( edep/p <= upper_cut && edep/p >= lower_cut ) return true	    
	    }
	}
	return false
    }

    def passElectronAntiPionCut = { event, index ->
	if(  event.p[index] > 4.5 ){
	    if( event.ecal_inner_status.contains(index) && event.pcal_status.contains(index) ){	  
		return  -event.pcal_energy[index]/event.p[index] + anti_pion_threshold < event.ecal_inner_energy[index]/event.p[index]
	    }
	    else{
		return false
	    }
	}
	return true
    }


    def passElectronMinMomentum = { event, index ->
	return (event.p[index] > p_min )
    }

    def passElectronEIEOCut = { event, index ->
	return (event.pcal_status.contains(index)) ? event.pcal_energy[index] > min_ecal_inner_dep : false
    }

    //detector layer r1-6, r2-18, r3-36
    //rotate hit position based on sector
    def rotateDCHitPosition(hit, sec) {
        def ang = Math.toRadians(sec * sect_angle_coverage)
        def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
        def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
        return [x1_rot, y1_rot]
    }

    //define left right 
    def borderDCHitPosition(y_rot, height) {
        def slope = 1 / Math.tan(Math.toRadians(0.5 * sect_angle_coverage))
        def left = (height - slope * y_rot)
        def right = (height + slope * y_rot)
        return [left, right]
    }

    // local sector for electrons based on DCR2
    def getLocalSector(x,y,z){
	def phi = Math.toDegrees( Math.atan2( y/ (Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2) )), x/Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2) ) ) )
	
	if(phi < 30 && phi >= -30){        return 1}
	else if(phi < 90 && phi >= 30){    return 2}
	else if(phi < 150 && phi >= 90){   return 3}
	else if(phi >= 150 || phi < -150){ return 4}
	else if(phi < -90 && phi >= -150){ return 5}
	else if(phi < -30 && phi >= -90){  return 6}
	
	return 0    
    }


    def passElectronDCFiducialCutXY(event, index){
	def region = -1
	if( event.dc1_status.contains(index) ){
	    region=1
	}
	else if( event.dc2_status.contains(index) ){
	    region=2
	}
	else if( event.dc3_status.contains(index) ){
	    region=3
	}
	else{
	    // no region then return false
	    return false
	}

	def X = -10000
	def Y = -10000
	def hit = null
	switch (region){
	    case 1: 		
 		hit = event.dc1.get(index).find{ hh -> hh.layer == 6}
		X = hit.x
		Y = hit.y
		break
	    case 2:
		hit = event.dc1.get(index).find{ hh -> hh.layer == 18}
		X = hit.x
		Y = hit.y
		break
	    case 3:
		hit = event.dc1.get(index).find{ hh -> hh.layer == 36}
		X = hit.x
		Y = hit.y
		break
	    default:
		X=0
		Y=0
		break
	}
	
	def dc_sector = -1
	//determine sector in stupid way based on if sector 2 is hit
	if (event.dc2_status.contains(index) ){	    
	    def hitsR2 = event.dc2[index]
	    def hitR2 = hitsR2.find{it.layer==18}	
	    if( hitR2 ){
		dc_sector = getLocalSector(hitR2.x,hitR2.y,hitR2.z)
	    }
	}
	else{
	    return false
	}

	if(dc_sector == 2){
	    def X_new = X * Math.cos(-60 * Math.PI / 180) - Y * Math.sin(-60 * Math.PI / 180)
	    Y = X * Math.sin(-60 * Math.PI / 180) + Y * Math.cos(-60 * Math.PI / 180)
	    X = X_new
	}
	if (dc_sector == 3){
	    def X_new = X * Math.cos(-120 * Math.PI / 180) - Y * Math.sin(-120 * Math.PI / 180)
	    Y = X * Math.sin(-120 * Math.PI / 180) + Y * Math.cos(-120 * Math.PI / 180)
	    X = X_new
	}
	if (dc_sector == 4){
	    def X_new = X * Math.cos(-180 * Math.PI / 180) - Y * Math.sin(-180 * Math.PI / 180);
	    Y = X * Math.sin(-180 * Math.PI / 180) + Y * Math.cos(-180 * Math.PI / 180);
	    X = X_new;
	}
	
	if (dc_sector == 5){
	    def X_new = X * Math.cos(120 * Math.PI / 180) - Y * Math.sin(120 * Math.PI / 180);
	    Y = X * Math.sin(120 * Math.PI / 180) + Y * Math.cos(120 * Math.PI / 180);
	    X = X_new;
	}

	if (dc_sector == 6){
	    def X_new = X * Math.cos(60 * Math.PI / 180) - Y * Math.sin(60 * Math.PI / 180);
	    Y = X * Math.sin(60 * Math.PI / 180) + Y * Math.cos(60 * Math.PI / 180);
	    X = X_new;
	}
	region=region-1 // shift to match indexing for arrays
	def calc_min = dc_fiducial_cutXY_min[dc_sector - 1][region][0] + dc_fiducial_cutXY_min[dc_sector - 1][region][1] * X;
	def calc_max = dc_fiducial_cutXY_max[dc_sector - 1][region][0] + dc_fiducial_cutXY_max[dc_sector - 1][region][1] * X;
	
	return ((Y > calc_min) && (Y < calc_max))
	
    }
	
    // apply for electrons with outbending field and to hadrons for both inb and outb field
    // this is only for electrons here


    def passElectronDCFiducialCutThetaPhi(event, index){
	
	def region = -1
	if( event.dc1_status.contains(index) ){
	    region=1
	}
	else if( event.dc2_status.contains(index) ){
	    region=2
	}
	else if( event.dc3_status.contains(index) ){
	    region=3
	}
	else{
	    // no region then return false
	    return false
	}
	
	//def dc_fiducial_cutThetaPhi_min
	def theta_DCr = 5000;
	def phi_DCr_raw = 5000;
		
	switch (region){		
	    case 1:
		def hit = event.dc1.get(index).find{ hh -> hh.layer == 6}
		theta_DCr = (180/Math.PI) * Math.acos(hit.z / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) )
		phi_DCr_raw = (180/Math.PI) * Math.atan2( hit.y / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) , hit.x/Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2)))
		break;
		
	    case 2:
		def hit = event.dc2.get(index).find{ hh -> hh.layer == 18}
		theta_DCr = (180/Math.PI) * Math.acos(hit.z / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) ) 
		phi_DCr_raw = (180/Math.PI) * Math.atan2( hit.y / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) , hit.x/Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2)))
		break;
		
	    case 3:
		def hit = event.dc3.get(index).find{ hh -> hh.layer == 36}
		theta_DCr = (180/Math.PI) * Math.acos(hit.z / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) ) 
		phi_DCr_raw = (180/Math.PI) * Math.atan2( hit.y / Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2) ) , hit.x/Math.sqrt( Math.pow(hit.x,2) + Math.pow(hit.y,2) + Math.pow(hit.z,2)))
		break;
		
	    default: return false;
		break;
		
	}

	// rotate coordinates
	def phi_DCr = 5000;
	def sector = -1
	//determine sector in stupid way based on if sector 2 is hit
	if (event.dc2_status.contains(index) ){	    
	    def hitsR2 = event.dc2[index]
	    def hitR2 = hitsR2.find{it.layer==18}	
	    if( hitR2 ){
		sector = getLocalSector(hitR2.x,hitR2.y,hitR2.z)
	    }
	}
	else{
	    return false
	}
	
	if (sector == 1) phi_DCr = phi_DCr_raw;
	if (sector == 2) phi_DCr = phi_DCr_raw - 60;
	if (sector == 3) phi_DCr = phi_DCr_raw - 120;
	if (sector == 4 && phi_DCr_raw > 0) phi_DCr = phi_DCr_raw - 180;
	if (sector == 4 && phi_DCr_raw < 0) phi_DCr = phi_DCr_raw + 180;
	if (sector == 5) phi_DCr = phi_DCr_raw + 120;
	if (sector == 6) phi_DCr = phi_DCr_raw + 60;
	
	
	def pid = 0;
	
	switch (event.pid[index]){		
	    case 11: pid = 0;
		break;
	    default: 
		return false;
		break;		
	}
	//shift to match array syntax
	region = region - 1
	
	def calc_phi_min = dc_fiducial_cutThetaPhi_min[pid][sector - 1][region][0] + dc_fiducial_cutThetaPhi_min[pid][sector - 1][region][1] * Math.log(theta_DCr) 
        + dc_fiducial_cutThetaPhi_min[pid][sector - 1][region][2] * theta_DCr + dc_fiducial_cutThetaPhi_min[pid][sector - 1][region][3] * theta_DCr * theta_DCr;
	
	def calc_phi_max = dc_fiducial_cutThetaPhi_max[pid][sector - 1][region][0] + dc_fiducial_cutThetaPhi_max[pid][sector - 1][region][1] * Math.log(theta_DCr) 
        + dc_fiducial_cutThetaPhi_max[pid][sector - 1][region][2] * theta_DCr + dc_fiducial_cutThetaPhi_max[pid][sector - 1][region][3] * theta_DCr * theta_DCr;
		
	return ((phi_DCr > calc_phi_min) && (phi_DCr < calc_phi_max));
	
    }
    
    def passElectronDCFiducialCuts = { event, index ->
	if( mag_field == 0 ) // inbending, number faster than string comparison
	    return passElectronDCFiducialCutXY( event, index )
	else{
	    return passElectronDCFiducialCutThetaPhi( event, index )
	}
    }

    /*
     // old and dated do not use

    def passElectronDCR1 = { event, index ->
	if (event.dc1_status.contains(index)){
		def sec = event.dc_sector.get(index)-1
 		def hit = event.dc1.get(index).find{ hit -> hit.layer == 12}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR1[sec])		
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR1 )
            } else {
                return false
            }
	}
	return false
    }

    def passElectronDCR2 = { event, index ->
	if (event.dc2_status.contains(index)){
	    def sec = event.dc_sector[index]-1
	    def hit = event.dc2.get(index).find{ hit -> hit.layer == 24}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR2[sec])
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR2 )
            } else {
                return false
            }
	}
	return false
    }

    def passElectronDCR3 = { event, index ->
	if (event.dc3_status.contains(index)){
	    def sec = event.dc_sector[index]-1
	    def hit = event.dc3.get(index).find{ hit -> hit.layer == 36}
	    if (hit) {
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR3[sec])

                 // from lower line: && x1_rot**2 > radius2_DCr1){ return true }
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR3 )
            } else {
                return false
            }
	}
	return false
    }
     */
    

}
