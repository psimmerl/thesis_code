package mon
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.LinkedHashMap
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import utils.KinTool

class TimingPhiEventMon{

    def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def rfBucketLength = 4.0 //ns
    //def rf_offset = [5164:
    //5032: 
    
    ///////////////////////////////////////
    // proton, kp, km - beta cut parameters
    def prot_mean_p0 = [1.00236, 0.999663, 0.999094, 0.999417, 1.00021, 1.00034]
    def prot_mean_p1 = [0.933781, 0.919325, 0.908676, 0.919357, 0.90765, 0.908293]
    def prot_sigma_p0 = [0.0100589, 0.010606, 0.0101618, 0.0113726, 0.0098664, 0.00972719]
    def prot_sigma_p1 = [-0.0119585, -0.0125292, -0.0115819, -0.0140509, -0.0117452, -0.0109518]
    def prot_sigma_p2 = [0.0208288, 0.020747, 0.0202862, 0.0216464, 0.0207485, 0.0202303]

    def kp_mean_p0 = [1.00, 1.00, 1.00, 1.00, 1.00, 1.00]         // literature
    def kp_mean_p1= [0.244, 0.244, 0.244, 0.244, 0.244, 0.244]   // literature
    def kp_sigma_p0 =[0.00429074, 0.00515927, 0.0050134, 0.00562043, 0.00558013, 0.00524988]          // copied from pip
    def kp_sigma_p1 = [0.000723371, -0.000719322, -0.000348278, -0.00176441, -0.00157653, -0.00079676]  // copied from pip
    def kp_sigma_p2 = [0.00530764, 0.0059365, 0.00564377, 0.00663122, 0.00615767, 0.00606232]   // copied from pip
    
    def km_mean_p0 = [1.00, 1.00, 1.00, 1.00, 1.00, 1.00]        // literature
    def km_mean_p1 = [0.244, 0.244, 0.244, 0.244, 0.244, 0.244]  // literature
    def km_sigma_p0 = [0.0048138, 0.00520169, 0.00570919, 0.0062253, 0.00617771, 0.00578081]            // copied from pim
    def km_sigma_p1 = [-0.000382083, -0.000771668, -0.00179812, -0.00299382, -0.00282061, -0.00191396]  // copied from pim
    def km_sigma_p2 = [0.00617906, 0.00592537, 0.00671878, 0.00752452, 0.00707839, 0.00691364]          // copied from pim

    //call electron cut constructor
    def field_config = ''
    def electron = new ElectronFromEvent();

    def setCutParameters(){
	println(' [epkpkm::setCutParameters] -> field configuration is ' + field_config)
	println(' [epkpkm::setCutParameters] -> setting cut values for field ' + field_config)
	electron.setElectronCutParameters(field_config)
    }
    
    def everythingElsePassed(results, i) {
	for ( def j=0; j<results.size(); j++){	 
	    if (i != j && results[j] == false){		
		return false	
	    }	    
	}	
	return true	
    }
    
    def everythingElseFailed(results, i){	
	def reversed = results.collect{ !it }	
	return everythingElsePassed(reversed,i)	
    }

    def getVertexTime( time_at_detector, path_to_detector, calculated_beta ){
	def vertex_time = time_at_detector - path_to_detector/(calculated_beta * speed_of_light) 
	return vertex_time	  
    }
    
    def getStartTime( vertex_time, vertex_pos ){
	def vzCorr = (-target_position + vertex_pos)/speed_of_light
	//def deltatr = - vertex_time - vzCorr 
	//def rfCorr = 0//deltatr % rfBucketLength - rfBucketLength/2.0
	return vertex_time - vzCorr
    }
    
    def myElectronCutStrategies = [
	electron.passElectronChargeCut,
	//electron.passElectronMinMomentum,
	//electron.passElectronEIEOCut,
	//electron.passElectronNpheCut,
	//electron.passElectronVertexCut,
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3,
	//electron.passElectronSamplingFractionCut,
	electron.passElectronStatus,
	//electron.passElectronTrackQualityCut
	electron.passElectronEBPIDCut
    ]

    def charge_cut_index=0
    def min_mntm_cut_index=1
    def eieo_cut_index=2
    def nphe_cut_index=3
    def vz_cut_index=4
    def pcalfid_cut_index=5
    def dcr1_cut_index=6
    def dcr2_cut_index=7
    def dcr3_cut_index=8
    def ecsf_cut_index=9
    def status_cut_index=10
    def eb_el_pid_cut_index=11

    /*
     def trck_qual_cut_index=6
    def ecsf_cut_index=11
     */


    // Combinations of Cut Pass Configurations to Plot
    // charge, status and pass DCr1, r2, r3 fid cuts
    def pass_dc_comb=[true,true,true,true,false,false,false,true,true,true]
    // charge, status and pass PCAL fid DCr1, r2, r3 fid cuts
    def pass_dc_pcal_comb=[true,true,true,true,false,true,false,true,true,true]
    def pass_cut_comb = [0:pass_dc_comb,
			 1:pass_dc_pcal_comb]
    
    
    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	hcombs : {title-> new H1F("$title","$title",9,1,10)},
	hphim : {title -> new H1F("$title","$title",75,0.85,1.9)},
	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.5) },
	calsf : { title -> new H2F("$title","$title", 500, 0.0, 2.0, 500, 0.0, 2.0) },
	uvwsf : { title -> new H2F("$title","$title", 500, 0.0, 500.0, 500, 0.0, 0.3) },
	uvwhit : { title -> new H1F("$title","$title", 500, 0.0, 100.0) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 50, 0.0, 50 ) },
	nphe : { title -> new H1F("$title","$title", 40, 0.0, 40.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 2.0, 1000, 0.0, 2.0 ) },
	vz   : { title -> new H1F("$title","$title", 600, -40, 40) },
	vxy : { title -> new H2F("$title","$title", 1000, -1, 1, 1000, -1, 1) },

	time : { title -> new H1F("$title","$title", 2000, -450, 450) },
	time_small : { title -> new H1F("$title","$title",1000, -25, 25) },
	iron_cross : { title -> new H2F("$title","$title", 500, -20.0, 20.0, 500, -20.0, 20.0 ) },

	time2d_small : { title -> new H2F("$title","$title", 100, 0.0, 10.0, 500, -25, 25) },
	time2d_large : { title -> new H2F("$title","$title", 61, 0.0, 61.0, 2000, -450, 450) },
	betatime : { title -> new H2F("$title","$title", 200, 0.0, 1.20, 200, 50, 200) },
	path : { title -> new H1F("$title","$title", 2000, -400, 400) },//
	sectors : { title -> new H1F("$title","$title", 7, 0, 7) },
	helicity : { title -> new H1F("$title","$title", 6, -3 , 3) },
	w    : { title -> new H1F("$title","$title", 500, 0.5, beam+0.3 ) },
	q2w  : { title -> new H2F("$title","$title", 500, 0.5, beam+0.3, 500, 0, 2 * beam ) },  
	passrate : {title -> new H1F("$title","$title", 8, 0.0, 8.0 ) },
	wtheta : { title -> new H2F("$title","$title", 100, 0.0, 2.5, 100, 0.0, 20.0 ) },
	deltakin : {title -> new H1F("$title","$title", 500, -2, 2  ) },
	deltakin2D : {title -> new H2F("$title","$title", 1000, 0.0, 500, 500, -2, 2  ) },
	deltakin2D_beam : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) },
	hadron_multiplicity : {title -> new H1F("$title","$title", 10, 0.0, 10.0 ) },
	recpart_chi2pid : {title -> new H1F("$title","$title", 200, -5.0, 5.0 ) },
	recpart_chi2pid2D : {title -> new H2F("$title","$title", 200, 0.0,  beam+0.3, 200, -10.0, 10.0 ) },
	betap : { title -> new H2F("$title","$title", 200, 0.0, 8.0, 200, 0.0, 1.20 ) }

    ]

    def processEvent( event ){
	
	//println("[ElectronMon::processEvent] - processing event now")
	def run_number = event.run_number
	def my_el_cuts = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()
	//(0..<event.npart)?.findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	
	def el_list = []
	def el = my_el_cuts.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ index-> 
	    def el_lv = LorentzVector.withPID(11, event.px[index], event.py[index], event.pz[index])
	    def p = el_lv.p()
	    def vz = event.vz[index]
	    def theta = Math.toDegrees(el_lv.theta())
	    def phi = Math.toDegrees(el_lv.phi())
	    def good_timing_box_cut_ftof2 = null
	    def good_timing_corner_cut_ftof2 = false
	    def good_timing_box_cut_ftof1 = null
	    def good_timing_corner_cut_ftof1 = false
	    def good_timing_box_cut_ftof0 = null
	    def good_timing_corner_cut_ftof0 = false
	    def layer_to_use = -1

	    //find first instance of a proton present in an event
	    def hadron_in_event = (0..<event.npart).find{hh -> ( event.pid[hh] == 2212 && (event.status[hh] >= 2000 && event.status[hh] < 4000))}
	    
	    //println("[ElectronMon::processEvent] - at least a hadron (proton) in the event " + hadron_in_event)
	    	    
	    //println("[ElectronMon::processEvent] - Doing timing analysis now")
	    def el_trigger_htcc_vertex_time = null
	    def el_trigger_ftof2_vertex_time = null
	    def el_trigger_ftof1_vertex_time = null
	    def el_trigger_ftof0_vertex_time = null
	    def el_trigger_pcal_vertex_time = null

	    /////////////////////////////////////////
	    // compare htcc, ftof, and pcal timing
	    if( event.cherenkov_status.contains(index) && event.tof_status.contains(index) && event.pcal_status.contains(index) ){
		///////////////////////////
		//get htcc info
		def htcc_t = event.cherenkov_time[index]
		def htcc_p = event.cherenkov_path[index]/10.0
		def htcc_sector = event.cherenkov_sector[index]
		
		//////////////////////////
		//get ftof info
		//layer_to_use = -1
		//println(' layers available ' + event.tof[index].layer)
		def layer_avail = event.tof[index].layer
		def my_layer = []
		
		def ftof2_t = -1000
		def ftof2_p = -1000
		def ftof2_sector = -1

		def ftof1_t = -1000
		def ftof1_p = -1000
		def ftof1_sector = -1
		for( ii in layer_avail){
		    my_layer.add(Integer.toString(ii))
		}

		if( my_layer.contains('2') && !my_layer.contains('3') && !my_layer.contains('1') ){
		  //  println('use layer 2 for timing ' )
		    layer_to_use = 2
		}
		else if ( my_layer.contains('1') && !my_layer.contains('2') &&  !my_layer.contains('3')  ){
		    //println('use layer 1 for timing ' )
		    layer_to_use = 1
		}		    
		else if (my_layer.contains('1') && my_layer.contains('2') ){
		    layer_to_use = 0
		}

		//println(' layer to use '  + layer_to_use)
		//println(' layers available ' + layer_avail)

		event.tof[index].eachWithIndex{ hh, ii ->  
		    if(layer_to_use == 2 &&  hh.layer == 2 ){
 			//println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			ftof2_sector = hh.sector
			ftof2_t = hh.time
			ftof2_p = hh.path
		    }
		    if(layer_to_use == 1 && hh.layer == 1 ){
 			//println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			ftof1_sector = hh.sector
			ftof1_t = hh.time
			ftof1_p = hh.path
		    }

		    if( layer_to_use == 0 && hh.layer == 1 ){
			ftof1_sector = hh.sector
			ftof1_t = hh.time
			ftof1_p = hh.path
		    }		    
		    if( layer_to_use == 0 && hh.layer == 2 ){
			ftof2_sector = hh.sector
			ftof2_t = hh.time
			ftof2_p = hh.path
		    }		    
		}

		
		///////////////////////////
		//get pcal info
		def pcal_sect = event.pcal_sector[index]
		def pcal_t = event.pcal_time[index]
		def pcal_p = event.pcal_path[index]

		def measured_t_htcc = htcc_t - event.start_time // - event.vt[index]
		def delta_t_htcc = htcc_p/29.9792458 - measured_t_htcc
		def measured_t_pcal = pcal_t - event.start_time//event.vt[index]
		def delta_t_pcal = pcal_p/29.9792458 - measured_t_pcal

		def diff_t_htcc_pcal = delta_t_htcc - delta_t_pcal
		
		el_trigger_htcc_vertex_time = getStartTime(getVertexTime(htcc_t, htcc_p, 1), event.vz[index])
		el_trigger_pcal_vertex_time = getStartTime(getVertexTime(pcal_t, pcal_p, 1), event.vz[index])
		 		
 		if( layer_to_use == 2 &&  ftof2_sector != -1 && ftof2_t != -1000 && ftof2_sector == htcc_sector && ftof2_sector == pcal_sect ){
		    
		    def measured_t_ftof2 = ftof2_t - event.start_time //event.vt[index]
		    def delta_t_ftof2 = ftof2_p/29.9792458 - measured_t_ftof2
		    		    
		    def diff_t_ftof2_pcal = delta_t_ftof2 - delta_t_pcal
		    def diff_t_htcc_ftof2 = delta_t_htcc - delta_t_ftof2
		    def diff_t_ftof2_htcc = delta_t_ftof2 - delta_t_htcc

		    el_trigger_ftof2_vertex_time = getStartTime(getVertexTime(ftof2_t, ftof2_p, 1), event.vz[index])
		
		    
		    //println("PCAL TIME " + pcal_t + " PATH LENGTH " + pcal_p + " VERTEX POSITION " + event.vz[index])
		    //println("[ElectronMon::processEvent] - HTCC vertex time " + el_trigger_htcc_vertex_time )
		    //println("[ElectronMon::processEvent] - FTOFL2 vertex time " + el_trigger_ftof2_vertex_time )
		    //println("[ElectronMon::processEvent] - PCAL vertex time " + el_trigger_pcal_vertex_time )
		    //println("[ElectronMon::processEvent] - event.vt[index] " + event.vt[index])
		    
 		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof2_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_pcal)
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof2)
		    histos.computeIfAbsent("iron_cross_ftof2_minus_pcal_vs_ftof2_minus_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,diff_t_ftof2_htcc)
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_delta_t_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,delta_t_ftof2)
		    histos.computeIfAbsent("iron_cross_htcc_minus_ftof2_vs_delta_t_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_ftof2,delta_t_pcal)
		    histos.computeIfAbsent("iron_cross_ftof2_minus_pcal_vs_delta_t_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,delta_t_htcc)
		    
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof2_minus_pcal_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_pcal)
		
			    
		    //y - cut @ 3.9 +- 3*0.7
		    //x - cut @ -0.2459 +- 3*0.65
		    def min_ftof_pcal_cut = -0.2459 - 4.5*0.65
		    def max_ftof_pcal_cut = -0.2459 + 4.5*0.65

		    def min_htcc_ftof_cut = 3.9 - 4.5*0.7
		    def max_htcc_ftof_cut = 3.9 + 4.5*0.7
		    

 		    def n_protons_matched_box_cut=0
		    def n_protons_unmatched_box_cut=0
		    //def good_timing_box_cut = null
		    def n_protons_matched_corner_cut=0
		    def n_protons_unmatched_corner_cut=0

 		    if ( diff_t_ftof2_pcal < max_ftof_pcal_cut && diff_t_ftof2_pcal > min_ftof_pcal_cut && diff_t_htcc_ftof2 < max_htcc_ftof_cut && diff_t_htcc_ftof2 > min_htcc_ftof_cut ) {
			good_timing_box_cut_ftof2=true
		    }
		    else if( diff_t_ftof2_pcal < max_ftof_pcal_cut || diff_t_ftof2_pcal > min_ftof_pcal_cut || diff_t_htcc_ftof2 < max_htcc_ftof_cut || diff_t_htcc_ftof2 > min_htcc_ftof_cut ) {
			good_timing_box_cut_ftof2=false
		    }

		    //3.3, -0.37
		    if( (diff_t_ftof2_pcal > 3.3 && diff_t_htcc_ftof2 < -0.37) || ( diff_t_ftof2_pcal < -3.3 && diff_t_htcc_ftof2 > 7.1 ) ){
			//println('here')
			//println(my_el_cuts.size())
			//println(index)
			good_timing_corner_cut_ftof2=true			    
			//println(good_timing_corner_cut)
			//println(el_trigger_ftof_vertex_time)
		    }
		    else{
			good_timing_corner_cut_ftof2=false
		    }
	 	}
 		else if( layer_to_use == 1 &&  ftof1_sector != -1 && ftof1_t != -1000 && ftof1_sector == htcc_sector && ftof1_sector == pcal_sect ){

		    def measured_t_ftof1 = ftof1_t - event.start_time //event.vt[index]
		    def delta_t_ftof1 = ftof1_p/29.9792458 - measured_t_ftof1
		    
		    el_trigger_ftof1_vertex_time = getStartTime(getVertexTime(ftof1_t, ftof1_p, 1), event.vz[index])
		    
		    def diff_t_ftof1_pcal = delta_t_ftof1 - delta_t_pcal
		    def diff_t_htcc_ftof1 = delta_t_htcc - delta_t_ftof1
		    def diff_t_ftof1_htcc = delta_t_ftof1 - delta_t_htcc


		    
		    //println("PCAL TIME " + pcal_t + " PATH LENGTH " + pcal_p + " VERTEX POSITION " + event.vz[index])
		    //println("[ElectronMon::processEvent] - HTCC vertex time " + el_trigger_htcc_vertex_time )
		    //println("[ElectronMon::processEvent] - FTOFL1 vertex time " + el_trigger_ftof1_vertex_time )
		    //println("[ElectronMon::processEvent] - PCAL vertex time " + el_trigger_pcal_vertex_time )
		    //println("[ElectronMon::processEvent] - event.vt[index] " + event.vt[index])
		    

		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof1_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof1_pcal)
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof1)
		    histos.computeIfAbsent("iron_cross_ftof1_minus_pcal_vs_ftof1_minus_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,diff_t_ftof1_htcc)
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_delta_t_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,delta_t_ftof1)
		    histos.computeIfAbsent("iron_cross_htcc_minus_ftof1_vs_delta_t_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_ftof1,delta_t_pcal)
		    histos.computeIfAbsent("iron_cross_ftof1_minus_pcal_vs_delta_t_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,delta_t_htcc)
		    
		    histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof1_minus_pcal_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof1_pcal)


		    //y - cut @ 3.9 +- 3*0.7
		    //x - cut @ -0.2459 +- 3*0.65
		    def min_ftof_pcal_cut = -0.2459 - 4.5*0.65
		    def max_ftof_pcal_cut = -0.2459 + 4.5*0.65

		    def min_htcc_ftof_cut = 3.9 - 4.5*0.7
		    def max_htcc_ftof_cut = 3.9 + 4.5*0.7
		    

 		    def n_protons_matched_box_cut=0
		    def n_protons_unmatched_box_cut=0
		    //def good_timing_box_cut = null
		    def n_protons_matched_corner_cut=0
		    def n_protons_unmatched_corner_cut=0

 		    if ( diff_t_ftof1_pcal < max_ftof_pcal_cut && diff_t_ftof1_pcal > min_ftof_pcal_cut && diff_t_htcc_ftof1 < max_htcc_ftof_cut && diff_t_htcc_ftof1 > min_htcc_ftof_cut ) {
			good_timing_box_cut_ftof1=true
		    }
		    else if( diff_t_ftof1_pcal < max_ftof_pcal_cut || diff_t_ftof1_pcal > min_ftof_pcal_cut || diff_t_htcc_ftof1 < max_htcc_ftof_cut || diff_t_htcc_ftof1 > min_htcc_ftof_cut ) {
			good_timing_box_cut_ftof1=false
		    }

		    //3.3, -0.37
		    if( (diff_t_ftof1_pcal > 3.3 && diff_t_htcc_ftof1 < -0.37) || ( diff_t_ftof1_pcal < -3.3 && diff_t_htcc_ftof1 > 7.1 ) ){
			//println('here')
			//println(my_el_cuts.size())
			//println(index)
			good_timing_corner_cut_ftof1=true			    
			//println(good_timing_corner_cut)
			//println(el_trigger_ftof_vertex_time)
		    }
		    else{
			good_timing_corner_cut_ftof1=false
		    }

		}
		else if( layer_to_use == 0 &&  ftof1_sector != -1 && ftof1_t != -1000 && ftof1_sector == htcc_sector && ftof1_sector == pcal_sect && ftof2_sector != -1 && ftof2_t != -1000 && ftof2_sector == htcc_sector && ftof2_sector == pcal_sect){

		    def measured_t_ftof1 = ftof1_t - event.start_time //event.vt[index]
		    def delta_t_ftof1 = ftof1_p/29.9792458 - measured_t_ftof1
		    def measured_t_ftof2 = ftof2_t - event.start_time //event.vt[index]
		    def delta_t_ftof2 = ftof2_p/29.9792458 - measured_t_ftof2
		    
		    el_trigger_ftof1_vertex_time = getStartTime(getVertexTime(ftof1_t, ftof1_p, 1), event.vz[index])
		    el_trigger_ftof2_vertex_time = getStartTime(getVertexTime(ftof2_t, ftof2_p, 1), event.vz[index])
		    
		    def diff_t_ftof1_pcal = delta_t_ftof1 - delta_t_pcal
		    def diff_t_htcc_ftof1 = delta_t_htcc - delta_t_ftof1
		    def diff_t_ftof1_htcc = delta_t_ftof1 - delta_t_htcc
		   		    
		    def diff_t_ftof2_pcal = delta_t_ftof2 - delta_t_pcal
		    def diff_t_htcc_ftof2 = delta_t_htcc - delta_t_ftof2
		    def diff_t_ftof2_htcc = delta_t_ftof2 - delta_t_htcc

		    def diff_t_ftof2_ftof1 = delta_t_ftof2 - delta_t_ftof1

		
		    //println(" BOTH FTOF LAYERS HIT " )
		    //println("PCAL TIME " + pcal_t + " PATH LENGTH " + pcal_p + " VERTEX POSITION " + event.vz[index])
		    //println("[ElectronMon::processEvent] - HTCC vertex time " + el_trigger_htcc_vertex_time )
		    //println("[ElectronMon::processEvent] - FTOFL2 vertex time " + el_trigger_ftof2_vertex_time )
		    //println("[ElectronMon::processEvent] - FTOFL1 vertex time " + el_trigger_ftof1_vertex_time )
		    //println("[ElectronMon::processEvent] - PCAL vertex time " + el_trigger_pcal_vertex_time )
		    //println("[ElectronMon::processEvent] - event.vt[index] " + event.vt[index])
		    

 		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof2_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_pcal)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_htcc_minus_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof2)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_pcal_vs_ftof2_minus_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,diff_t_ftof2_htcc)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_delta_t_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,delta_t_ftof2)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_ftof2_vs_delta_t_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_ftof2,delta_t_pcal)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_pcal_vs_delta_t_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,delta_t_htcc)		    
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof2_minus_pcal_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_pcal)

		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof1_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof1_pcal)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_htcc_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_pcal_vs_ftof1_minus_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,diff_t_ftof1_htcc)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_delta_t_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,delta_t_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_ftof1_vs_delta_t_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_ftof1,delta_t_pcal)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_pcal_vs_delta_t_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,delta_t_htcc)		    
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof1_minus_pcal_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof1_pcal)

		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof2_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof2_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_ftof1)
		    //minus sign to keep procedure uniform
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_pcal_vs_ftof1_minus_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,-diff_t_ftof2_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_pcal_vs_ftof2_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,diff_t_ftof2_ftof1)
		    //minus sign to keep procedure uniform
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_htcc_vs_ftof1_minus_ftof2_passall",histoBuilders.iron_cross).fill(diff_t_ftof1_htcc,-diff_t_ftof2_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_htcc_vs_ftof2_minus_ftof1_passall",histoBuilders.iron_cross).fill(diff_t_ftof2_htcc,diff_t_ftof2_ftof1)
		    histos.computeIfAbsent("iron_cross_bothftofs_htcc_minus_pcal_vs_ftof2_minus_ftof1_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof2_ftof1)


		    if( hadron_in_event ){
			histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_pcal_vs_ftof2_minus_ftof1_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,diff_t_ftof2_ftof1)
			histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_htcc_vs_ftof2_minus_ftof1_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof2_htcc,diff_t_ftof2_ftof1)
			histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_htcc_vs_ftof1_minus_ftof2_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof1_htcc,-diff_t_ftof2_ftof1)
			histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_pcal_vs_ftof1_minus_ftof2_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,-diff_t_ftof2_ftof1)
			histos.computeIfAbsent("iron_cross_bothftofs_ftof1_minus_pcal_vs_ftof1_minus_htcc_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof1_pcal,diff_t_ftof1_htcc)
			histos.computeIfAbsent("iron_cross_bothftofs_ftof2_minus_pcal_vs_ftof2_minus_htcc_passall_at_least_1proton",histoBuilders.iron_cross).fill(diff_t_ftof2_pcal,diff_t_ftof2_htcc)
		    }
		    //y - cut @ 3.9 +- 3*0.7
		    //x - cut @ -0.2459 +- 3*0.65
		    def min_ftof_pcal_cut = -0.2459 - 4.5*0.65
		    def max_ftof_pcal_cut = -0.2459 + 4.5*0.65

		    def min_htcc_ftof_cut = 3.9 - 4.5*0.7
		    def max_htcc_ftof_cut = 3.9 + 4.5*0.7
		    

 		    def n_protons_matched_box_cut=0
		    def n_protons_unmatched_box_cut=0
		    //def good_timing_box_cut = null
		    def n_protons_matched_corner_cut=0
		    def n_protons_unmatched_corner_cut=0

 		    if ( diff_t_ftof2_pcal < max_ftof_pcal_cut && diff_t_ftof2_pcal > min_ftof_pcal_cut && diff_t_htcc_ftof2 < max_htcc_ftof_cut && diff_t_htcc_ftof2 > min_htcc_ftof_cut
			&& diff_t_ftof1_pcal < max_ftof_pcal_cut && diff_t_ftof1_pcal > min_ftof_pcal_cut && diff_t_htcc_ftof1 < max_htcc_ftof_cut && diff_t_htcc_ftof1 > min_htcc_ftof_cut) {
			good_timing_box_cut_ftof0=true
		    }
		    else if( diff_t_ftof2_pcal < max_ftof_pcal_cut || diff_t_ftof2_pcal > min_ftof_pcal_cut || diff_t_htcc_ftof2 < max_htcc_ftof_cut || diff_t_htcc_ftof2 > min_htcc_ftof_cut 
		    && diff_t_ftof1_pcal < max_ftof_pcal_cut || diff_t_ftof1_pcal > min_ftof_pcal_cut || diff_t_htcc_ftof1 < max_htcc_ftof_cut || diff_t_htcc_ftof1 > min_htcc_ftof_cut ) {
			good_timing_box_cut_ftof0=false
		    }

		    //3.3, -0.37
		    /*
		    if( (diff_t_ftof2_pcal > 3.3 && diff_t_htcc_ftof2 < -0.37) || ( diff_t_ftof2_pcal < -3.3 && diff_t_htcc_ftof2 > 7.1 )
		       && (diff_t_ftof1_pcal > 3.3 && diff_t_htcc_ftof1 < -0.37) || ( diff_t_ftof1_pcal < -3.3 && diff_t_htcc_ftof1 > 7.1 )){
			//println('here')
			//println(my_el_cuts.size())
			//println(index)
			good_timing_corner_cut_ftof0=true			    
			println(good_timing_corner_cut)
			//println(el_trigger_ftof_vertex_time)
		    }
		    else{
			good_timing_corner_cut_ftof0=false
		    }
		     */
		    //iron_cross_bothftofs_ftof2_minus_pcal_vs_ftof2_minus_ftof1_passall
		    //diff_t_ftof2_pcal,diff_t_ftof2_ftof1
		    //corner cut on ftof2 - ftof1 when both layers are hit		    
		    if( ( diff_t_ftof2_pcal > 0.839552 && diff_t_ftof2_ftof1  > 1.16099 )  || ( diff_t_ftof2_pcal < -1.58582 && diff_t_ftof2_ftof1 < -1.31579 ) ){
			//println('here')
			//println(my_el_cuts.size())
			//println(index)
			good_timing_corner_cut_ftof0=true			    
			//println(good_timing_corner_cut_ftof0)
			//println(el_trigger_ftof_vertex_time)
		    }
		    else{
			good_timing_corner_cut_ftof0=false
		    }


		}


		
			
		//println("[ElectronMon::processEvent] - HTCC and FTOF sector are same")
		//println("[ElectronMon::processEvent] - counting number of protons in event for when timing is matched")			

	    }
	
	    ///println(good_timing_corner_cut)
	    def m_ele = [eindex:index, 
			 corner_time_result2: good_timing_corner_cut_ftof2, 
			 corner_time_result1: good_timing_corner_cut_ftof1, 
			 corner_time_result0: good_timing_corner_cut_ftof0, 
			 ftof_layer:layer_to_use, 
			 lv: el_lv, 
			 htcctime: el_trigger_htcc_vertex_time, 
			 ftof0time: el_trigger_ftof0_vertex_time, 
			 ftof1time: el_trigger_ftof1_vertex_time, 
			 ftof2time: el_trigger_ftof2_vertex_time, 
			 pcaltime:el_trigger_pcal_vertex_time ]
	    el_list.add(m_ele)
	}
	
	

	//hadron map for pid
	def pr_map = []
	def kp_map = []
	def km_map = []
	def hadron_list = []


	if( el_list.size() > 0  ){
	    el_list.each{ ele -> 
   		//println("[ElectronMon::processEvent] - electron passed corner time result " )
		(0..<event.npart).findAll{ ( (event.charge[it] != 0)) && (it!=ele['eindex']) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{hindx->
		    def hadron_pid = event.pid[hindx]		
		    def hadron_charge = event.charge[hindx]
		    def rec_hadron_chi2 = event.chi2pid[hindx]
		    def hadron_p = event.p[hindx]
		    def hadron_vz = event.vz[hindx]
		    def hadron_vt = event.vt[hindx]
		    def hadron_beta = event.beta[hindx]
		    def hadron_status = event.status[hindx]
		    def theory_beta = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.938*0.938)
		    
		    def hadron_vertex_corr_htcc_start_time = 0// getStartTime(ele['htcctime'], hadron_vz)
 		    def hadron_vertex_corr_ftof2_start_time = 0//getStartTime(ele['ftoftime'], hadron_vz)
 		    def hadron_vertex_corr_ftof1_start_time = 0//getStartTime(ele['ftoftime'], hadron_vz)
 		    def hadron_vertex_corr_ftof0_start_time = 0//getStartTime(ele['ftoftime'], hadron_vz)
 		    def hadron_vertex_corr_pcal_start_time = 0//getStartTime(ele['pcaltime'], hadron_vz)
		    
		    if ( Math.abs(rec_hadron_chi2) < 4.2 ){
			def pr_ftof_t = -1000
			def pr_ftof_p = -1000
			def pr_ftof_sector = -1
			def pr_layer_to_use=null
			if( event.tof_status.contains(hindx) ){
			    //////////////////////////
			    //get ftof info
			    def layer_to_use = -1
			    //println(' layers available ' + event.tof[index].layer)
			    def pr_layer_avail = event.tof[hindx].layer
			    def pr_my_layer = []
			    
			    for( ii in pr_layer_avail){
				pr_my_layer.add(Integer.toString(ii))
			    }
			    if( pr_my_layer.contains('2') && !pr_my_layer.contains('3') && !pr_my_layer.contains('1') ){
				//println('use layer 2 for timing ' )
				layer_to_use = 2
			    }
			    else if ( !pr_my_layer.contains('2') && pr_my_layer.contains('1') && !pr_my_layer.contains('3')  ){
				//println('use layer 1 for timing ' )
				layer_to_use = 1
			    }
			    else if ( pr_my_layer.contains('1') && pr_my_layer.contains('2') && !pr_my_layer.contains('3')  ){
				//println('use layer 1 for timing ' )
				layer_to_use = 0
			    }		    
			    
			    event.tof[hindx].eachWithIndex{ hh, ii ->  
				if( layer_to_use == 2 ){//&& hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 				    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
				    pr_ftof_sector = hh.sector
				    pr_ftof_t = hh.time
				    pr_ftof_p = hh.path	    
				}
				if(layer_to_use == 1){// && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 				    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
				    pr_ftof_sector = hh.sector
 				    pr_ftof_t = hh.time
				    pr_ftof_p = hh.path	    
				}
				if( layer_to_use == 0 && hh.layer == 2 ){
				    pr_ftof_sector = hh.sector
 				    pr_ftof_t = hh.time
				    pr_ftof_p = hh.path	    
				}
			    }
			    pr_layer_to_use = layer_to_use				
			}
			
			// counts the number of protons for when electron passes cut


			def ftof_sector = pr_ftof_sector-1
			//	    println('pr tof ' + pr_ftof_t)


			if( ele['htcctime'] != null ){ hadron_vertex_corr_htcc_start_time = getStartTime(ele['htcctime'], hadron_vz) }
			if( ele['ftof2time'] != null ) { hadron_vertex_corr_ftof2_start_time = getStartTime(ele['ftof2time'], hadron_vz) }
   			if( ele['ftof1time'] != null ) { hadron_vertex_corr_ftof1_start_time = getStartTime(ele['ftof1time'], hadron_vz) }
			if( ele['ftof0time'] != null ) { hadron_vertex_corr_ftof0_start_time = getStartTime(ele['ftof0time'], hadron_vz) }
			if( ele['pcaltime'] != null ){ hadron_vertex_corr_pcal_start_time = getStartTime(ele['pcaltime'], hadron_vz) }
 			if( (pr_layer_to_use == 0 || pr_layer_to_use == 1 || pr_layer_to_use == 2 ) && pr_ftof_t > 0 && ftof_sector >= 0 ){ 			       			      	
			    
			    
 			    if( hadron_pid == 2212 ){
				histos.computeIfAbsent("charged_beta_p_protons_ftofl$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_beta)
 				histos.computeIfAbsent("charged_beta_p_protons",histoBuilders.betap).fill(hadron_p,hadron_beta)
			    }

			    def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
			    def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
			    def hadron_ftof2_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof2_start_time)/speed_of_light
			    def hadron_ftof1_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof1_start_time)/speed_of_light
			    //def hadron_ftof1_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof1_start_time)/speed_of_light
			    def hadron_ftof0_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof0_start_time)/speed_of_light
			    
			    def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
			    def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light // dont use 
			
			
			    //				println("[ElectronMon::processEvent] - calculated beta with vt " +  calc_rec_beta)
			    //
			    ///////////////////////
			    // Determine hadron pid
			    //println(ftof_sector)
			    def pr_beta_p_mean = prot_mean_p0[ftof_sector]*hadron_p/Math.sqrt( hadron_p*hadron_p + prot_mean_p1[ftof_sector])
			    def pr_beta_p_sig = prot_sigma_p0[ftof_sector] + prot_sigma_p1[ftof_sector]/Math.sqrt(hadron_p) + prot_sigma_p2[ftof_sector]/(hadron_p*hadron_p)
			    
			    def kp_beta_p_mean = kp_mean_p0[ftof_sector] * hadron_p/Math.sqrt(hadron_p*hadron_p + kp_mean_p1[ftof_sector])
			    def kp_beta_p_sig = kp_sigma_p0[ftof_sector] + kp_sigma_p1[ftof_sector]/Math.sqrt(hadron_p) + kp_sigma_p2[ftof_sector]/(hadron_p*hadron_p)
			    
			    def km_beta_p_mean = km_mean_p0[ftof_sector] * hadron_p/Math.sqrt(hadron_p*hadron_p + km_mean_p1[ftof_sector])
			    def km_beta_p_sig = km_sigma_p0[ftof_sector] + km_sigma_p1[ftof_sector]/Math.sqrt(hadron_p) + km_sigma_p2[ftof_sector]/(hadron_p*hadron_p)
			    
			    def hadron_type = -10000
			    def pr_present=false
			    def kp_present=false
			    def km_present=false

			    def ns = 4.5
			    def nns = 4.5

			    
			    // t0 - [tftof - l/beta_calc]
			    def beta_pr_cal = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.938*0.938)
			    def beta_kp_cal = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.493*0.493)
 			    def beta_pip_cal = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.139*0.139)
			    def beta_deut_cal = hadron_p/Math.sqrt(hadron_p*hadron_p + 1.8756*1.8756)

			    def pr_delta_t_ftof1 = hadron_vertex_corr_ftof1_start_time - (pr_ftof_t - pr_ftof_p/(beta_pr_cal*speed_of_light))
			    def kp_delta_t_ftof1 = hadron_vertex_corr_ftof1_start_time - (pr_ftof_t - pr_ftof_p/(beta_kp_cal*speed_of_light))
			    def pip_delta_t_ftof1 = hadron_vertex_corr_ftof1_start_time - (pr_ftof_t - pr_ftof_p/(beta_pip_cal*speed_of_light))
			    def deut_delta_t_ftof1 = hadron_vertex_corr_ftof1_start_time - (pr_ftof_t - pr_ftof_p/(beta_deut_cal*speed_of_light))
			    
			    
			    def pr_delta_t_ftof2 = hadron_vertex_corr_ftof2_start_time - (pr_ftof_t - pr_ftof_p/(beta_pr_cal*speed_of_light))
			    def kp_delta_t_ftof2 = hadron_vertex_corr_ftof2_start_time - (pr_ftof_t - pr_ftof_p/(beta_kp_cal*speed_of_light))
			    def pip_delta_t_ftof2 = hadron_vertex_corr_ftof2_start_time - (pr_ftof_t - pr_ftof_p/(beta_pip_cal*speed_of_light))
			    def deut_delta_t_ftof2 = hadron_vertex_corr_ftof2_start_time - (pr_ftof_t - pr_ftof_p/(beta_deut_cal*speed_of_light))

 			    def timing_list_ftof1 = [Math.abs(pr_delta_t_ftof1), Math.abs(kp_delta_t_ftof1), Math.abs(pip_delta_t_ftof1), Math.abs(deut_delta_t_ftof1)]
 			    def timing_list_ftof2 = [Math.abs(pr_delta_t_ftof2), Math.abs(kp_delta_t_ftof2), Math.abs(pip_delta_t_ftof2), Math.abs(deut_delta_t_ftof2)]

 			    int indexOfMinimum_ftof1 = timing_list_ftof1.indexOf(Collections.min(timing_list_ftof1));
 			    int indexOfMinimum_ftof2 = timing_list_ftof2.indexOf(Collections.min(timing_list_ftof2));


			    if( hadron_charge > 0 && indexOfMinimum_ftof2 == 0 &&  ele['ftof2time'] != null ){
				histos.computeIfAbsent("proton_beta_p_redefined_ftof2_st_test",histoBuilders.betap).fill(hadron_p,hadron_ftof2_beta)				
				histos.computeIfAbsent("proton_beta_eventbuilder_ftof2_to_compare_against",histoBuilders.betap).fill(hadron_p,hadron_beta)

			    }
			    
			    if( ele['corner_time_result0'] == true ){			    					

				if( ele['ftof1time'] != null) histos.computeIfAbsent("charged_beta_p_redfined_ftof1_st_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
				if( ele['ftof1time'] != null) histos.computeIfAbsent("charged_beta_p_redfined_ftof1_st_allhadronftofs",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
				
				
				histos.computeIfAbsent("charged_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_beta)
 				if( ele['pcaltime'] != null ) histos.computeIfAbsent("charged_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
				if( ele['htcctime'] != null ) histos.computeIfAbsent("charged_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 				if( ele['ftof2time'] != null) histos.computeIfAbsent("charged_beta_p_redfined_ftof2_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_ftof2_beta)
 				if( ele['ftof1time'] != null) histos.computeIfAbsent("charged_beta_p_redfined_ftof1_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
				histos.computeIfAbsent("charged_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
				histos.computeIfAbsent("charged_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched_hadronAtFtof_$pr_layer_to_use",histoBuilders.betap).fill(hadron_p,calc_rec_beta)					
				
			
				if (hadron_charge > 0 ){
				    if( indexOfMinimum_ftof1 == 0 ){
					//draw proton beta vs p 
					histos.computeIfAbsent("proton_beta_redefined_ftof1_st_mybeta",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					//check the case when my pid matches EB pid
					if( hadron_pid == 2212 ){
 					    histos.computeIfAbsent("proton_beta_redefined_ftof1_st_mybeta_match_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					}
					else{
					    //check the case when my pid does not match EB pid
					    histos.computeIfAbsent("proton_beta_redefined_ftof1_st_mybeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					    histos.computeIfAbsent("proton_beta_redefined_ftof1_st_ebbeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_beta)
					}				    
				    }
				    if( indexOfMinimum_ftof1 == 1 ){
					histos.computeIfAbsent("kp_beta_redefined_ftof1_st_mybeta",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					//draw kp beta vs p 
					//check the case when my pid matches EB pid
					if( hadron_pid == 321 ){
 					    histos.computeIfAbsent("kp_beta_redefined_ftof1_st_mybeta_match_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					}
					else{
					    //check the case when my pid does not match EB pid
					    histos.computeIfAbsent("kp_beta_redefined_ftof1_st_mybeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					    histos.computeIfAbsent("kp_beta_redefined_ftof1_st_ebbeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_beta)
					}				    
				    }
				    if( indexOfMinimum_ftof1 == 2 ){
					histos.computeIfAbsent("pip_beta_redefined_ftof1_st_mybeta",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					//check the case when my pid matches EB pid
					if( hadron_pid == 211 ){
 					    histos.computeIfAbsent("pip_beta_redefined_ftof1_st_mybeta_match_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					}
					else{
					    //check the case when my pid does not match EB pid
					    histos.computeIfAbsent("pip_beta_redefined_ftof1_st_mybeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					    histos.computeIfAbsent("pip_beta_redefined_ftof1_st_ebbeta_notmatch_eb",histoBuilders.betap).fill(hadron_p,hadron_beta)
					}				    
					//draw pip beta vs p 
				    }				    
				}				
				//println(indexOfMinimum)

				if( event.charge[hindx] > 0 ){
				    if(hadron_p > 4 ) { ns = 2.5 }
				    def pr_betap_upper_lim = pr_beta_p_mean + ns*pr_beta_p_sig
				    def pr_betap_lower_lim = pr_beta_p_mean - ns*pr_beta_p_sig
				    
				    def kp_betap_upper_lim = kp_beta_p_mean + ns*kp_beta_p_sig
				    def kp_betap_lower_lim = kp_beta_p_mean - ns*kp_beta_p_sig
				    
 				    // println("[ElectronMon::processEvent] - rec beta " + calc_rec_beta )
 				    //println("[ElectronMon::processEvent] - pr beta min " + pr_betap_lower_lim + " pr beta max " + pr_betap_upper_lim )
 				    //println("[ElectronMon::processEvent] - kp beta min " + kp_betap_lower_lim + " kp beta max " + kp_betap_upper_lim )

				    if( calc_rec_beta < pr_betap_upper_lim && calc_rec_beta > pr_betap_lower_lim){
					hadron_type=2212
					pr_present=true
				    }				   
				    else if( calc_rec_beta < kp_betap_upper_lim && calc_rec_beta > kp_betap_lower_lim){
					hadron_type=321
					kp_present=true
				    }
				    //println("status of proton " + pr_present + " status of kaon " + kp_present)
				    if( hadron_type == 2212 && pr_present == true && kp_present==false){
					pr_map=[pid: 2212, indx:hindx, lv: LorentzVector.withPID(2212, event.px[hindx], event.py[hindx], event.pz[hindx]), beta:calc_rec_beta]
					hadron_list.add(pr_map)
					//	println(hadron_list)

				    }
				    
				    if( hadron_type == 321 && kp_present == true && pr_present==false){
					kp_map=[pid: 321, indx:hindx, lv: LorentzVector.withPID(321, event.px[hindx], event.py[hindx], event.pz[hindx]), beta:calc_rec_beta]
					hadron_list.add(kp_map)
					//	println(hadron_list)

				    }				    
				}
				else if( event.charge[hindx] < 0 ){
				    if(hadron_p > 4 ) { nns = 2.5 }
				    def km_betap_upper_lim = km_beta_p_mean + nns*km_beta_p_sig
				    def km_betap_lower_lim = km_beta_p_mean - nns*km_beta_p_sig
				    if( calc_rec_beta < km_betap_upper_lim && calc_rec_beta > km_betap_lower_lim){
					hadron_type=-321					
					km_map=[pid: -321, indx:hindx, lv: LorentzVector.withPID(321, event.px[hindx], event.py[hindx], event.pz[hindx]), beta:calc_rec_beta]
					hadron_list.add(km_map)
					histos.computeIfAbsent("km_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_beta)
 					histos.computeIfAbsent("km_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
					histos.computeIfAbsent("km_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 					histos.computeIfAbsent("km_beta_p_redfined_ftof2_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_ftof2_beta)
 					histos.computeIfAbsent("km_beta_p_redfined_ftof1_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_ftof1_beta)
					histos.computeIfAbsent("km_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
					histos.computeIfAbsent("km_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,calc_rec_beta)				
				    }
				}
				
			    }
		    	    
			}
		    }


		    ///////////loop over to make combos
		    def combs = el_list.collectMany{ mele->
			(0..<hadron_list.size).findAll{ hh -> hadron_list[hh]['pid'] == 2212 }.findResults{ipr->
			    //println('got proton, index ' +ipr)
			    histos.computeIfAbsent("pr_betap",histoBuilders.betap).fill(hadron_list[ipr]['lv'].p(), hadron_list[ipr]['beta'])
			    return [mele, hadron_list.get(ipr)]
			}
		    }.collectMany{mele,mpro->
			(0..<hadron_list.size).findAll{ hh -> hadron_list[hh]['pid'] == 321 }.findResults{ikp->
			    //println('got kp, index ' +ikp)
			    histos.computeIfAbsent("kp_betap",histoBuilders.betap).fill(hadron_list[ikp]['lv'].p(), hadron_list[ikp]['beta'])
			    return [mele, mpro, hadron_list.get(ikp)]
			}
		    }.collectMany{mele,mpro,mkp->
			(0..<hadron_list.size).findAll{ hh -> hadron_list[hh]['pid'] == -321 }.findResults{ikm->
			    //println('got km, index ' +ikm)
			    histos.computeIfAbsent("km_betap",histoBuilders.betap).fill(hadron_list[ikm]['lv'].p(), hadron_list[ikm]['beta'])
			    return [mele, mpro, mkp, hadron_list.get(ikm)]
			}
		    }.findResults{elec,pr,kp,km->
			return [elec,pr,kp,km]
		    }
		    
		    //println("[TimingPhiEventMon::processEvent] - number of combinations " + combs.size())
		    histos.computeIfAbsent('combsize',histoBuilders.hcombs).fill(combs.size())
		    for(comb in combs){
			def (m_elec,m_pro,m_kp,m_km) = comb
			//println(" kp " + m_kp['lv'] + " km " + m_km['lv'] )
			def electron = m_elec['lv']
			def pro = m_pro['lv']
			def kp = m_kp['lv']
			def km = m_km['lv']
			
			def vphi = kp+km
			//println('[TimingPhiEventMon::processEvent] - phi mass ' + vphi.mass())
			histos.computeIfAbsent('phim',histoBuilders.hphim).fill(vphi.mass())	
			histos.computeIfAbsent("pr_betap_final",histoBuilders.betap).fill(m_pro['lv'].p(), m_pro['beta'])
			histos.computeIfAbsent("kp_betap_final",histoBuilders.betap).fill(m_kp['lv'].p(), m_kp['beta'])
			histos.computeIfAbsent("km_betap_final",histoBuilders.betap).fill(m_km['lv'].p(), m_km['beta'])

			
		    }
		    
		}
		//return hadron_list
	    }	    


	}

	
	
	//println( hadron_list )
	//println(hadron_list.size())
	/////////// DONE REDEFINING EVENTS ///////////////
	

	//println("[ElectronMon::processEvent] - Number of matched HTCC/FTOF protons " + n_protons_matched)
	//println("[ElectronMon::processEvent] - Number of UNmatched HTCC/FTOF protons " + n_protons_unmatched)
	/*
	 if(good_timing_box_cut && n_protons_matched_box_cut >= 0 ){
	 //add extra 1 bc started at -1
	 histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_box_cut_matched",histoBuilders.hadron_multiplicity).fill(n_protons_matched_box_cut)
	 histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_box_cut_matched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
 	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_matched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_matched",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall_box_cut_matched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_ftof_htcc)

     }
	 if( !good_timing_box_cut &&  n_protons_unmatched_box_cut >= 0 ){
	 histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_box_cut_unmatched",histoBuilders.hadron_multiplicity).fill(n_protons_unmatched_box_cut)
	 histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_box_cut_unmatched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall_box_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_ftof_htcc)
     }
	 
	 if(good_timing_corner_cut && n_protons_matched_corner_cut >= 0 ){
	 //add extra 1 bc started at -1
	 histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_corner_cut_matched",histoBuilders.hadron_multiplicity).fill(n_protons_matched_corner_cut)
	 histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_corner_cut_matched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_matched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_matched",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_ftof_minus_htcc_passall_corner_cut_matched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_ftof_htcc)
     }
	 if( !good_timing_corner_cut &&  n_protons_unmatched_corner_cut >= 0 ){
	 histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_corner_cut_unmatched",histoBuilders.hadron_multiplicity).fill(n_protons_unmatched_corner_cut)
	 histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_corner_cut_unmatched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof)
	 histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall_corner_cut_unmatched",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_ftof_htcc)

     }
	 */
    	
	
	
	
    }
}
