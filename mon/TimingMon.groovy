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

class TimingMon{

    def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def rfBucketLength = 4.0 //ns
    //def rf_offset = [5164:
    //5032: 
    
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
	def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	//println(my_el_cuts)
	
	my_el_cuts.each{ index, value ->
	    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	    def p = lv.mag()
	    def vz = event.vz[index]
	    def theta = Math.toDegrees(lv.theta())
	    def phi = Math.toDegrees(lv.phi())
	    
	    value.eachWithIndex{ cut_res, cut_ind -> cut_res ? histos.computeIfAbsent("passResults",histoBuilders.passrate).fill(cut_ind) : null }
	    
	    //println("[ElectronMon::processEvent] - Doing timing analysis now")
	    //println("[ElectronMon::processEvent] - value result and index ")
	    //println(index + " " + value)
	    if(!value.contains(false)){
		if( event.tof_status.contains(index) && event.pcal_status.contains(index) ){
		    def pcal_sect = event.pcal_sector[index]
		    def pcal_t = event.pcal_time[index]
		    def pcal_p = event.pcal_path[index]
		    
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[index].layer
		    def my_layer = []
		    
		    def ftof_t = -1000
		    def ftof_p = -1000
		    def ftof_sector = -1
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }
		    if( my_layer.contains('2') && !my_layer.contains('3') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    
		    event.tof[index].eachWithIndex{ hh, ii ->  
			if(layer_to_use!=-1 && 2 == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 			    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			    ftof_sector = hh.sector
			    ftof_t = hh.time
			    ftof_p = hh.path
			}
		    }

		    if( ftof_sector != -1 && ftof_t != -1000 && ftof_sector == pcal_sect){
			def delta_t_scint_pcal = pcal_t - ftof_t
			def measured_t_scint  = ftof_t - event.vt[index]
			def measured_t_pcal  = pcal_t - event.vt[index]
			def delta_t_pcal  = pcal_p/29.9792458 - measured_t_pcal
			def delta_t_scint  = ftof_p/29.9792458 - measured_t_scint
			histos.computeIfAbsent("delta_t_scint_pcal_s$ftof_sector"+"_passall",histoBuilders.time_small).fill(delta_t_scint_pcal)
			histos.computeIfAbsent("el_p_vs_delta_t_scint_pcal_s$ftof_sector"+"_passall",histoBuilders.time2d_small).fill(event.p[index], delta_t_scint_pcal)
			histos.computeIfAbsent("iron_cross_xscint_ypcal_passall",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xscint_ypcal_r$run_number"+"_passall",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xscint_ypcal_s$ftof_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xscint_ypcal_r$run_number"+"_s$ftof_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)

			def n_protons_matched=0
			def n_protons_unmatched=0
			def good_delta_t_pcal = null
			if(delta_t_pcal < 3 && delta_t_pcal > -3){
			    good_delta_t_pcal=true
			}

			if(delta_t_pcal <= -3 || delta_t_pcal >= 3){
			    good_delta_t_pcal=false
			}
//println("[ElectronMon::processEvent] - ftof time measured and FTOF and PCAL sector are same")
			//println("[ElectronMon::processEvent] - counting number of protons in event for when timing is matched")			
 			(0..<event.npart).findAll{(event.charge[it] > 0) && (it!=index) && (event.pid[it]==2212) && (event.status[it] >= 2000 && event.status[it] < 4000)}.findResults{hindx->				
			    def hadron_pid = event.pid[hindx]			    
			    def rec_hadron_chi2 = event.chi2pid[hindx]
			    def hadron_p = event.p[hindx]

			    //println("[ElectronMon::processEvent] - found proton")
			   // println("[ElectronMon::processEvent] - " + hadron_pid )
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron_r$run_number",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron2D",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron2D_r$run_number",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    //good ones
			    if( (delta_t_pcal < 3 && delta_t_pcal > -3) && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_matched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron_tmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron2D_tmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)

			    }
			    if( (delta_t_pcal <= -3 || delta_t_pcal >= 3) && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_unmatched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron_tunmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_ftof_pcal_goodelectron2D_tunmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)

			    }
			}

			//println("[ElectronMon::processEvent] - Number of matched protons " + n_protons_matched)
			//println("[ElectronMon::processEvent] - Number of UNmatched protons " + n_protons_unmatched)
			//println("[ElectronMon::processEvent] - status of good_delta_t_pcal " + good_delta_t_pcal)
			if( good_delta_t_pcal && n_protons_matched >= 0 ){
			    //add extra 1 bc started at -1
			    histos.computeIfAbsent("hadron_multiplicities_proton_ftof_pcal_timing_matched",histoBuilders.hadron_multiplicity).fill(n_protons_matched)
			    histos.computeIfAbsent("iron_cross_xftof_ypcal_passall_proton_matched",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)
			}
			
			if( !good_delta_t_pcal && n_protons_unmatched >= 0 ){
			    //println("[ElectronMon::processEvent] - Filling unmatched protons ")
			    //println("[ElectronMon::processEvent] - " + good_delta_t_pcal + " " + n_protons_unmatched + " " + delta_t_scint + " " + delta_t_pcal )
			    histos.computeIfAbsent("hadron_multiplicities_proton_ftof_pcal_timing_unmatched",histoBuilders.hadron_multiplicity).fill(n_protons_unmatched)			    
			    //println("[ElectronMon::processEvent] - Number of entries in hadron multiplicities proton ftof pcal timing unmatched " + histos["hadron_multiplicities_proton_ftof_pcal_timing_unmatched"].getEntries())
			    histos.computeIfAbsent("iron_cross_xftof_ypcal_passall_proton_unmatched",histoBuilders.iron_cross).fill(delta_t_scint, delta_t_pcal)
			    //println("[ElectronMon::processEvent] - Number of entries in hadron iron cross proton ftof pcal timing unmatched " + histos["iron_cross_xftof_ypcal_passall_proton_unmatched"].getEntries())
			    //println("[ElectronMon::processEvent] - Number of OVERFLOW entries in hadron iron cross proton ftof pcal timing unmatched " + histos["iron_cross_xftof_ypcal_passall_proton_unmatched"].projectionX().getOverflow())
			    //println("[ElectronMon::processEvent] - Number of UNDERFLOW entries in hadron iron cross proton ftof pcal timing unmatched " + histos["iron_cross_xftof_ypcal_passall_proton_unmatched"].projectionX().getUnderflow())

			}


		    }
		}
	    }
			
	    //check timing between htcc and pcal
 	    if(!value.contains(false)){
		if( event.cherenkov_status.contains(index) && event.pcal_status.contains(index) ){
		    def pcal_sect = event.pcal_sector[index]
		    def pcal_t = event.pcal_time[index]
		    def pcal_p = event.pcal_path[index]

		    def htcc_t = event.cherenkov_time[index]
		    def htcc_p = event.cherenkov_path[index]/10.0
		    def htcc_sector = event.cherenkov_sector[index]
		    
		    //e_HTCC_vt = e_HTCC_t - e_HTCC_path/29.98f - STT; 
		    // path - time - start 
		    if( pcal_sect == htcc_sector){
			def delta_t_htcc_pcal = pcal_t - htcc_t
			def measured_t_htcc = htcc_t - event.start_time // - event.vt[index]
			def measured_t_pcal = pcal_t - event.vt[index]
			def delta_t_pcal = pcal_p/29.9792458 - measured_t_pcal
			def delta_t_htcc = htcc_p/29.9792458 - measured_t_htcc
			histos.computeIfAbsent("delta_t_htcc_pcal_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(delta_t_htcc_pcal)
			histos.computeIfAbsent("delta_t_htcc_pcal_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(delta_t_htcc_pcal)
			histos.computeIfAbsent("htcc_t_s$htcc_sector"+"_passall",histoBuilders.time).fill(htcc_t)
			histos.computeIfAbsent("htcc_t_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time).fill(htcc_t)
			histos.computeIfAbsent("htcc_t_minus_stt_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(htcc_t-event.start_time)
			histos.computeIfAbsent("htcc_t_minus_stt_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(htcc_t-event.start_time)
			histos.computeIfAbsent("el_p_vs_delta_t_htcc_pcal_s$htcc_sector"+"_passall",histoBuilders.time2d_small).fill(event.p[index], delta_t_htcc_pcal)
			histos.computeIfAbsent("iron_cross_xhtcc_ypcal_passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xhtcc_ypcal_r$run_number"+"_passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xhtcc_ypcal_r$run_number"+"_s$htcc_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)
			histos.computeIfAbsent("iron_cross_xhtcc_ypcal_s$htcc_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)

			def n_protons_matched=0
			def n_protons_unmatched=0
			def good_delta_t_htcc_pcal = null
			if (delta_t_pcal < 3 && delta_t_pcal > -3   && delta_t_htcc < 6 && delta_t_htcc > 1){
			    good_delta_t_htcc_pcal=true
			}
			else if((delta_t_pcal >= 3 || delta_t_pcal <= -3) || (delta_t_htcc >= 6 || delta_t_htcc <= 1) ){
			    good_delta_t_htcc_pcal=false
			}
			//println("[ElectronMon::processEvent] - HTCC and PCAL sector are same")
			//println("[ElectronMon::processEvent] - counting number of protons in event for when timing is matched")			
 			(0..<event.npart).findAll{(event.charge[it] > 0) && (it!=index) && (event.pid[it]==2212) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{hindx->				
			    def hadron_pid = event.pid[hindx]			    
			    def rec_hadron_chi2 = event.chi2pid[hindx]
			    def hadron_p = event.p[hindx]
			  //  println("[ElectronMon::processEvent] - found proton")
			  //  println("[ElectronMon::processEvent] - " + hadron_pid )
 			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron_r$run_number",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron2D",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron2D_r$run_number",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    //good ones
			    if( delta_t_pcal < 3 && delta_t_pcal > -3 && delta_t_htcc < 6 && delta_t_htcc > 1 && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_matched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron_tmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron2D_tmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    }
			    else if( ((delta_t_pcal >= 3 || delta_t_pcal <= -3) || (delta_t_htcc >= 6 || delta_t_htcc <= 1)) &&  Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_unmatched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron_tunmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_pcal_goodelectron2D_tunmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    }
			}

			//println("[ElectronMon::processEvent] - Number of matched HTCC/PCAL protons " + n_protons_matched)
			//println("[ElectronMon::processEvent] - Number of UNmatched HTCC/PCAL protons " + n_protons_unmatched)

			if( good_delta_t_htcc_pcal && n_protons_matched >= 0 ){
			    //add extra 1 bc started at -1
			    histos.computeIfAbsent("hadron_multiplicities_proton_htcc_pcal_timing_matched",histoBuilders.hadron_multiplicity).fill(n_protons_matched)
			    histos.computeIfAbsent("iron_cross_xhtcc_ypcal_passall_proton_matched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)
			}
			if( !good_delta_t_htcc_pcal && n_protons_unmatched >= 0 ){
			    histos.computeIfAbsent("hadron_multiplicities_proton_htcc_pcal_timing_unmatched",histoBuilders.hadron_multiplicity).fill(n_protons_unmatched)
			    histos.computeIfAbsent("iron_cross_xhtcc_ypcal_passall_proton_unmatched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_pcal)
			}
		    }
		}
	    }

	    //check timing between htcc and ftof
 	    if(!value.contains(false)){
		if( event.cherenkov_status.contains(index) && event.tof_status.contains(index) ){
		    def htcc_t = event.cherenkov_time[index]
		    def htcc_p = event.cherenkov_path[index]/10.0
		    def htcc_sector = event.cherenkov_sector[index]
		    
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[index].layer
		    def my_layer = []
		    
		    def ftof_t = -1000
		    def ftof_p = -1000
		    def ftof_sector = -1
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }
		 
		    if( my_layer.contains('2') && !my_layer.contains('3') && !my_layer.contains('1') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    
		    event.tof[index].eachWithIndex{ hh, ii ->  
			if(layer_to_use!=-1 && 2 == hh.layer && hh.layer!=3 && hh.layer!=0 && layer_to_use == 2){// hh.layer == 2 ){
 			    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			    ftof_sector = hh.sector
			    ftof_t = hh.time
			    ftof_p = hh.path
			}
		    }
		    if( ftof_t != -1000 && ftof_sector == htcc_sector){

			def delta_t_htcc_ftof = ftof_t - htcc_t
			def measured_t_htcc = htcc_t - event.start_time // - event.vt[index]
			def measured_t_ftof = ftof_t - event.vt[index]
		 	def delta_t_ftof = ftof_p/29.9792458 - measured_t_ftof
			def delta_t_htcc = htcc_p/29.9792458 - measured_t_htcc
			histos.computeIfAbsent("delta_t_htcc_ftof_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(delta_t_htcc_ftof)
			histos.computeIfAbsent("delta_t_htcc_ftof_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(delta_t_htcc_ftof)
			histos.computeIfAbsent("htcc_t_s$htcc_sector"+"_passall",histoBuilders.time).fill(htcc_t)
			histos.computeIfAbsent("htcc_t_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time).fill(htcc_t)
			histos.computeIfAbsent("htcc_t_minus_stt_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(htcc_t-event.start_time)
			histos.computeIfAbsent("htcc_t_minus_stt_r$run_number"+"_s$htcc_sector"+"_passall",histoBuilders.time_small).fill(htcc_t-event.start_time)
			histos.computeIfAbsent("el_p_vs_delta_t_htcc_ftof_s$htcc_sector"+"_passall",histoBuilders.time2d_small).fill(event.p[index], delta_t_htcc_ftof)
			histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
			histos.computeIfAbsent("iron_cross_xhtcc_yftof_r$run_number"+"_passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
			histos.computeIfAbsent("iron_cross_xhtcc_yftof_r$run_number"+"_s$htcc_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
			histos.computeIfAbsent("iron_cross_xhtcc_yftof_s$htcc_sector"+"passall",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)


			def n_protons_matched=0
			def n_protons_unmatched=0
			def good_delta_t_htcc = null
			if (delta_t_htcc < 6 && delta_t_htcc > 1){
			    good_delta_t_htcc=true
			}
			else if( (delta_t_htcc >= 6 || delta_t_htcc <= 1) ){
			    good_delta_t_htcc=false
			}

			//println("[ElectronMon::processEvent] - HTCC and FTOF sector are same")
			//println("[ElectronMon::processEvent] - counting number of protons in event for when timing is matched")			
 			(0..<event.npart).findAll{(event.charge[it] > 0) && (it!=index) && (event.pid[it]==2212) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{hindx->
			    def hadron_pid = event.pid[hindx]			    
			    def rec_hadron_chi2 = event.chi2pid[hindx]
			    def hadron_p = event.p[hindx]
			    def hadron_status = event.status[hindx]
			   
			    //println("[ElectronMon::processEvent] - found proton")
			    //println("[ElectronMon::processEvent] - " + hadron_pid )
 			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_r$run_number",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D_r$run_number",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    //good ones
			    if( delta_t_htcc >= 1 && delta_t_htcc <= 6 && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_matched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_tmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D_tmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)

			    }
			    else if( (delta_t_htcc < 1 || delta_t_htcc > 6) && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_unmatched+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_tunmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D_tunmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)

			    }
			}

			//println("[ElectronMon::processEvent] - Number of matched HTCC/FTOF protons " + n_protons_matched)
			//println("[ElectronMon::processEvent] - Number of UNmatched HTCC/FTOF protons " + n_protons_unmatched)
			
			if(good_delta_t_htcc && n_protons_matched >= 0 ){
			    //add extra 1 bc started at -1
			    histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_matched",histoBuilders.hadron_multiplicity).fill(n_protons_matched)
			    histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_matched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
			}
			if( !good_delta_t_htcc &&  n_protons_unmatched >= 0 ){
			    histos.computeIfAbsent("hadron_multiplicities_proton_htcc_ftof_timing_unmatched",histoBuilders.hadron_multiplicity).fill(n_protons_unmatched)
			    histos.computeIfAbsent("iron_cross_xhtcc_yftof_passall_proton_unmatched",histoBuilders.iron_cross).fill(delta_t_htcc, delta_t_ftof)
			}
		    }
		}
	    }

	    
	    /////////////////////////////////////////
	    // compare htcc, ftof, and pcal timing
	    if(!value.contains(false)){
		if( event.cherenkov_status.contains(index) && event.tof_status.contains(index) && event.pcal_status.contains(index) ){
		    ///////////////////////////
		    //get htcc info
		    def htcc_t = event.cherenkov_time[index]
		    def htcc_p = event.cherenkov_path[index]/10.0
		    def htcc_sector = event.cherenkov_sector[index]
		    

		    //////////////////////////
		    //get ftof info
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[index].layer
		    def my_layer = []
		    
		    def ftof_t = -1000
		    def ftof_p = -1000
		    def ftof_sector = -1
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }


		    if( my_layer.contains('2') && !my_layer.contains('3') && !my_layer.contains('1') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    

		    event.tof[index].eachWithIndex{ hh, ii ->  
			if(layer_to_use!=-1 && 2 == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 			    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			    ftof_sector = hh.sector
			    ftof_t = hh.time
			    ftof_p = hh.path
			}
		    }

		    

		    ///////////////////////////
		    //get pcal info
		    def pcal_sect = event.pcal_sector[index]
		    def pcal_t = event.pcal_time[index]
		    def pcal_p = event.pcal_path[index]
	
		    if( ftof_sector != -1 && ftof_t != -1000 && ftof_sector == htcc_sector && ftof_sector == pcal_sect){
			def measured_t_htcc = htcc_t - event.start_time // - event.vt[index]
			def delta_t_htcc = htcc_p/29.9792458 - measured_t_htcc
		
			def measured_t_ftof = ftof_t - event.vt[index]
		 	def delta_t_ftof = ftof_p/29.9792458 - measured_t_ftof
		
			def measured_t_pcal = pcal_t - event.vt[index]
			def delta_t_pcal = pcal_p/29.9792458 - measured_t_pcal
		
			def diff_t_htcc_pcal = delta_t_htcc - delta_t_pcal
			def diff_t_ftof_pcal = delta_t_ftof - delta_t_pcal
			def diff_t_htcc_ftof = delta_t_htcc - delta_t_ftof
			def diff_t_ftof_htcc = delta_t_ftof - delta_t_htcc

			def el_trigger_htcc_vertex_time = getStartTime(getVertexTime(htcc_t, htcc_p, 1), event.vz[index])
			def el_trigger_ftof_vertex_time = getStartTime(getVertexTime(ftof_t, ftof_p, 1), event.vz[index])
			def el_trigger_pcal_vertex_time = getStartTime(getVertexTime(pcal_t, pcal_p, 1), event.vz[index])

			//println("PCAL TIME " + pcal_t + " PATH LENGTH " + pcal_p + " VERTEX POSITION " + event.vz[index])
			//println("[ElectronMon::processEvent] - HTCC vertex time " + el_trigger_htcc_vertex_time )
			///println("[ElectronMon::processEvent] - FTOFL2 vertex time " + el_trigger_ftof_vertex_time )
			//println("[ElectronMon::processEvent] - PCAL vertex time " + el_trigger_pcal_vertex_time )
			//println("[ElectronMon::processEvent] - event.vt[index] " + event.vt[index])

 			histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof_pcal)
			histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_htcc_ftof)
			histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,diff_t_ftof_htcc)
			histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_delta_t_ftof_passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,delta_t_ftof)
			histos.computeIfAbsent("iron_cross_htcc_minus_ftof_vs_delta_t_pcal_passall",histoBuilders.iron_cross).fill(diff_t_htcc_ftof,delta_t_pcal)
			histos.computeIfAbsent("iron_cross_ftof_minus_pcal_vs_delta_t_htcc_passall",histoBuilders.iron_cross).fill(diff_t_ftof_pcal,delta_t_htcc)
						
			histos.computeIfAbsent("iron_cross_htcc_minus_pcal_vs_ftof_minus_pcal_s$pcal_sect"+"passall",histoBuilders.iron_cross).fill(diff_t_htcc_pcal,diff_t_ftof_pcal)
					
			//y - cut @ 3.9 +- 3*0.7
			//x - cut @ -0.2459 +- 3*0.65
			def min_ftof_pcal_cut = -0.2459 - 3*0.65
			def max_ftof_pcal_cut = -0.2459 + 3*0.65

			def min_htcc_ftof_cut = 3.9 - 3*0.7
			def max_htcc_ftof_cut = 3.9 + 3*0.7
			

			def n_protons_matched_box_cut=0
			def n_protons_unmatched_box_cut=0
			def good_timing_box_cut = null
			def n_protons_matched_corner_cut=0
			def n_protons_unmatched_corner_cut=0
			def good_timing_corner_cut = null
			
			if ( diff_t_ftof_pcal < max_ftof_pcal_cut && diff_t_ftof_pcal > min_ftof_pcal_cut && diff_t_htcc_ftof < max_htcc_ftof_cut && diff_t_htcc_ftof > min_htcc_ftof_cut ) {
			    good_timing_box_cut=true
			}
			else if( diff_t_ftof_pcal < max_ftof_pcal_cut || diff_t_ftof_pcal > min_ftof_pcal_cut || diff_t_htcc_ftof < max_htcc_ftof_cut || diff_t_htcc_ftof > min_htcc_ftof_cut ) {
			    good_timing_box_cut=false
			}

			//3.3, -0.37
			if( (diff_t_ftof_pcal > 3.3 && diff_t_htcc_ftof < -0.37) || ( diff_t_ftof_pcal < -3.3 && diff_t_htcc_ftof > 7.1 ) ){
			    good_timing_corner_cut=true			    
			}
			else{
			    good_timing_corner_cut=false
			}
			    			
			//println("[ElectronMon::processEvent] - HTCC and FTOF sector are same")
			//println("[ElectronMon::processEvent] - counting number of protons in event for when timing is matched")			

			// REDEFINE THE BAD ELECTRON EVENTS HERE //
			if( good_timing_corner_cut ){
			    (0..<event.npart).findAll{(event.charge[it] > 0) && (it!=index) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{hindx->
				def hadron_pid = event.pid[hindx]			    
				def rec_hadron_chi2 = event.chi2pid[hindx]
				def hadron_p = event.p[hindx]
				def hadron_vz = event.vz[hindx]
				def hadron_vt = event.vt[hindx]
				def hadron_beta = event.beta[hindx]
				def hadron_status = event.status[hindx]
				def theory_beta = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.938*0.938)
				
 				def hadron_vertex_corr_htcc_start_time = getStartTime(el_trigger_htcc_vertex_time, hadron_vz)
 				def hadron_vertex_corr_ftof_start_time = getStartTime(el_trigger_ftof_vertex_time, hadron_vz)
 				def hadron_vertex_corr_pcal_start_time = getStartTime(el_trigger_pcal_vertex_time, hadron_vz)
								
				if ( Math.abs(rec_hadron_chi2) < 4.2 ){
				    def pr_ftof_t = -1000
				    def pr_ftof_p = -1000
				    def pr_ftof_sector = -1
				    
				    if( event.tof_status.contains(hindx) ){
					//////////////////////////
					//get ftof info
					def pr_layer_to_use = -1
					//println(' layers available ' + event.tof[index].layer)
					def pr_layer_avail = event.tof[hindx].layer
					def pr_my_layer = []
					
					for( ii in pr_layer_avail){
					    pr_my_layer.add(Integer.toString(ii))
					}
					if( pr_my_layer.contains('2') && !pr_my_layer.contains('3') ){
					    //println('use layer 2 for timing ' )
					    pr_layer_to_use = 2
					}
					else if ( !pr_my_layer.contains('2') && pr_my_layer.contains('1') && !pr_my_layer.contains('3')  ){
					    //println('use layer 1 for timing ' )
					    pr_layer_to_use = 1
					}		    
					event.tof[hindx].eachWithIndex{ hh, ii ->  
					    if(pr_layer_to_use!=-1 && 2 == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 						//println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
						pr_ftof_sector = hh.sector
						pr_ftof_t = hh.time
						pr_ftof_p = hh.path	    
					    }
					}
					
				    }
				    if( pr_ftof_t > 0 ){
					//println("[ElectronMon::processEvent] - hadron vertex time HTCC " + hadron_vertex_corr_htcc_start_time )
					//println("[ElectronMon::processEvent] - hadron vertex time FTOFL2 " + hadron_vertex_corr_ftof_start_time )
					//println("[ElectronMon::processEvent] - hadron vertex time PCAL " + hadron_vertex_corr_pcal_start_time )
				
 					def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
					def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
					def hadron_ftof_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof_start_time)/speed_of_light
					def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
					def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light

					histos.computeIfAbsent("pos_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_beta)
 					histos.computeIfAbsent("pos_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
					histos.computeIfAbsent("pos_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 					histos.computeIfAbsent("pos_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_ftof_beta)
					histos.computeIfAbsent("pos_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
					histos.computeIfAbsent("pos_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,calc_rec_beta)				
				    }
				    
				}    
			    }
			}
			/////////// DONE REDEFINING EVENTS ///////////////
			
 			(0..<event.npart).findAll{(event.charge[it] > 0) && (it!=index) && (event.pid[it]==2212) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{hindx->
			    def hadron_pid = event.pid[hindx]			    
			    def rec_hadron_chi2 = event.chi2pid[hindx]
			    def hadron_p = event.p[hindx] 
			    def hadron_vz = event.vz[hindx]
			    def hadron_vt = event.vt[hindx]
			    def hadron_beta = event.beta[hindx]
			    def hadron_status = event.status[hindx]
			    def theory_beta = hadron_p/Math.sqrt(hadron_p*hadron_p + 0.938*0.938)

 			    def hadron_vertex_corr_htcc_start_time = getStartTime(el_trigger_htcc_vertex_time, hadron_vz)
 			    def hadron_vertex_corr_ftof_start_time = getStartTime(el_trigger_ftof_vertex_time, hadron_vz)
 			    def hadron_vertex_corr_pcal_start_time = getStartTime(el_trigger_pcal_vertex_time, hadron_vz)
			    
			    
			    
			    //////////////////////////
			    //get ftof info
			    def pr_layer_to_use = -1
			    //println(' layers available ' + event.tof[index].layer)
			    def pr_layer_avail = event.tof[hindx].layer
			    def pr_my_layer = []
			    
			    def pr_ftof_t = -1000
			    def pr_ftof_p = -1000
			    def pr_ftof_sector = -1
			    for( ii in pr_layer_avail){
				pr_my_layer.add(Integer.toString(ii))
			    }
			    if( pr_my_layer.contains('2') && !pr_my_layer.contains('3') ){
				//println('use layer 2 for timing ' )
				pr_layer_to_use = 2
			    }
			    else if ( !pr_my_layer.contains('2') && pr_my_layer.contains('1') && !pr_my_layer.contains('3')  ){
				//println('use layer 1 for timing ' )
				pr_layer_to_use = 1
			    }		    
			    event.tof[hindx].eachWithIndex{ hh, ii ->  
				if(pr_layer_to_use!=-1 && 2 == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
 				    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
				    pr_ftof_sector = hh.sector
				    pr_ftof_t = hh.time
				    pr_ftof_p = hh.path
				}
			    }
			    			  
			    //println("[ElectronMon::processEvent] - found proton")
			    //println("[ElectronMon::processEvent] - " + hadron_pid )
 			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_r$run_number",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D_r$run_number",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    //good ones
			    if( good_timing_box_cut &&  Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_matched_box_cut+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron2D_tmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)
			    
				//println("[ElectronMon::processEvent] - HTCC vertex time " + el_trigger_htcc_vertex_time )
				//println("[ElectronMon::processEvent] - FTOFL2 vertex time " + el_trigger_ftof_vertex_time )
				//println("[ElectronMon::processEvent] - PCAL vertex time " + el_trigger_pcal_vertex_time )
				//println("[ElectronMon::processEvent] - event.vt[index] " + event.vt[index])
				
				//println("[ElectronMon::processEvent] - hadron vertex time HTCC " + hadron_vertex_corr_htcc_start_time )
				//println("[ElectronMon::processEvent] - hadron vertex time FTOFL2 " + hadron_vertex_corr_ftof_start_time )
				//println("[ElectronMon::processEvent] - hadron vertex time PCAL " + hadron_vertex_corr_pcal_start_time )
				//println("[ElectronMon::processEvent] - hadron event.vt[index] " + event.vt[hindx])

 				if( pr_ftof_t > 0 ){
			
				    def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
				    def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
				    def hadron_ftof_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof_start_time)/speed_of_light
				    def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
				    def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light
				    
				    //println("[ElectronMon::processEvent] - hadron new beta with pcal start time " + hadron_new_beta)
				    //println("[ElectronMon::processEvent] - hadron new beta with start time " + hadron_new_beta_start_time)
				    
				    histos.computeIfAbsent("proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,hadron_beta)
 				    histos.computeIfAbsent("proton_beta_p_redefined_pcal_st_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_htcc_st_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 				    histos.computeIfAbsent("proton_beta_p_redefined_ftof_st_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,hadron_ftof_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
				    histos.computeIfAbsent("proton_beta_p_redefined_vt_st_htcc_ftof_pcal_goodelectron_tmatched",histoBuilders.betap).fill(hadron_p,calc_rec_beta)
				    
				    
				    
				    //println("[ElectronMon::processEvent] - calc rec beta - " + calc_rec_beta)
				    //println("[ElectronMon::processEvent] - rec beta " + hadron_beta)
				    //println("[ElectronMon::processEvent] - theory beta " + theory_beta)
				}
				
			    }
			    else if( !good_timing_box_cut && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_unmatched_box_cut+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron2D_tunmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)		
				
 				if( pr_ftof_t > 0 ){								    
				    def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
				    def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
				    def hadron_ftof_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof_start_time)/speed_of_light
				    def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
				    def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light
				    
				    //println("[ElectronMon::processEvent] - hadron new beta with pcal start time " + hadron_new_beta)
				    //println("[ElectronMon::processEvent] - hadron new beta with start time " + hadron_new_beta_start_time)
				    
				    histos.computeIfAbsent("proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_beta)
				    histos.computeIfAbsent("proton_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
				    histos.computeIfAbsent("proton_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 				    histos.computeIfAbsent("proton_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_ftof_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
				    histos.computeIfAbsent("proton_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched",histoBuilders.betap).fill(hadron_p,calc_rec_beta)
				}

			    }
			
			    
			    if( good_timing_corner_cut && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_matched_corner_cut+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron_corner_cut_tmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron2D_corner_cut_tmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)

 				if( pr_ftof_t > 0 ){
				    
				    def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
				    def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
				    def hadron_ftof_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof_start_time)/speed_of_light
				    def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
				    def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light
				    
				    //println("[ElectronMon::processEvent] - hadron new beta with pcal start time " + hadron_new_beta)
				    //println("[ElectronMon::processEvent] - hadron new beta with start time " + hadron_new_beta_start_time)
				    
				    histos.computeIfAbsent("proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_pcal_st_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_htcc_st_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
 				    histos.computeIfAbsent("proton_beta_p_redefined_ftof_st_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_ftof_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
				    histos.computeIfAbsent("proton_beta_p_redefined_vt_st_htcc_ftof_pcal_goodelectron_tmatched_corner",histoBuilders.betap).fill(hadron_p,calc_rec_beta)
				}
			    }
			    else if( !good_timing_corner_cut && Math.abs(rec_hadron_chi2) < 4.2 ){
				n_protons_unmatched_corner_cut+=1
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron_corner_cut_tunmatched",histoBuilders.recpart_chi2pid).fill(rec_hadron_chi2)
				histos.computeIfAbsent("recpart_chi2pid_proton_pid_htcc_ftof_pcal_goodelectron2D_corner_cut_tunmatched",histoBuilders.recpart_chi2pid2D).fill(hadron_p,rec_hadron_chi2)		

 				if( pr_ftof_t > 0 ){
				    
				    def hadron_pcal_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_pcal_start_time)/speed_of_light
				    def hadron_htcc_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_htcc_start_time)/speed_of_light
				    def hadron_ftof_beta = pr_ftof_p/(pr_ftof_t - hadron_vertex_corr_ftof_start_time)/speed_of_light
				    def hadron_new_beta_start_time = pr_ftof_p/(pr_ftof_t - event.start_time)/speed_of_light
				    def calc_rec_beta = pr_ftof_p/(pr_ftof_t - event.vt[hindx])/speed_of_light
				    
				    //println("[ElectronMon::processEvent] - hadron new beta with pcal start time " + hadron_new_beta)
				    //println("[ElectronMon::processEvent] - hadron new beta with start time " + hadron_new_beta_start_time)
				    
				    histos.computeIfAbsent("proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_beta)
				    histos.computeIfAbsent("proton_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_pcal_beta)
				    histos.computeIfAbsent("proton_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_htcc_beta)
				    histos.computeIfAbsent("proton_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_ftof_beta)
				    histos.computeIfAbsent("proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,hadron_new_beta_start_time)
				    histos.computeIfAbsent("proton_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched_corner",histoBuilders.betap).fill(hadron_p,calc_rec_beta)
				}			
			    }


			}

			//println("[ElectronMon::processEvent] - Number of matched HTCC/FTOF protons " + n_protons_matched)
			//println("[ElectronMon::processEvent] - Number of UNmatched HTCC/FTOF protons " + n_protons_unmatched)
			
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
		
    
		    }

		}
	    }
			/*
	    if( event.tof_status.contains(index) ){
		//def ftof_sector = event.tof_sector[index]
		//def ftof_t = event.tof_time[index]
		//def ftof_p = event.tof_path[index]
		//def ftof_paddle = event.tof_paddle[index]
		//def ftof_layer = event.tof_layer[index]
		//def ftof_t_corr = ftof_t - ftof_p/(29.9792458)
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){		   
			//histos.computeIfAbsent("ftoft_cut$cut_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
			//histos.computeIfAbsent("ftoft_cut$cut_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)
			//histos.computeIfAbsent("ftof_t_corr_cut$cut_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
 			//histos.computeIfAbsent("ftof_t_corr_cut$cut_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t_corr)
			//histos.computeIfAbsent("ftoft_cut$cut_index"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_layer,ftof_t)
			//histos.computeIfAbsent("ftof_t_corr_cut$cut_index"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_layer,ftof_t_corr)
		    }
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("ftof_betatime_passall_s$ftof_sector",histoBuilders.betatime).fill(event.beta[index],ftof_t)
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector",histoBuilders.time).fill(ftof_t)
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)		    
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector"+"_paddle$ftof_paddle",histoBuilders.time).fill(ftof_t)		    
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector"+"_per_paddle",histoBuilders.time2d_large).fill(ftof_paddle,ftof_t)
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector"+"_per_layer",histoBuilders.time2d_small).fill(ftof_layer,ftof_t)
	         }
		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
 		    if( cut_comb_value.equals(value) ){
			//histos.computeIfAbsent("ftoft_cut_comb$cut_comb_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
			//histos.computeIfAbsent("ftoft_cut_comb$cut_comb_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)
			//histos.computeIfAbsent("ftof_t_corr_cut_comb$cut_comb_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
 			//histos.computeIfAbsent("ftof_t_corr_cut_comb$cut_comb_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t_corr)
			//histos.computeIfAbsent("ftoft_cut_comb_$cut_comb_index"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_layer,ftof_t)
			//histos.computeIfAbsent("ftof_t_corr_cut_comb_$cut_comb_index"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_layer,ftof_t_corr)
		    }
		}
	    }
	    */
	    	    	    
	    //generic p, theta, phi, vx, vy, vz information
	    //also CHI2 and NDF
	    if( event.npart>0 && event.dc1_status.contains(index)) {
		def dc_sector  = event.dc_sector[index]		    
		def dc_chi2 = event.dc_chi2[index]
		def dc_ndf = event.dc_ndf[index]
		def vx = event.vx[index]
		def vy = event.vy[index]
		def measured_el = LorentzVector.withPID(11,event.px[index],event.py[index],event.pz[index])
		def delta_energy = KinTool.delta_meas_energy(beam, measured_el) 
		def delta_theta = KinTool.delta_meas_theta(beam, measured_el)
		
 		value.eachWithIndex{ cut_result, cut_index ->
		    if( cut_result && cut_index==0 ){
			histos.computeIfAbsent("vz_cut$cut_index"+"_s$dc_sector",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent("ptheta_cut$cut_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent("vztheta_cut$cut_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent("vxy_cut$cut_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent("vzphi_cut$cut_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			histos.computeIfAbsent("chi_cut$cut_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 			histos.computeIfAbsent("ndf_cut$cut_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    
		    }
		}
		if(!value.contains(false)){		    
		    histos.computeIfAbsent("vz_passall"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_passall"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_passall",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_passall"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall",histoBuilders.vzphi).fill(phi,vz)
 		    /*histos.computeIfAbsent("chi_passall"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 		    histos.computeIfAbsent("ndf_passall"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		
 		    histos.computeIfAbsent("delta_energy_passall_$dc_sector",histoBuilders.deltakin).fill(delta_energy)
		    histos.computeIfAbsent("delta_theta_passall_$dc_sector",histoBuilders.deltakin).fill(delta_theta)
 		    histos.computeIfAbsent("delta_energy_vs_chi2_passall_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_energy)
		    histos.computeIfAbsent("delta_theta_vs_chi2_passall_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_theta)
		    histos.computeIfAbsent("delta_energy_vs_p_passall_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_energy)
		    histos.computeIfAbsent("delta_theta_vs_p_passall_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_theta)
		*/
		}

		if( everythingElsePassed(value,vz_cut_index) ){
		    histos.computeIfAbsent("vz_passall_passall_but_vz"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_passall_but_vz"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vz_passall_passall_but_vz",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_passall_but_vz",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall_but_vz",histoBuilders.vzphi).fill(phi,vz)
		}
		
		if( everythingElsePassed(value,min_mntm_cut_index) ){
		    histos.computeIfAbsent("ptheta_passall_but_minmntm"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		}
		
		if( value[10] && value[11] ){		    		   
		    histos.computeIfAbsent("vz_pass_eventbuilder"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_pass_eventbuilder"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vz_pass_eventbuilder",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_pass_eventbuilder",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_pass_eventbuilder",histoBuilders.vzphi).fill(phi,vz)
		    histos.computeIfAbsent("ptheta_pass_eventbuilder"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)				
		}



		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
 		    if( cut_comb_value.equals(value) ){
/*			histos.computeIfAbsent("vz_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent("ptheta_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			histos.computeIfAbsent("phitheta_cut_comb$cut_comb_index",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent("vztheta_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent("vxy_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent("vzphi_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			histos.computeIfAbsent("chi_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 			histos.computeIfAbsent("ndf_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    
			histos.computeIfAbsent("delta_energy_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin).fill(delta_energy)
			histos.computeIfAbsent("delta_theta_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin).fill(delta_theta)
 			histos.computeIfAbsent("delta_energy_vs_chi2_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_energy)
			histos.computeIfAbsent("delta_theta_vs_chi2_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_theta)
			histos.computeIfAbsent("delta_energy_vs_p_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_energy)
			histos.computeIfAbsent("delta_theta_vs_p_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_theta)
*/		    }
		}
	    }	
	}

	if(event.npart>0){		
	    histos.computeIfAbsent("event_helicity",histoBuilders.helicity).fill(event.helicity)
	    histos.computeIfAbsent("event_start_time",histoBuilders.time).fill(event.start_time)
	    histos.computeIfAbsent("event_rf_time",histoBuilders.time).fill(event.rf_time)		    		
	    event.charge.findResults{ it.value < 0 ? it.value : null}?.each{histos.computeIfAbsent("event_charge_neg",histoBuilders.helicity).fill(it)}
	    event.charge.findResults{ it.value > 0 ? it.value : null}?.each{histos.computeIfAbsent("event_charge_pos",histoBuilders.helicity).fill(it)}
	}		
    }

}
