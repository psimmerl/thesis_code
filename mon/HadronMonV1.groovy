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

class HadronMon{

    def beam = 10.604

    //call electron cut constructor
    def electron = new ElectronFromEvent();

    def myElectronCutStrategies = [
	electron.passElectronStatus,
	electron.passElectronChargeCut,
	electron.passElectronEBPIDCut,
	//electron.passElectronNpheCut,
	//electron.passElectronVertexCut
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronEIEOCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3
    ]


    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.3) },
	uvwsf : { title -> new H2F("$title","$title", 500, 0.0, 500.0, 500, 0.0, 0.3) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 50, 0.0, 50 ) },
	nphe : { title -> new H1F("$title","$title", 40, 0.0, 40.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 1.0, 1000, 0.0, 1.0 ) },
	vz   : { title -> new H1F("$title","$title", 400, -40, 40) },
	vxy : { title -> new H2F("$title","$title", 400, -1, 1,400, -1, 1) },
	
	time : { title -> new H1F("$title","$title", 2000, -450, 450) },
	time2d : { title -> new H2F("$title","$title", 61, 0.0, 61.0, 2000, -450, 450) },
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
	deltat : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) },
	deltab : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -1, 1  ) },
	betap : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 1.4  ) },
	deltab_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -1, 1  ) },
	deltat_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -2, 2  ) },
	calib_time : { title -> new H2F("$title","$title", 1000, -100.0, 100, 1000, -100, 100  ) },
	betacomp : { title -> new H2F("$title","$title", 100, 0.80, 1.2, 100, 0.80, 1.2  ) },
	vz_vt : { title -> new H2F("$title","$title", 100, -70, 70, 100, -1.1, 101  ) },
	rf_vt : { title -> new H2F("$title","$title", 100, -70, 70, 100, -1.1, 101  ) },

	sum_time : { title -> new H1F("$title","$title", 5000, 0.0, 5000  ) }

    ]

    def processEvent( event, hevent ){

	//println( ' processing event ' )
	def start_time = event.start_time
	def rf_time = event.rf_time
	if( start_time != null && rf_time != null ){
	    histos.computeIfAbsent("start_time",histoBuilders.time).fill(start_time)
	    histos.computeIfAbsent("rf_time",histoBuilders.time).fill(rf_time)
	}
		
	def electron_candidates = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  
	
	def ev_sum_tof_time = 0
 	def time_hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    //look at positive charged hadrons
	    (0..<event.npart).findAll{event.charge[it] > 0 && ( it != indx )}.findResults{index->		
		def vt = event.vt[index]
		ev_sum_tof_time+=vt
		
	    }
	}

	if ( ev_sum_tof_time < 86.616 + 5){


 	def hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    //look at positive charged hadrons
	    (0..<event.npart).findAll{event.charge[it] > 0 && ( it != indx )}.findResults{index->		
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def p = lv.mag()
		def beta = event.beta[index]
		def vz = event.vz[index]
		def start_time_rf_corrected = event.vt[index] // start time with vertex correction
		def elvz = event.vz[indx]
		def vt = event.vt[index]


		def run_number = event.run_number
		
		if( event.dc1_status.contains(index) ){		    
		    def hits = event.dc1[index]
		    def hit = hits.find{it.layer==12}.each{
 			histos.computeIfAbsent("dcr1hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc2_status.contains(index) ){		    
		    def hits = event.dc2[index]
		    def hit = hits.find{it.layer==24}.each{
 			histos.computeIfAbsent("dcr2hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc3_status.contains(index) ){		    
		    def hits = event.dc3[index]
		    def hit = hits.find{it.layer==36}.each{
 			histos.computeIfAbsent("dcr3hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}
		
		//get electron tof info
		if(event.tof_status.contains(index) ){
		   
		    //println('sector ' + event.tof[index].sector + ' paddle ' + event.tof[index].paddle + ' layer ' + event.tof[index].layer )


		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[index].layer
		    def my_layer = []
		    
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }
			

		    //println(layer_avail)
		    //println(my_layer )
		    //println(' -------> checking if layer_list contains number 2: ' + my_layer.contains('2') )
		    if( my_layer.contains('2') && !my_layer.contains('3') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }

		    

		    event.tof[index].eachWithIndex{ hh, ii ->  
			//println(' index ' + ii )
			//println(' --> USING LAYER ' + layer_to_use + ' FOR TIMING ' )
			if ( layer_to_use == -1 ){ 			    
			    //println(' --> NO GOOD LAYER TO USE -- SKIP THIS TRACK -- ' ) 
			}

			if(layer_to_use!=-1 && layer_to_use == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
			    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			    def measured_beta = beta// ftof_p/(ftof_t - start_time)/(29.9792458)
 			    if( measured_beta > 1.1 && measured_beta < 1.23 && p < 2.0){
			    }

			    if( Math.abs(vz-elvz) < 5.0 ) { //p < 1.0 && measured_beta > 0.91 ){
				
  				//histos.computeIfAbsent("ftof_betap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 				//histos.computeIfAbsent("ftof_betap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_layer$layer_to_use"+"_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 				//histos.computeIfAbsent("ftof_betap_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
  				histos.computeIfAbsent("vz_vt_$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  				histos.computeIfAbsent("rf_vt_$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)

				if( vt > 80.0 && vt < 94 ){
 				    histos.computeIfAbsent("ftof_RECbetap_pos_vtcut_layer$layer_to_use"+"_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				    histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut_vtcut",histoBuilders.betap).fill(p, measured_beta)
 				    histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut_vtcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
  				    histos.computeIfAbsent("vz_vt_vtcut_$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  				    histos.computeIfAbsent("rf_vt_vtcut_$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)

				}

				if( vt > 80.0 && vt < 94 && vz < 10 && vz > -12 ){
 				    histos.computeIfAbsent("ftof_RECbetap_pos_vtcut_vzcut_layer$layer_to_use"+"_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				    histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut_vtcut_vzcut",histoBuilders.betap).fill(p, measured_beta)
 				    histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut_vtcut_vzcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
  				    histos.computeIfAbsent("vz_vt_vtcut_vzcut_$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  				    histos.computeIfAbsent("rf_vt_vtcut_vzcut_$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)

				}

			    }
			    else{
				histos.computeIfAbsent("ftof_RECbetap_pos_layer$layer_to_use"+"_BADelvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_BADelvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_BADelvzcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
			    }
				
 			    histos.computeIfAbsent("ftof_betap_pos_layer$layer_to_use",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_betap_pos",histoBuilders.betap).fill(p, measured_beta)
 			    
			    

			    def shiftp= 0
			    def shift=0
			    if(hh.layer == 2 ){  
				shiftp=0
				shift=0 // -4
			    }
			    def ftof_sector = hh.sector
			    def ftof_t = hh.time
			    def ftof_p = hh.path+shiftp
			    def ftof_paddle = hh.paddle
			    def ftof_layer = hh.layer
			    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    		   
			    		    
			    def measured_t = ftof_t - start_time 
			    def measured_beta2 = ftof_p/(ftof_t - start_time_rf_corrected + shift)/(29.9792458)
			    def measured_beta2shift = ftof_p/(ftof_t - start_time - 2)/(29.9792458)
			    if( true ){// p < 1.0 && measured_beta2 > 0.91 ){
				
 				if( measured_beta > 1.1 && measured_beta < 1.23 && p < 2.0 && Math.abs(measured_beta - measured_beta2) > 0.1 && Math.abs(vz-elvz) < 5.0 ){
				    //print(' ---------> PARTICLE INDEX ' + index )
				    //hevent.getBank("RUN::config").show()
				    //hevent.getBank("REC::Event").show()
				    //hevent.getBank("REC::Particle").show()
				    //hevent.getBank("REC::Scintillator").show()

				    //println(' using layer ' + layer_to_use )
				    //println(' measured beta 2 '  + measured_beta2 )
				    //println(' measured beta2 shifted by -2ns ' + measured_beta2shift )
				    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy) 
				}

				
   				histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)

 				if( Math.abs(vz-elvz) < 5.0 ){
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_elvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_elvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_mybetap_pos_elvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)

				}
				else{
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
  				    histos.computeIfAbsent("ftof_mybetap_pos_BADelvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)

				}
				histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)

				
			    }

/*
			    ///////////////////////////////////
			    //  look at events with beta > 1.1 here for this run to apply a time simple shift
			    if( measured_beta > 1.1 && measured_beta < 1.25){
				def my_beta_shift = ftof_p/(ftof_t - start_time - 2.00)/(29.9792458) 
								
  				histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
 				histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
				histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_s$ftof_sector",histoBuilders.betap).fill(p, my_beta_shift)
 				histos.computeIfAbsent("ftof_betap_pos_timeshift",histoBuilders.betap).fill(p, my_beta_shift)
 				histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
			    }
			    
*/			    
			    
			}
		    }
			

		}
		    /*def ftof_sectors = event.tof_sector[indx]
		    def ftof_t = event.tof_time[indx]
		    def ftof_p = event.tof_path[indx]
		    def ftof_paddle = event.tof_paddle[indx]
		    def ftof_layer = event.tof_layer[indx]
		    
		    def tdiff = event.tof_time[indx] - ftof_p/(29.9792458) - rf_time
		    def tdmod = tdiff % 2
		    def t_rf_corr  =  (event.tof_time[indx] - ftof_p/(29.9792458)) - event.vz[indx]/(29.9792458)
		    def rf_corr = rf_time - t_rf_corr
		     
		    //println('rf time is ' + rf_time )
		   // println( ' t_rf_corr is ' + t_rf_corr)
		   // println( ' rf_corr is ' + rf_corr)
		    
		    //histos.computeIfAbsent("ftoft_calib_time"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.calib_time).fill(tdiff,tdmod)
		    //histos.computeIfAbsent("ftoft_calib_time"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.calib_time).fill(rf_time,tdmod)		    
		    //histos.computeIfAbsent("ftof_t_rf_corr"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.time).fill(t_rf_corr)
		    //histos.computeIfAbsent("ftof_t_rf_corr_final"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.time).fill(rf_corr)

		    
		}

		if( event.tof_status.contains(index) ){
		    def ftof_sector = event.tof_sector[index]
		    def ftof_t = event.tof_time[index]
		    def ftof_p = event.tof_path[index]
		    def ftof_paddle = event.tof_paddle[index]
		    def ftof_layer = event.tof_layer[index]
		    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    		   


		    def measured_beta = beta// ftof_p/(ftof_t - start_time)/(29.9792458)
		    def measured_t = ftof_t - start_time
		    def measured_beta2 = ftof_p/(ftof_t - start_time)/(29.9792458)

 		    def calc_beta_pr = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(2212)**2)
		    def calc_time_pr = ftof_p/(calc_beta_pr*(29.9792458))

		    def calc_beta_kp = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(321)**2)
		    def calc_time_kp = ftof_p/(calc_beta_kp*(29.9792458))

		    def calc_beta_pip = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(211)**2)
		    def calc_time_pip = ftof_p/(calc_beta_pip*(29.9792458))

		    def delta_beta_pr = calc_beta_pr - measured_beta
		    def delta_t_pr = calc_time_pr - measured_t

		    def delta_beta_kp = calc_beta_kp - measured_beta
		    def delta_t_kp = calc_time_kp - measured_t

		    def delta_beta_pip = calc_beta_pip - measured_beta
		    def delta_t_pip = calc_time_pip - measured_t

		    //histos.computeIfAbsent("ftoft_pos"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
		    //histos.computeIfAbsent("ftof_t_corr_pos"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
		    //histos.computeIfAbsent("ftoft_t_pos"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_paddle,ftof_t)
		    //histos.computeIfAbsent("ftof_t_corr_pos"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_paddle,ftof_t_corr)
		    
		    

		    if( p < 1.0 && measured_beta > 0.91 ){
			
  			histos.computeIfAbsent("ftof_betap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 			histos.computeIfAbsent("ftof_betap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta)
 			histos.computeIfAbsent("ftof_betap_pos",histoBuilders.betap).fill(p, measured_beta)
 			histos.computeIfAbsent("ftof_betap_pos"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 			histos.computeIfAbsent("ftof_betap_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)

			
  			histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 			histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta2)
 			histos.computeIfAbsent("ftof_mybetap_pos",histoBuilders.betap).fill(p, measured_beta2)
 			histos.computeIfAbsent("ftof_mybetap_pos"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 			histos.computeIfAbsent("ftof_mybetap_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)

 			histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 			histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)


			///////////////////////////////////
			//  look at events with beta > 1.1 here for this run to apply a time simple shift
			if( measured_beta > 1.1 && measured_beta < 1.25){
			    def my_beta_shift = ftof_p/(ftof_t - start_time - 2.00)/(29.9792458) 
			    

  			    histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
 			    histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
			    histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_s$ftof_sector",histoBuilders.betap).fill(p, my_beta_shift)
 			    histos.computeIfAbsent("ftof_betap_pos_timeshift",histoBuilders.betap).fill(p, my_beta_shift)
 			    histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_r$run_number",histoBuilders.betap).fill(p, my_beta_shift)
			}

 			//histos.computeIfAbsent("ftof_betap_pos"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.betap).fill(p, measured_beta)
 			
			//histos.computeIfAbsent("ftof_deltab_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pr)
			//histos.computeIfAbsent("ftof_deltat_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pr)
			
 			//histos.computeIfAbsent("ftof_deltab_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_kp)
			//histos.computeIfAbsent("ftof_deltat_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_kp)
			
 			//histos.computeIfAbsent("ftof_deltab_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pip)
			//histos.computeIfAbsent("ftof_deltat_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pip)
 			//histos.computeIfAbsent("ftof_deltab_pip"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltab).fill(p, delta_beta_pip)
			//histos.computeIfAbsent("ftof_deltat_pip"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pip)

 			//histos.computeIfAbsent("ftof_betap_pos_timeshift"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.betap).fill(p, measured_beta)
 			
			//histos.computeIfAbsent("ftof_deltab_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pr)
			//histos.computeIfAbsent("ftof_deltat_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pr)
			
 			//histos.computeIfAbsent("ftof_deltab_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_kp)
			//histos.computeIfAbsent("ftof_deltat_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_kp)
			
 			//histos.computeIfAbsent("ftof_deltab_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pip)
			//histos.computeIfAbsent("ftof_deltat_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pip)
 			//histos.computeIfAbsent("ftof_deltab_pip"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltab).fill(p, delta_beta_pip)
			//histos.computeIfAbsent("ftof_deltat_pip"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pip)

		    }
		    //////////////// 
		    // per paddle 
		    //histos.computeIfAbsent("ftof_t_deltab_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta_pr)
		    //histos.computeIfAbsent("ftof_t_deltab_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta_kp)
		    //histos.computeIfAbsent("ftof_t_deltab_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta_pip)

		    //histos.computeIfAbsent("ftof_t_deltat_pr"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t_pr)
		    //histos.computeIfAbsent("ftof_t_deltat_kp"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t_kp)
		    //histos.computeIfAbsent("ftof_t_deltat_pip"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t_pip)
		    
		}
		     */	   

	    }	    
	    

	    (0..<event.npart).findAll{event.charge[it] < 0 && ( it != indx )}.findResults{index->		
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def p = lv.mag()
		def beta = event.beta[index]

		def run_number = event.run_number

		if( event.dc1_status.contains(index) ){		    
		    def hits = event.dc1[index]
		    def hit = hits.find{it.layer==12}.each{
 			histos.computeIfAbsent("dcr1hit_neg",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc2_status.contains(index) ){		    
		    def hits = event.dc2[index]
		    def hit = hits.find{it.layer==24}.each{
 			histos.computeIfAbsent("dcr2hit_neg",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc3_status.contains(index) ){		    
		    def hits = event.dc3[index]
		    def hit = hits.find{it.layer==36}.each{
 			histos.computeIfAbsent("dcr3hit_neg",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}
			

		//if( event.tof_status.contains(index) ){
		    /*def ftof_sector = event.tof_sector[index]
		    def ftof_t = event.tof_time[index]
		    def ftof_p = event.tof_path[index]
		    def ftof_paddle = event.tof_paddle[index]
		    def ftof_layer = event.tof_layer[index]
		    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)

 		    def measured_beta = beta//ftof_p/(ftof_t - start_time)/(29.9792458)
		    def measured_t = ftof_t - start_time

		    def calc_beta_km = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(-321)**2)
		    def calc_time_km = ftof_p/(calc_beta_km*(29.9792458))

		    def calc_beta_pim = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(-211)**2)
		    def calc_time_pim = ftof_p/(calc_beta_pim*(29.9792458))

		    def delta_beta_km = calc_beta_km - measured_beta
		    def delta_t_km = calc_time_km - measured_t

		    def delta_beta_pim = calc_beta_pim - measured_beta
		    def delta_t_pim = calc_time_pim - measured_t

		    //histos.computeIfAbsent("ftoft_neg"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
		    //histos.computeIfAbsent("ftof_t_corr_neg"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
		    //histos.computeIfAbsent("ftoft_t_neg"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_paddle,ftof_t)
 		    //histos.computeIfAbsent("ftof_t_corr_neg"+"_s$ftof_sector",histoBuilders.time2d).fill(ftof_paddle,ftof_t_corr)

 		    histos.computeIfAbsent("ftof_betap_neg"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta)
 		    histos.computeIfAbsent("ftof_betap_neg"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 		    histos.computeIfAbsent("ftof_betap_neg",histoBuilders.betap).fill(p, measured_beta)
 		    histos.computeIfAbsent("ftof_betap_neg"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 		    //histos.computeIfAbsent("ftof_betap_neg"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.betap).fill(p, measured_beta)

 		    //histos.computeIfAbsent("ftof_deltab_km"+"_s$ftof_sector",histoBuilders.deltab).fill(p, delta_beta_km)
		    //histos.computeIfAbsent("ftof_deltat_km"+"_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_km)

 		   // histos.computeIfAbsent("ftof_deltab_pim"+"_s$ftof_sector",histoBuilders.deltab).fill(p, delta_beta_pim)
		   // histos.computeIfAbsent("ftof_deltat_pim"+"_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pim)


		    //////////////// 
		    // per paddle 
		    //histos.computeIfAbsent("ftof_t_deltab_km"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta_km)
		    //histos.computeIfAbsent("ftof_t_deltab_pim"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta_pim)

		    //histos.computeIfAbsent("ftof_t_deltat_km"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t_km)
		    //histos.computeIfAbsent("ftof_t_deltat_pim"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t_pim)

 		    //histos.computeIfAbsent("ftof_deltab_km"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_km)
		    //histos.computeIfAbsent("ftof_deltat_km"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_km)

 		    //histos.computeIfAbsent("ftof_deltab_pim"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pim)
		    //histos.computeIfAbsent("ftof_deltat_pim"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pim)
		    
		}
*/
	   
	    }
	   
	}
    }
}
