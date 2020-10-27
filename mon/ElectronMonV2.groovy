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


class ElectronMonV2{

    def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def rfBucketLength = 4.0 //ns
    
    //call electron cut constructor
    def field_config = ''
    def el_cut_strictness = null
    def electron = new ElectronFromEvent();

    def sf_jlab = [1: [1] ]
    
    def sigma_range = 3.5// to match rga note 08/2020 before it was 4.0
    def p0mean = [0.116405, 0.118783, 0.115558, 0.112439, 0.120496, 0.117107]
    def p1mean = [-0.0155433, 0.109417, -0.0088641, -0.0704795, 0.100674, -0.0164561]
    def p2mean = [0.00722919, 0.0064115, 0.00773618, 0.00769267, 0.00393716, 0.00702487]
    def p3mean = [-0.000859442, -0.000704544, -0.000853684, -0.00073044, -0.000403603, -0.000825457]
    def p0sigma = [0.0164086, 0.00448347, 0.00770184, 0.00412631, 0.00687372, 0.00747944]
    def p1sigma = [0.00771924, 0.020591, 0.0164028, 0.0217636, 0.0183004, 0.0184456]
    def p2sigma = [-0.00125094, 0.000362296, -0.000381947, 0.000232047, -0.000439471, -0.000560982]
    def p3sigma = [8.43022e-05, -2.94663e-05, 2.47408e-05, -5.70412e-06, 5.14683e-05, 5.80948e-05]

    

    def setCutParameters(){
	println(' [epkpkm::setCutParameters] -> field configuration is ' + field_config)
	println(' [epkpkm::setCutParameters] -> setting cut values for field ' + field_config)
	electron.setElectronCutParameters(field_config)
    }
    
    def setCutStrictness(){
	println(' [epkpkm::setCutStrictness] -> cut strictness level is ' + el_cut_strictness )
	electron.setElectronCutStrictness(el_cut_strictness)
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
        
    def myElectronCutStrategies = [
	electron.passElectronChargeCut,
	electron.passElectronEBPIDCut,
	electron.passElectronStatus,
	electron.passElectronMinMomentum,
	electron.passElectronEIEOCut,
	electron.passElectronNpheCut,
	electron.passElectronVertexCut,
	electron.passElectronPCALFiducialCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3,
	electron.passElectronSamplingFractionCut
 	//electron.passElectronTrackQualityCut

    ]

    def charge_cut_index=0
    def eb_el_pid_cut_index=1
    def status_cut_index=2
    def min_mntm_cut_index=3
    def eieo_cut_index=4
    def nphe_cut_index=5
    def vz_cut_index=6
    def pcalfid_cut_index=7
    //def dcr1_cut_index=6
    //def dcr2_cut_index=7
    //def dcr3_cut_index=8
    def ecsf_cut_index=8


    
    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 1000, -450.0, 450.0, 1000, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.5) },
	ecsf_eb : { title -> new H2F("$title","$title", 500, 0.0, 2.5, 500, 0.0, 0.4) },
	ecsf_chi : { title -> new H1F("$title","$title", 100, -5.0, 5.0) },	
	calsf : { title -> new H2F("$title","$title", 500, 0.0, 0.35, 500, 0.0, 0.35) },
	
	calsf_pion_remover : { title -> new H2F("$title","$title", 500, 0.0, 0.3, 500, 0.0, 0.35) },

	uvwsf : { title -> new H2F("$title","$title", 1000, 0.0, 500.0, 1000, 0.0, 0.36) },
	uvwhit : { title -> new H1F("$title","$title", 500, 0.0, 100.0) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 75, 0.0, 75 ) },
	nphe : { title -> new H1F("$title","$title", 75, 0.0, 75.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 2.0, 1000, 0.0, 2.0 ) },
 	eieo_small : { title -> new H2F("$title","$title", 1000, 0.0, 0.4, 1000, 0.0, 0.4 ) },

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

    
    def passElectronSamplingFractionCut( event, index ){
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
		def mean = p0mean[sector] * (1 + p/Math.sqrt(p*p + p1mean[sector])) + p2mean[sector]*p + p3mean[sector]*p*p
		def sigma = p0sigma[sector] + p1sigma[sector]/Math.sqrt(p) + p2sigma[sector]*p + p3sigma[sector]*p*p

		def upper_cut = mean + sigma_range * sigma
		def lower_cut = mean - sigma_range * sigma
		
		if( edep/p <= upper_cut && edep/p >= lower_cut ) return true	    
	    }
	}
	return false
    }



    def processEvent( event ){
	
	//println("[ElectronMon::processEvent] - processing event now")
	def run_number = event.run_number
	//def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	//println(my_el_cuts)
	
	//my_el_cuts.each{ index, value ->
	for( index in 0..<event.npart ){
	    if( event.charge[index]>=0 ) continue
	    if( event.charge[index] < 0 ){

		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def p = lv.mag()
		def vz = event.vz[index]
		def theta = Math.toDegrees(lv.theta())
		def phi = Math.toDegrees(lv.phi())
		
		//passing statistics
		//value.eachWithIndex{ cut_res, cut_ind -> cut_res ? histos.computeIfAbsent("passResults",histoBuilders.passrate).fill(cut_ind) : null }
		/*
		 //plot cherenkov number of photons
		 if (event.cherenkov_status.contains(index)) { 
		 def nphe = event.nphe[index]	    
		 def cherenkov_sector = event.cherenkov_sector[index]
		 value.eachWithIndex{ cut_result, cut_index ->
		 // look at nphe for all negative particles
		 if(cut_result && cut_index==0){
		 histos.computeIfAbsent("nphe_cut$cut_index",histoBuilders.nphe).fill(nphe)
		 histos.computeIfAbsent("nphe_s"+cherenkov_sector+"_cut$cut_index",histoBuilders.nphe).fill(nphe)
	     }
		 if(cut_result && cut_index==eb_el_pid_cut_index){
		 histos.computeIfAbsent("nphe_cut$cut_index",histoBuilders.nphe).fill(nphe)
		 histos.computeIfAbsent("nphe_s"+cherenkov_sector+"_cut$cut_index",histoBuilders.nphe).fill(nphe)
	     }
	     }
		 //plot nphe for tracks passing all cuts
		 if( !value.contains(false) ){
		 histos.computeIfAbsent("nphe_passall",histoBuilders.nphe).fill(nphe)
		 histos.computeIfAbsent("nphe_passall_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  
	     }		    
		 
		 // check nphe with all cuts but this one applied
		 if( everythingElsePassed(value,nphe_cut_index) ){
		 histos.computeIfAbsent("nphe_pass_all_but_nphe",histoBuilders.nphe).fill(nphe)
		 histos.computeIfAbsent("nphe_pass_all_but_nphe_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  		    		   
	     }
		 
	     }
		 
		 

		 //plot PCAL sampling fraction distributions
		 if( event.pcal_status.contains(index) ){
		 def e_pcal = event.pcal_energy[index]
		 def pcal_sect = event.pcal_sector[index]
		 def pcal_u = event.pcal_u[index]
		 def pcal_v = event.pcal_v[index]
		 def pcal_w = event.pcal_w[index]
		 def pcal_x = event.pcal_x[index]
		 def pcal_y = event.pcal_y[index]
		 def pcal_z = event.pcal_z[index]
		 def pcal_t = event.pcal_time[index]
		 def pcal_l = event.pcal_path[index]
		 value.eachWithIndex{ cut_result, cut_index ->
		 //pcal UVW position vs pcal/P - tells use where to cut
		 if(cut_result && cut_index==0){
		 histos.computeIfAbsent("pcalsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 
		 histos.computeIfAbsent("pcalusf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		 histos.computeIfAbsent("pcalvsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		 histos.computeIfAbsent("pcalwsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
		 histos.computeIfAbsent("pcalxy_cut$cut_index"+"_s$pcal_sect",histoBuilders.hitxy).fill(pcal_x, pcal_y)			
 		 histos.computeIfAbsent("pcalxy_cut$cut_index",histoBuilders.hitxy).fill(pcal_x, pcal_y)			
	     }

		 if(cut_result && cut_index==eb_el_pid_cut_index){
		 histos.computeIfAbsent("pcalsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 
		 histos.computeIfAbsent("pcalusf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		 histos.computeIfAbsent("pcalvsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		 histos.computeIfAbsent("pcalwsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
		 histos.computeIfAbsent("pcalxy_cut$cut_index"+"_s$pcal_sect",histoBuilders.hitxy).fill(pcal_x, pcal_y)
		 histos.computeIfAbsent("pcalxy_cut$cut_index",histoBuilders.hitxy).fill(pcal_x, pcal_y)
	     }
	     }
		 //plot pcal UVW hit positions vs PCAL / P
		 if(!value.contains(false)){
		 //histos.computeIfAbsent("pcalsf_passall"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		 histos.computeIfAbsent("pcalusf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		 histos.computeIfAbsent("pcalvsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
 		 histos.computeIfAbsent("pcalwsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
		 histos.computeIfAbsent("pcalxy_passall_s$pcal_sect",histoBuilders.hitxy).fill(pcal_x, pcal_y)
		 histos.computeIfAbsent("pcalxy_passall",histoBuilders.hitxy).fill(pcal_x, pcal_y)
	     }		

		 //now check dist without cutting on pcal fiducial edge but with all other cuts in place
		 if( everythingElsePassed(value,pcalfid_cut_index) ){
 		 histos.computeIfAbsent("pcalusf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		 histos.computeIfAbsent("pcalvsf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		 histos.computeIfAbsent("pcalwsf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)		    
		 histos.computeIfAbsent("pcalxy_passall_but_pcalfid_s$pcal_sect",histoBuilders.hitxy).fill(pcal_x, pcal_y)
		 histos.computeIfAbsent("pcalxy_passall_but_pcalfid",histoBuilders.hitxy).fill(pcal_x, pcal_y)
	     }		

	     }

		 //pcal vs ecal
		 // check the cut on the total sampling fraction
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
		 def pcal_sect = event.pcal_sector[index]
		 //manually calculate the chi2 for sampling fraction based on sector, p and edep/p
		 def sf_chi = KinTool.getSamplingFractionFromMomentum(sector, p , edep/p)
		 value.eachWithIndex{ cut_result, cut_index ->
		 //sampling fraction for all negatives 
		 if(cut_result && cut_index==0){		
		 //println(pcaldep)
		 //println(eodep + eidep)
		 histos.computeIfAbsent("calsfchi_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 histos.computeIfAbsent("calsfchi_cut$cut_index"+"_s$sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		 //traditional ecal/p distributions
		 histos.computeIfAbsent("calsf_cut$cut_index"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_cut$cut_index"+"_s$sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 	 histos.computeIfAbsent("pcalvsecal_cut$cut_index",histoBuilders.eieo).fill(pcaldep, eidep + eodep)

		 
		 histos.computeIfAbsent("calsf_pcalvsecal_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		 
		 if( p > 4.0 ){
		 histos.computeIfAbsent("calsf_pcalvsecal_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 histos.computeIfAbsent("calsfchi_bigp_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 histos.computeIfAbsent("calsfchi_bigp_cut$cut_index"+"_s$sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_bigp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
	     }
		 else if ( p < 4.0 ){
		 histos.computeIfAbsent("calsfchi_smallp_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 histos.computeIfAbsent("calsfchi_smallp_cut$cut_index"+"_s$sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		 histos.computeIfAbsent("calsf_pcalvsecal_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_lowp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
	     }			    
	     }
		 //check what only applying EB PID does
		 if(cut_result && cut_index==eb_el_pid_cut_index){
		 histos.computeIfAbsent("calsfchi_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 histos.computeIfAbsent("calsfchi_cut$cut_index"+"_s$sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		 //traditional ecal/p distributions
		 histos.computeIfAbsent("calsf_cut$cut_index"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_cut$cut_index"+"_s$sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 	 histos.computeIfAbsent("pcalvsecal_cut$cut_index",histoBuilders.eieo).fill(pcaldep, eidep + eodep)

		 
		 histos.computeIfAbsent("calsf_pcalvsecal_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
	     }
	     }

		 if(!value.contains(false)){
		 histos.computeIfAbsent("calsf_passall"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_passall"+"_s$sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 	 histos.computeIfAbsent("pcalvsecal_passall",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
		 histos.computeIfAbsent("calsfchi_passall",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 histos.computeIfAbsent("calsfchi_passall_s$sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		 histos.computeIfAbsent("calsf_pcalvsecal_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		 
		 if( p > 4.0 ){
		 histos.computeIfAbsent("calsf_pcalvsecal_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_bigp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
	     }
		 else if ( p < 4.0 ){
		 histos.computeIfAbsent("calsf_pcalvsecal_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_lowp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
	     }


	     }

		 if( everythingElsePassed(value,eieo_cut_index) ){
		 histos.computeIfAbsent("calsf_passall_but_eieo"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_passall_but_eieo_"+"_s$sector",histoBuilders.eieo).fill(eidep + eodep,pcaldep)
	 	 histos.computeIfAbsent("pcalvsecal_passall_but_eieo",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
		 histos.computeIfAbsent("calsf_pcalvsecal_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)

		 histos.computeIfAbsent("calsf_pcalvsecal_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_eieo",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, eodep/p)

	     }


		 if( everythingElsePassed(value,ecsf_cut_index) ){
		 //pass all SF cut, see why we need to apply it
		 histos.computeIfAbsent("calsf_passall_but_ecsf"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_passall_but_ecsf_"+"_s$sector",histoBuilders.eieo).fill(eidep + eodep,pcaldep)
	 	 histos.computeIfAbsent("pcalvsecal_passall_but_ecsf",histoBuilders.eieo).fill(eidep + eodep, pcaldep)

		 histos.computeIfAbsent("calsf_pcalvsecal_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_ecsf",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, eodep/p)

	     }


		 if( everythingElsePassed(value,min_mntm_cut_index) ){
		 //pass all but low mntm electrons - check SF -> why we need it
		 histos.computeIfAbsent("cal_sf_passall_but_minmntm"+"_s$sector",histoBuilders.ecsf).fill(p, edep/p) 	
		 histos.computeIfAbsent("pcalvsecal_passall_but_minmntm_"+"_s$sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 	 histos.computeIfAbsent("pcalvsecal_passall_but_minmntm",histoBuilders.eieo).fill(eidep + eodep, pcaldep)

		 histos.computeIfAbsent("calsf_pcalvsecal_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		 if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_minmntm",histoBuilders.calsf).fill(eidep/p, eodep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		 if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, eodep/p)	
	     }

	     }
		 
	     }

		 // EI vs EO energy deposited cut
		 if( event.ecal_inner_status.contains(index) &&  event.ecal_outer_status.contains(index)){
		 def e_eical = event.ecal_inner_energy[index]
 		 def e_eocal = event.ecal_outer_energy[index]
		 def ei_sect = event.ecal_inner_sector[index]
		 def eo_sect = event.ecal_outer_sector[index]
		 
		 //make sure sector is the same!
		 if( ei_sect == eo_sect ){			
		 value.eachWithIndex{ cut_result, cut_index ->
		 // EI vs EO for all negative tracks
		 if(cut_result && cut_index==0){
		 histos.computeIfAbsent("eieo_cut$cut_index",histoBuilders.eieo).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_small_cut$cut_index",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_small_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
		 if( e_eocal > 0 && e_eical > 0 ) histos.computeIfAbsent("ecalsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
	     }
		 if(cut_result && cut_index==eb_el_pid_cut_index){
		 histos.computeIfAbsent("eieo_cut$cut_index",histoBuilders.eieo).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_small_cut$cut_index",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
		 histos.computeIfAbsent("eieo_small_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
		 if( e_eocal > 0 && e_eical > 0 ) histos.computeIfAbsent("ecalsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
	     }



	     }		    
		 // EI vs EO passing all cuts
		 if(!value.contains(false)){
		 histos.computeIfAbsent("eieo_passall",histoBuilders.eieo).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_passall"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_small_passall",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_small_passall"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		 if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent("ecalsf_eivseo_passall",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
	     }		   
		 // check ei vs eo without passing all cuts
		 if( everythingElsePassed(value,eieo_cut_index) ){
		 histos.computeIfAbsent("eieo_passall_but_eieo",histoBuilders.eieo).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_passall_but_eieo"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_small_passall_but_eieo",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		 histos.computeIfAbsent("eieo_small_passall_but_eieo"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		 if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent("ecalsf_eivseo_passall_but_eieo",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
	     }		    
		 
	     }		    
	     }
		 
		 */
		/// these are specialty plots for the rga analysis note
		// first check the SF distribution with Nphe and PCAL edep cut
		
		def pass_eb_ele_or_pim =  ((event.status[index] < 0 && event.pid[index] == 11 ) || (event.pid[index] == -211 ))
 		def pass_sf = passElectronSamplingFractionCut( event, index )
		

		if (event.cherenkov_status.contains(index) &&  event.pcal_status.contains(index) && 
		    (event.ecal_inner_status.contains(index) ||  event.ecal_outer_status.contains(index) ) ){
		    
		    def nphe = event.nphe[index]	    
		    def cherenkov_sector = event.cherenkov_sector[index]
		    
		    def eidep = 0
		    def eodep = 0
		    def pcaldep = 0
		    def sector = -1
		    if( event.ecal_inner_status.contains(index) ){
			eidep = event.ecal_inner_energy[index]
		    }
		    if( event.ecal_outer_status.contains(index) ){
			eodep = event.ecal_outer_energy[index]
		    }
		    if( event.pcal_status.contains(index) ){
			pcaldep = event.pcal_energy[index]
			sector = event.pcal_sector[index]-1
		    }
		    // now cut on the nphe and the pcal energy
		    def sf_edep_p = (eidep + eodep + pcaldep)/p
		    def edep = eidep + eodep + pcaldep


		    histos.computeIfAbsent("calsf_pion_remover_all_sect",histoBuilders.calsf_pion_remover).fill(eidep/p, pcaldep/p)
		    def pass_pcal_diagonal = false
		    if( pass_sf ){
			if( p > 4.5 ){
			    pass_pcal_diagonal = eidep/p < 0.2 - pcaldep/p
			}
			else{
			    pass_pcal_diagonal=true
			}
		    }
		
		    if( pass_pcal_diagonal ){
			histos.computeIfAbsent("calsf_pion_remover_pass_p_pass_sf_pass_diag_all_sect",histoBuilders.calsf_pion_remover).fill(eidep/p, pcaldep/p)
		    }
			


		    histos.computeIfAbsent("sampling_fraction_all_sect",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
		    histos.computeIfAbsent("sampling_fraction_all_sect${sector}",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
		    
 		    histos.computeIfAbsent("sampling_fraction_all_sect_nom",histoBuilders.ecsf).fill(p, sf_edep_p)
		    histos.computeIfAbsent("sampling_fraction_all_sect${sector}_nom",histoBuilders.ecsf).fill(p, sf_edep_p)		    

		    histos.computeIfAbsent("pcalvsecal_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
		    histos.computeIfAbsent("nphe_all_sect",histoBuilders.nphe).fill(nphe)
		    
		    if( nphe > 2 && pcaldep > 0.06 ){
			// get sf as edep / momentum 
			// vs edep ( since this is how jlab does it )
			histos.computeIfAbsent("sampling_fraction_pass_nphe_pass_pcaldep_all_sect",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
			histos.computeIfAbsent("sampling_fraction_pass_nphe_pass_pcaldep_sect${sector}",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
			
 			histos.computeIfAbsent("sampling_fraction_pass_nphe_pass_pcaldep_all_sect_nom",histoBuilders.ecsf).fill(p, sf_edep_p)
			histos.computeIfAbsent("sampling_fraction_pass_nphe_pass_pcaldep_sect${sector}_nom",histoBuilders.ecsf).fill(p, sf_edep_p)		    
		    }		
		    else{
			histos.computeIfAbsent("sampling_fraction_fail_at_least_nphe_pcaldep_all_sect",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
			histos.computeIfAbsent("sampling_fraction_fail_at_least_nphe_pcaldep_sect${sector}",histoBuilders.ecsf_eb).fill(edep, sf_edep_p)
			
 			histos.computeIfAbsent("sampling_fraction_fail_at_least_nphe_pcaldep_all_sect_nom",histoBuilders.ecsf).fill(p, sf_edep_p)
			histos.computeIfAbsent("sampling_fraction_fail_at_least_nphe_pcaldep_sect${sector}_nom",histoBuilders.ecsf).fill(p, sf_edep_p)		    
		    }

		    if( Math.abs(event.chi2pid[index]) < 5 && nphe > 2 ){
			histos.computeIfAbsent("pcalvsecal_pass_nphe_pass_calsf_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			histos.computeIfAbsent("pcalvsecal_pass_nphe_pass_calsf_all_sect_small",histoBuilders.eieo_small).fill(pcaldep, eidep + eodep)
		    }
		    else{
			histos.computeIfAbsent("pcalvsecal_fail_at_least_nphe_calsf_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			histos.computeIfAbsent("pcalvsecal_fail_at_least_nphe_calsf_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)

		    }

		    if( Math.abs(event.chi2pid[index]) < 5 && pcaldep > 0.06 ){
			histos.computeIfAbsent("nphe_pass_calsf_pass_pcaldep_all_sect",histoBuilders.nphe).fill(nphe)
		    }	
		    else{
			histos.computeIfAbsent("nphe_fail_at_least_calsf_pcaldel_all_sect",histoBuilders.nphe).fill(nphe)
		    }

		    ////////////////////
		    // modified pid using function for sf
		    if( pass_sf && nphe > 2 ){
			histos.computeIfAbsent("pcalvsecal_pass_nphe_pass_modcalsf_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
		    }
		    else{
			histos.computeIfAbsent("pcalvsecal_fail_at_least_nphe_modcalsf_all_sect",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
		    }

		    if( pass_sf && pcaldep > 0.06 ){
			histos.computeIfAbsent("nphe_pass_modcalsf_pass_pcaldep_all_sect",histoBuilders.nphe).fill(nphe)
		    }	
		    else{
			histos.computeIfAbsent("nphe_fail_at_least_modcalsf_pcaldel_all_sect",histoBuilders.nphe).fill(nphe)
		    }

		
		}


		

		// check the distributions to justify fiducial volume cuts
		/*
 		 if( event.dc1_status.contains(index) ){		    
		 def hits = event.dc1[index]
		 def hit = hits.find{it.layer==12}.each{			
		 value.eachWithIndex{ cut_result, cut_index ->
		 if(cut_result && cut_index==0){
  		 if( cut_result ){ histos.computeIfAbsent("dcr1hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
	     }
	     }		    
		 if(!value.contains(false)){histos.computeIfAbsent("dcr1hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		 
		 if( everythingElsePassed(value,dcr1_cut_index) ){
		 histos.computeIfAbsent("dcr1hit_pass_all_but_dcr1cut",histoBuilders.hitxy).fill(it.x, it.y) 
	     }		
		 
	     }
		 
	     }
 		 
		 if( event.dc2_status.contains(index) ){		    
		 def hits = event.dc2[index]
		 def hit = hits.find{it.layer==24}.each{
		 value.eachWithIndex{ cut_result, cut_index ->
		 //println(cut_result + ' ' + cut_index)
		 if(cut_result && cut_index==0){
		 if( cut_result ){ histos.computeIfAbsent("dcr2hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y)}
	     }
		 if(!value.contains(false)){histos.computeIfAbsent("dcr2hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}		
	     }
		 if( everythingElsePassed(value,dcr2_cut_index) ){
		 histos.computeIfAbsent("dcr2hit_pass_all_but_dcr2cut",histoBuilders.hitxy).fill(it.x, it.y)
	     }
	     }
	     }
		 
		 if( event.dc3_status.contains(index) ){		    
		 def hits = event.dc3[index]
		 def hit = hits.find{it.layer==36}.each{
		 value.eachWithIndex{ cut_result, cut_index ->
		 if(cut_result && cut_index==0){
		 if( cut_result ){ histos.computeIfAbsent("dcr3hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
	     }
		 if(!value.contains(false)){histos.computeIfAbsent("dcr3hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
	     }
		 if( everythingElsePassed(value,dcr3_cut_index) ){
		 histos.computeIfAbsent("dcr3hit_pass_all_but_dcr3cut",histoBuilders.hitxy).fill(it.x, it.y)
	     }
	     }
	     }
		 */
	    	
		//generic p, theta, phi, vx, vy, vz information
		//also CHI2 and NDF
		if( event.npart>0 && event.dc1_status.contains(index) && event.pid[index] == 11 ) {
		    def dc_sector  = event.dc_sector[index]		    
		    def dc_chi2 = event.dc_chi2[index]
		    def dc_ndf = event.dc_ndf[index]
		    def vx = event.vx[index]
		    def vy = event.vy[index]
		    def measured_el = LorentzVector.withPID(11,event.px[index],event.py[index],event.pz[index])
		    def delta_energy = KinTool.delta_meas_energy(beam, measured_el) 
		    def delta_theta = KinTool.delta_meas_theta(beam, measured_el)
		    histos.computeIfAbsent("vz_cut_$dc_sector",histoBuilders.vz).fill(vz)

		    /*
 		    value.eachWithIndex{ cut_result, cut_index ->
			if( cut_result && cut_index==0 ){
			    histos.computeIfAbsent("ptheta_cut$cut_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			    histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
			    histos.computeIfAbsent("vztheta_cut$cut_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			    histos.computeIfAbsent("vxy_cut$cut_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			    histos.computeIfAbsent("vzphi_cut$cut_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			    histos.computeIfAbsent("chi_cut$cut_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 			    histos.computeIfAbsent("ndf_cut$cut_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    
			}

			if(cut_result && cut_index==eb_el_pid_cut_index){
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

		    if( value[status_cut_index] && value[eb_el_pid_cut_index]){
			histos.computeIfAbsent("vz_cut_pass_eb_pass_stat"+"_s$dc_sector",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent("ptheta_cut_eb_pass_stat"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			histos.computeIfAbsent("phitheta_cut_eb_pass_stat",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent("vztheta_cut_eb_pass_stat"+"s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent("vxy_cut_eb_pass_stat"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent("vzphi_cut_eb_pass_stat"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			histos.computeIfAbsent("chi_cut_eb_pass_stat"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 			histos.computeIfAbsent("ndf_cut_eb_pass_stat"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    		    
		    }
		    */
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
}
