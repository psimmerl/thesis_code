package mon
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.LinkedHashMap
import pid.electron.Electron

def beam = 10.604

//call electron cut constructor
electron = new Electron();

def myElectronCutStrategies = [
    electron.passElectronStatus,
    electron.passElectronChargeCut,
    electron.passElectronEBPIDCut,
    electron.passElectronNpheCut,
    electron.passElectronVertexCut,
    electron.passElectronPCALFiducialCut,
    electron.passElectronEIEOCut,
    electron.passElectronDCR1
]

reqBanks = [ 
    part : "REC::Particle",
    ec   : "REC::Calorimeter",
    cc   : "REC::Cherenkov",    
    traj : "REC::Traj",
    trck : "REC::Track"
]

def out = new TDirectory()

def histos = new ConcurrentHashMap()
histoBuilders = [

    ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
    phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
    vzphi : { title -> new H2F("$title","$title", 140, 0, 140, 1000, -100, 100 ) },
    vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
    hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
    ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.5) },
    nphe : { title -> new H1F("$title","$title", 40, 0.0, 40 ) },
    eieo : { title -> new H2F("$title","$title", 1000, 1.0, 1000, 1.0 ) },
    vz   : { title -> new H1F("$title","$title", 400, -40, 40) },
    w    : { title -> new H1F("$title","$title", 500, 0.5, beam+0.3 ) },
    q2w  : { title -> new H2F("$title","$title", 500, 0.5, beam+0.3, 500, 0, 2 * beam ) },  
    passrate : {title -> new H1F("$title","$title", 8, 0.0, 8.0 ) },
    wtheta : { title -> new H2F("$title","$title", 100, 0.0, 2.5, 100, 0.0, 20.0 ) }

]

for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    def cc = 0
    
    while(reader.hasEvent() && cc < 20000) {
	
	def event = reader.getNextEvent()

	def banks  = reqBanks.collect{name, bank -> [name, event.getBank(bank)] }.collectEntries()
	def banks_pres = reqBanks.every{name, bank -> event.hasBank(bank) }

	
	if (banks_pres) { 
	    cc+=1
	    
	    //loop over only negative tracks first
	    def my_el_cuts = (0..<banks.part.rows()).findAll{banks.part.getByte('charge',it)<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(banks,ii) } ] }.collectEntries()
	
 	    //def ical = banks.ec.getShort("pindex").findIndex{ it == ipt }
            //particle.sector = banks.ec.getByte("sector", ical)
		
	    //key is index in rec particle bank, value is the list of booleans on whether track passed ith cut
	    /*
	    my_el_cuts.each{ key, value -> 
		(0..<banks.cc.rows()).find{ ii -> banks.cc.getInt("pindex", ii) == key }?.with{ inphe -> 
		    def nphe = banks.cc.getFloat("nphe", inphe) }?.with{ nphe-> value.eachWithIndex{ val, indx ->
		    if(val) histos.computeIfAbset("el_nphe_pass_cut$indx",histoBuilders.nphe).fill(nphe) } } }
	   */

	    /*
	    my_el_cuts.each{ key, value ->
		banks.cc.getShort("pindex").find{it==key}?.with{ icc->
		    def nphe = banks.cc.getFloat("nphe",icc)
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){			    
			    println("passed cut $cut_index and nphe is $nphe")
			    histos.computeIfAbsent("el_nphe_pass_cut$cut_index",histoBuilders.nphe).fill(nphe)
			}
		    }
		}
	    }
	     */


	    //other method to load banks before looping over particles
	    cher = [banks.cc.getShort('pindex'), banks.cc.getFloat('nphe')].transpose().collectEntries()	   
	    e_pcal = (0..<banks.ec.rows()).findResults{(banks.ec.getByte('detector',it) == DetectorType.ECAL.getDetectorId() && 
							banks.ec.getInt('layer',it) == 1) ? it : null }.collectEntries{ index -> [(banks.ec.getShort('pindex',index)), banks.ec.getFloat('energy',index)]} 
	    
	    ei_ecal = (0..<banks.ec.rows()).findResults{(banks.ec.getByte('detector',it) == DetectorType.ECAL.getDetectorId() && 
							 banks.ec.getInt('layer',it) == 4) ? it : null }.collectEntries{ index -> [banks.ec.getShort('pindex',index), banks.ec.getFloat('energy',index)]} 
	    
	    eo_ecal = (0..<banks.ec.rows()).findResults{(banks.ec.getByte('detector',it) == DetectorType.ECAL.getDetectorId() && 
							 banks.ec.getInt('layer',it) == 7) ? it : null }.collectEntries{ index -> [banks.ec.getShort('pindex',index), banks.ec.getFloat('energy',index)]} 

	    //findResults returns list of the list indices satisfying closure. 
	    dc_r1_hit = (0..<banks.traj.rows()).findResults{(banks.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() && 
							     banks.traj.getInt('layer',it) == 12) ? it : null }.collectEntries{index -> [(banks.traj.getShort('pindex',index)), 'xy'.collect{ axis -> banks.traj.getFloat(axis,index)}]} 

	    dc_r2_hit = (0..<banks.traj.rows()).findResults{(banks.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() && 
							     banks.traj.getInt('layer',it) == 24) ? it : null }.collectEntries{index -> [banks.traj.getShort('pindex',index), 'xy'.collect{ axis -> banks.traj.getFloat(axis,index)}]} 

	    dc_r3_hit = (0..<banks.traj.rows()).findResults{(banks.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() && 
							     banks.traj.getInt('layer',it) == 36) ? it : null }.collectEntries{index -> [banks.traj.getShort('pindex',index), 'xy'.collect{ axis -> banks.traj.getFloat(axis,index)}]} 


	    my_el_cuts.each{ key, value ->
		def lv = new Vector3(*['px','py','pz'].collect{banks.part.getFloat(it,key)})
		def p = lv.mag()
		def vz = banks.part.getFloat('vz',key)
		def theta = Math.toDegrees(lv.theta())
		def phi = Math.toDegrees(lv.phi())

		//println("key type " )
 		//println(key.getClass())

		//println("$key in ")
		///println(e_pcal.keySet())
		//println(key in e_pcal.keySet())

		cher.find{
		    if( it.key == key ){ 
			//plot nphe for each cut
			value.eachWithIndex{ cut_result, cut_index ->
			    if(cut_result){histos.computeIfAbsent("nphe_cut$cut_index",histoBuilders.nphe).fill(it.value)}					    
			}
			//plot nphe for tracks passing all cuts
			if(!value.contains(false)){histos.computeIfAbsent("nphe_passall",histoBuilders.nphe).fill(it.value)}
		    }
		}
				
				
		e_pcal.find{ //containsKey(key)){
		    if( it.key == key ){
			value.eachWithIndex{ cut_result, cut_index ->
			    if(cut_result){histos.computeIfAbsent("pcalsf_cut$cut_index",histoBuilders.ecsf).fill(p, it.value/p)}
			}
			if(!value.contains(false)){histos.computeIfAbsent("pcalsf_passall",histoBuilders.ecsf).fill(p, it.value/p)}		    
		    }	
		}
		
		if( ei_ecal.containsKey(key)  && eo_ecal.containsKey(key) ){			
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){histos.computeIfAbsent("eieo_cut$cut_index",histoBuilders.eieo).fill(ei_ecal.get(key), eo_ecal.get(key))}
		    }
		    if(!value.contains(false)){histos.computeIfAbsent("eieo_passall",histoBuilders.eieo).fill(ei_ecal.get(key), eo_ecal.get(key))}		    
		}
		
		
		//println("key is $key")
		//println(dc_r1_hit)
		dc_r1_hit.find{
		    if( it.key == key ){ //dc_r1_hit.containsKey(key) ){		
			value.eachWithIndex{ cut_result, cut_index ->
			    if( cut_result ){ histos.computeIfAbsent("dcr1hit_cut$cut_index",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr1hit_passall",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}			
		    }
		}

		dc_r2_hit.find{
		    if( it.key == key ){ //dc_r2_hit.containsKey(key) ){		
			value.eachWithIndex{ cut_result, cut_index ->
			    if( cut_result ){ histos.computeIfAbsent("dcr2hit_cut$cut_index",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr2hit_passall",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}			
		    }
		}

		dc_r3_hit.find{
		    if( it.key == key ){ //dc_r3_hit.containsKey(key) ){		
			value.eachWithIndex{ cut_result, cut_index ->
			    if( cut_result ){ histos.computeIfAbsent("dcr3hit_cut$cut_index",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr3hit_passall",histoBuilders.hitxy).fill(it.value.get(0), it.value.get(1))}			
		    }
		}

		/*
		if( dc_r2_hit.containsKey(key) ){
		    value.eachWithIndex{ cut_result, cut_index ->
			if( cut_result ){ histos.computeIfAbsent("dcr2hit_cut$cut_index",histoBuilders.hitxy).fill(dc_r2_hit.get(key).get(0), dc_r2_hit.get(key).get(1))}
		    }
		    if(!value.contains(false)){histos.computeIfAbsent("dcr2hit_passall",histoBuilders.hitxy).fill(dc_r2_hit.get(key).get(0), dc_r2_hit.get(key).get(1))}			
		}

		if( dc_r3_hit.containsKey(key) ){
		    value.eachWithIndex{ cut_result, cut_index ->
			if( cut_result ){ histos.computeIfAbsent("dcr3hit_cut$cut_index",histoBuilders.hitxy).fill(dc_r3_hit.get(key).get(0), dc_r3_hit.get(key).get(1))}
		    }
		    if(!value.contains(false)){histos.computeIfAbsent("dcr3hit_passall",histoBuilders.hitxy).fill(dc_r3_hit.get(key).get(0), dc_r3_hit.get(key).get(1))}			
		}
		*/
		value.eachWithIndex{ cut_result, cut_index ->
		    histos.computeIfAbsent("vz_cut$cut_index",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_cut$cut_index",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_cut$cut_index",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_cut$cut_index",histoBuilders.vzphi).fill(vz, phi)
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("vz_passall",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_passall",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_passall",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_passall",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall",histoBuilders.vzphi).fill(vz, phi)
		}
	    		   
	    }   	    	   	    
	}
    }

    out.mkdir("/test/")
    out.cd("/test/")
    histos.values().each{out.addDataSet(it)}
    out.writeFile("electron_pid.hipo")

    reader.close();


}
