#!/opt/coatjava/6.3.1/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import java.util.LinkedHashMap
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import utils.KinTool

LorentzVector.metaClass.static.withPID = {pid,x,y,z->
    double m = PDGDatabase.getParticleMass(pid)
    double e = Math.sqrt(x*x+y*y+z*z+m*m)
    return new LorentzVector(x,y,z,e)
}

LorentzVector.metaClass.plus = {
    LorentzVector a = new LorentzVector(delegate)
    a.add(it)
    return a
}

LorentzVector.metaClass.minus = {
    LorentzVector a = new LorentzVector(delegate)
    a.sub(it)
    return a
}

def beam = 10.604

def histos = new ConcurrentHashMap()
def histoBuilders = [
    
    deltat : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) },
    deltab : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -1, 1  ) },
    betap : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 1.4  ) },
    deltab_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -1, 1  ) },
    deltat_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -2, 2  ) },
    calib_time : { title -> new H2F("$title","$title", 1000, -100.0, 100, 1000, -100, 100  ) },
    betacomp : { title -> new H2F("$title","$title", 100, 0.80, 1.2, 100, 0.80, 1.2  ) }

]

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }
	    
	    //Use converter to get data in the Event Class structure
            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)
	    
	    // Reconstructed (for data and simulation)
	    //
	    (0..<event.npart).find {
		event.pid[it] == 11 && event.status[it] < 0 && event.p[it] > 1.5
            }?.each { idx ->
	    	
		def ele = LorentzVector.withPID(11, event.px[idx], event.py[idx], event.pz[idx]) 	    	    	    	   
		//look at positive charged hadrons
		
		// status :  FT (1000 <= status < 2000), FD (2000 <= status < 3000) or CD (status >= 40000) 
		(0..<event.npart).findAll{event.charge[it] > 0 && ( it != idx ) && event.status[it] >= 2000 && event.status[it] < 4000 }.findResults{index->		
		    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		    def pid = event.pid[index]
 		    def p = lv.mag()
		    def beta = event.beta[index]
              	    def start_time_rf_corrected = event.vt[index] // start time with vertex correction
		    def vt = event.vt[index]
 		    def run_number = event.run_number
		    
		    if(event.tof_status.contains(index) ){		  
			def layer_to_use = -1
			def layer_avail = event.tof[index].layer
			def my_layer = []
			
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
			    if(layer_to_use!=-1 && layer_to_use == hh.layer && hh.layer!=3 && hh.layer!=0 ){
				//println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
				def measured_beta = beta
				def ftof_sector = hh.sector
				def ftof_t = hh.time
				def ftof_p = hh.path+shiftp
				def ftof_paddle = hh.paddle
				def ftof_layer = hh.layer
				def measured_t  = ftof_t - vt

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

				//////////////////////////////////
				// delta t 
				//
				histos.computeIfAbsent("ftof_deltat_pr",histoBuilders.deltat).fill(p, delta_t_pr)
				histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pr)
				histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pr)
				histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pr)

				histos.computeIfAbsent("ftof_deltat_kp",histoBuilders.deltat).fill(p, delta_t_kp)
				histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_kp)
				histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_kp)
				histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_kp)

				histos.computeIfAbsent("ftof_deltat_pip",histoBuilders.deltat).fill(p, delta_t_pip)
				histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pip)
				histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pip)
				histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pip)


			    }
			}
		    }
		}
	    }
	    eventIndex++
	}
    }
}


def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("mon-deltatime-hadrons.hipo")
