#!/opt/coatjava/6.3.1/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.clas.physics.Particle
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.jroot.ROOTFile

//groovy n java
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
//personal 
import examples.epkpkm
import mon.HadronMonV3
import event.EventConverter
import event.Event

//dmriser
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


Particle.metaClass.pindex = null
Particle.metaClass.sphi = null
Particle.metaClass.sector = null
LorentzVector.metaClass.vz = null
Vector3.metaClass.vz = null


def topology = 'epkpkm'
def field_setting = 'inbending'
// cut lvl meanings: 0 loose, 1 med, 2 tight
el_cut_strictness_lvl=["ecal_cut_lvl":1,
		       "nphe_cut_lvl":1,
		       "vz_cut_lvl":1,
		       "min_u_cut_lvl":1,
		       "min_v_cut_lvl":1,
		       "min_w_cut_lvl":1,
		       "max_u_cut_lvl":1,
		       "max_v_cut_lvl":1,
		       "max_w_cut_lvl":1,
		       "dcr1_cut_lvl":1,
		       "dcr2_cut_lvl":1,
		       "dcr3_cut_lvl":1,
		       "anti_pion_cut_lvl":1
]

def outname = args[0].split('/')[-1]
def out = new TDirectory()

def processors = [new HadronMonV3(field_config: field_setting, el_cut_strictness: el_cut_strictness_lvl)]
processors.each{it.setCutStrictness(); it.setCutParameters() }

processors.each{out.mkdir("/"+it.getClass().getSimpleName())}


GParsPool.withPool 16, {
    args.eachParallel{fname->
	def reader = new HipoDataSource()
	reader.open(fname)

	cc=0
	while(reader.hasEvent()) {
	    if (cc % 5000 == 0) {
                println("Processing " + cc)
            }
	    //if( cc == 100000 ) break

	    def data_event = reader.getNextEvent()	    
	    
	    def event = EventConverter.convert(data_event)	    
	    processors.each{it.processEvent(event,data_event)}
	    cc+=1
	}
	reader.close()
    }
}


//finally write output to a ROOT file
def outroot = new ROOTFile("monitor_hadronPIDv1_thesis_$outname"+".root")
processors.each{channel -> channel.histos.each{outroot.writeDataSet(it.value)}}
outroot.close()
