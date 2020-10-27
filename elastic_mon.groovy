#!/opt/coatjava/6.3.1/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap

import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple

import mon.ElasticMon
import event.EventConverter
import event.Event


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


//USE GATED!!!
def convertFCtoLumi(acc_charge_ug, acc_charge_g){
    def tot_beam_charge = acc_charge_g
    def target_density = 0.0701; // g/cm^3
    def atomic_mass_hydrogen = 1.00794; // g/mol
    def avogado_const = 6.0221367E23; // Number/mol	
    def target_length = 5.0; //cm
    def cm_to_microbarn = 1E30;
    def el_charge = 1.602177E-19; // Coulomb
    def nC_to_C = 1e-9;
    def ns_to_s = 1e-9;

    println(" >> target density          " + target_density );
    println(" >> atomic mass of hydrogen " + atomic_mass_hydrogen );
    println(" >> avogado constant        " + avogado_const );
    println(" >> targe length            " + target_length );
    println(" >> cm to microbarn         " + cm_to_microbarn );
    println(" >> electron charge         " + el_charge );
    println(" >> nC to C                 " + nC_to_C );	

    def n_el = (tot_beam_charge * nC_to_C)/el_charge;

    def n_el_ug = acc_charge_ug *nC_to_C / el_charge;  
    def n_el_g = acc_charge_g *nC_to_C / el_charge;  

    def n_pr = ( ( target_length * target_density * avogado_const ) / atomic_mass_hydrogen ) ;       	
    def lum_factor = (n_el_g*n_pr)/cm_to_microbarn;       	

    println(" >> Total Gated Charge " +  acc_charge_g );
    println(" >> Total UNGated Charge " +  acc_charge_ug );
    println( " >> Number of electron " + n_el );
    println( " >> Number of electron GATED " + n_el_g );
    println( " >> Number of electron UNGATED " + n_el_ug );
    println( " >> Number of protons  " + n_pr );
    println( " >> Luminosity factor  " + lum_factor );

    return lum_factor
}

// 1 for SIMULATION, 0 for DATA
def DATATYPE = 1
def proc_method = -1
data_prefix=null
if( DATATYPE == 1 ){
    data_prefix='mc'
    proc_method = 3
}
else{
    data_prefix='data'
    proc_method = 2
}


def field_setting = 'outbending'
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

def processors = [new ElasticMon(field_config: field_setting, el_cut_strictness: el_cut_strictness_lvl)]
processors.each{ it.setCutStrictness(); it.setCutParameters() }
processors.each{out.mkdir("/"+it.getClass().getSimpleName())}

//def proc_names = [data_prefix+"_elasticXS_"+field_setting]
//processors.eachWithIndex{ proc, indx -> proc.createROOTFile(proc.getClass().getSimpleName()+"_"+proc_names[indx]+".root") }

//def timer = new Timer()

//timer.schedule({
//  processors.each{
//	out.cd("/"+it.getClass().getSimpleName())
//	//it.histos.each{out.addDataSet(it.value)}
//    }
//    //out.writeFile(data_prefix+"_elasticXS_"+"_"+field_setting+"$outname")
//} as TimerTask, 30000,60000)

///////////////////////////////////////
// variables for getting beam charge
def total_fcup_ug=0
def total_fcup_g=0
def total_beam_charge=0
def config_event=0
def current_event=0


GParsPool.withPool 16, {
    args.eachParallel{fname->
	def reader = new HipoDataSource()
	reader.open(fname)
	cc=0
	while(reader.hasEvent()){
	    if( cc % 50000 == 0 ){
		println(' Processing ' + cc )

	    }
	    if( cc == 1000000 ) break

	    /*
	    def hev = reader.getNextEvent()

	    def has_event = hev.hasBank("REC::Event")
	    def has_scaler = hev.hasBank("RUN::scaler")
	    def has_config = hev.hasBank("RUN::config")
	    if( has_config ){
		config_event = hev.getBank("RUN::config").getInt("event",0)
	    }
	    if( has_event ){
		def eventBank = hev.getBank("REC::Event")
		def beamcharge = eventBank.getFloat("beamCharge",0)
		if( beamcharge > total_beam_charge ){
  		    total_beam_charge = beamcharge 
		}
	    }
	    if( has_scaler ){
		def scalerBank = hev.getBank("RUN::scaler")
		def fcup_ug = scalerBank.getFloat("fcup",0)
		def fcup_g = scalerBank.getFloat("fcupgated",0)
		if( config_event > current_event ){
		    //println(" current config event " + current_event + " total fcup ug " + total_fcup_ug )
		    //println(" current config event " + current_event + " total fcup g " + total_fcup_g )
		    //println(" current config event " + current_event + " total beam charge is " + total_beam_charge)
		    total_fcup_ug = fcup_ug
		    total_fcup_g = fcup_g
		    current_event = config_event
		    
		}
	    }
	    */


	    def event = EventConverter.convert(reader.getNextEvent())
	    processors.each{it.processEvent(event)}
	    cc+=1
	}
	reader.close()
    }
}

println(">> TOTAL BEAM CHARGE " + total_beam_charge )
//convertFCtoLumi(total_fcup_ug, total_fcup_g)


//timer.cancel()
println('Finished Event Loop')

//finally write output to a ROOT file

def outroot = new ROOTFile(data_prefix+"_elasticXS_$outname"+"_"+field_setting+".root")
processors.each{channel -> channel.histos.each{outroot.writeDataSet(it.value)}}
outroot.close()
println('Finished Writing Histogram File')

//println('Writing TTree to ROOT File')
//processors.each{channel -> channel.closeOutTree() }

println('DONE')
