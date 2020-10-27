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

import examples.epkpkn
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

// 1 for SIMULATION, 0 for DATA
def DATATYPE = 0


def outname = args[0].split('/')[-1]
def out = new TDirectory()

def processors = [new epkpkn()]
def proc_names = ["test_pid"]
processors.each{out.mkdir("/"+it.getClass().getSimpleName())}

//processors.eachWithIndex{ proc, indx -> proc.createROOTFile(proc.getClass().getSimpleName()+"_"+proc_names[indx]+".root") }

def timer = new Timer()

timer.schedule({
    processors.each{
	out.cd("/"+it.getClass().getSimpleName())
	it.hists.each{out.addDataSet(it.value)}
    }
    out.writeFile("test_pid$outname")
} as TimerTask, 30000,60000)

GParsPool.withPool 16, {
    args.eachParallel{fname->
	def reader = new HipoDataSource()
	reader.open(fname)

	cc=0
	while(reader.hasEvent()){
	    def event = EventConverter.convert(reader.getNextEvent())
	    processors.each{it.processEventTest(event)}
	    cc+=1
	}

	reader.close()
    }
}

timer.cancel()
println('Finished Event Loop')

//finally write output to a ROOT file
def outroot = new ROOTFile("test_pid_$outname"+".root")
processors.each{channel -> channel.hists.each{outroot.writeDataSet(it.value)}}
outroot.close()
println('Finished Writing File')

processors.each{channel -> channel.closeOutTree() }

println('DONE')
