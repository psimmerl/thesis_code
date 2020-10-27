#!/opt/coatjava/6.3.1/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource

def fname = args[0]

def hiporeader = new HipoDataSource()
hiporeader.open(fname)
def max_event=hiporeader.getSize()
println(max_event);


