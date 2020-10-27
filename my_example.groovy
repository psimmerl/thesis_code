import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import pid.electron.Electron

def hemm2th = [:].withDefault{new H2F("hemm2th_$it", "MM2 vs theta", 250,5,30, 250,0,5)}
def hemm1fi = [:].withDefault{new H2F("hemm2fi_$it", "MM2 vs phi", 250,140,190, 250,0,5)}

def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

def rotateDCHitPosition(hit,sec){
//    println(" in hit rotator xy $hit , sector $sec")
    def ang = Math.toRadians(sec*60.0)
    def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
    def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
  //  println([x1_rot,y1_rot])
    return [x1_rot,y1_rot]
}


def electronPIDStrategies = [   
    { partb , iele -> partb.getInt('pid',iele) == 11 },
    { partb , iele -> partb.getInt('status',iele)<0 }
]

def electronPIDStrategies3 = [   
    { banks , iele -> partb.getInt('pid',iele) == 11 },
    { banks , iele -> partb.getInt('status',iele)<0 }

]

def electronPIDStrategies2 = [   
    { partb -> partb.getInt('pid').collect{ it==11 } },
    { partb -> partb.getInt('status').collect{ it<0 } }
]


//call electron cut constructor
electron = new Electron();

//define cut strategy
//def myElectronCutStrategies = [
//  {banks, index -> electron.passElectronEBPIDCut(banks,index)}
//]
/*def myElectronCutStrategies = [
    { banks, index -> electron.passElectronChargeCut(banks,index) },
    { banks, index -> electron.passElectronEBPIDCut(banks,index) },
    { banks, index -> electron.passElectronNpheCut(banks,index) },
    { banks, index -> electron.passElectronVertexCut(banks,index) },
    { banks, index -> electron.passElectronPCALFiducialCut(banks,index) },
    { banks, index -> electron.passElectronEIEOCut(banks,index) },
    { banks, index -> electron.passElectronDCR1(banks,index) },
    //{ banks, index -> electron.passElectronDCR2(banks,index) },
    //{ banks, index -> electron.passElectronDCR3(banks,index) }
]
*/

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
//out.mkdir("/dcr1xy")

def histos = new ConcurrentHashMap()
//def dcxyr1 = {new H2F("$it","$it", 500, -400, 400, 500 -400, 400 ) }
histoBuilders = [

    dcxy : { title -> new H2F("$title","$title", 500, -400.0, 400.0, 500, -400.0, 400.0 ) }

]

for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    def cc = 0
    
    while(reader.hasEvent() && cc < 20000) {
	
	def event = reader.getNextEvent()

	def banks  = reqBanks.collect{name, bank -> [name, event.getBank(bank)] }.collectEntries()
	def banks_pres = reqBanks.every{name, bank -> event.hasBank(bank) }
	//println(" are all banks present ")
	//println(banks_pres)

	
	if (banks_pres) { //event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Cherenkov") ) {
	        cc+=1
	        // -method 1 we can directly loop over the banks here the same old way
	        // -method 2 loop over objects in banks 
	        def partb = event.getBank("REC::Particle")
	        def calb = event.getBank("REC::Calorimeter")
	        def chb = event.getBank("REC::Cherenkov")

	        //get rows in bank pass to findAll, a closure, iterate over rec bank status and print it 
	      //  println( ' get stat ' )
	        //def my_el_stat = (0..<partb.rows()).findAll{println(partb.getShort('status',it)) }
	       // println( ' print findall ' )
	        //println(my_el_stat)

	        //findAll returns list of indices satisfying the closure
	        //send to findResults
	        //use each whn iterating over list
	        //reserve findAll and findResults  find all elements matching critieria
	        //def my_el_pid  =  (0..<partb.rows()).each{partb.getInt('pid',it) == 11 && partb.getInt('status',it)<0}.each{ ii -> println( ['px','py','pz'].collect{partb.getFloat(it,ii) } ) }
	        //println( my_el_pid )

	        //println(" trying closure " )
	        //def test_el = (0..<partb.rows()).findAll{ electronPIDStrategies.each{ el_test -> el_test(partb,it) }.findResults{ii -> println(ii) }  
	        //loop over rec particle bank -> use findAll to return list of indices to loop over with 'it' -> pass the rec.particle index to method electronPIDStrategies -> find all instances with true and return the index
	        //def test_el = (0..<partb.rows()).findAll{ electronPIDStrategies.findResults{ el_test -> el_test(partb,it) }.findAll{println(it) } }//Result{ii= -> if(ii== true) }  }
	        //def test_el2 = (0..<partb.rows()).findAll{ println(electronPIDStrategies.findResults{ el_test -> el_test(partb,it) }) }//findResult{ii= -> if(ii== true) }  }
	        //println(test_el2)

	        //important note: when making dict do not put ':', or it will literally take the object as the key
	        def test_el_map = (0..<partb.rows()).collect{ ii -> [ii , electronPIDStrategies.collect{ el_test -> el_test(partb,ii) } ] }.collectEntries()
	        //println(test_el_map)
	        //println(test_el_map.getClass())
	        
	        
	        def test_el_final = (0..<partb.rows()).findIndexValues{ electronPIDStrategies.every{ el_test -> el_test(partb,it) }  }
	        
	        //def test_el_final = (0..<partb.rows()).findIndexValues{ electronPIDStrategies.every{ el_test -> el_test(partb,it) }  }
	        //println(' test final list of indices ')
	        //println( test_el_final)
	        
	        //println(" showing traj bank " )
	        //banks.traj.show()
	        //println(" showing trck bank " )
	        //banks.trck.show()
	    // 3rd method - no loop  
	        def test_no_loop = electronPIDStrategies2.collect{ el_strat2 -> el_strat2(partb) }.transpose()

	        //println('testing any method ')
	        //println([1, 2, 3].any{ (it > 2) }   )


	        def my_el_cuts = (0..<partb.rows()).collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(banks,ii) } ] }.collectEntries()
	        //println(' my electron cuts ' )
	        //println(my_el_cuts)
	        
	        //test how to get the sector from cal bank for correct particle
	        //partb.show()
	        //calb.show()
	        //def sect  = (0..<partb.rows()).each{ ii -> (calb.getShort('pindex')*.toInteger()).findResults{ println(it) == ii  } }//calb.getShort('pindex',ii) ) 


	        // now lets try to get the nphe for events 
	        //map nphe[pindex] = nphe
	        def nphe = [chb.getShort('pindex'), chb.getFloat('nphe')].transpose().collectEntries{pind,nph-> [(pind.toInteger()):nph] }
	        //println( ' nphe ')
	        //println( nphe )

	        //can we get it without invoking map
	        //put the index that you want the number of photoelect for.
	        //def phe = (0..<banks.cc.rows()).find{ pindex -> banks.cc.getFloat("nphe", pindex) > 2 && banks.cc.getInt("pindex", pindex) == index }
	    //    def nphe = (0..<partb.rows()).find{cher.getByte('detector',it)==7}
	    //.collectEntries{(cherb.guODetShort('pindex',it)):cherb.getFloat('energy',it) }
	    //    nphe[0]>5
	        
	        
	        //def nphe  = (0..<partb.rows()).findResult{ chb.getInt('nphe') 
							        
	        //println(' sector ' )
	        //println( sect  )
	        //println( calb.getShort('pindex')*.toInteger().getClass() )
	        
	        //nice way to check range conditions
	        //banks.part.getFloat('vz',0).with{println(it)}// < max_vz && it > min_vz}
	        //println((0..<banks.part.rows()).findAll{ banks.part.getFloat('vz',0).with{it < 10 && it > -10} })
	        
	        //println(' testing no loop method ')
	       // println(test_no_loop)


	        //testing any method 
	        
	        //println( ' test object type ' )
	        //println( my_el_pid.getClass() )

	        //println((0..<banks.cc.rows()).any{banks.cc.getInt('pindex',it) == 0 && banks.cc.getFloat('nphe',it) > 2})

	    //def temp = ((0..<banks.ec.rows()).any{ ( (banks.ec.getByte('detector',it).find{ ii-> ii == 7 }.find{ println(it) } ) ) } )//banks.ec.getInt('pindex',it)) } ) ) } )// == DetectorType.PCAL.getDetectorId() && bank.ec.getInt('pindex',it) == index && bank.ec.getFloat('lu',it).with{it < max_u && it > min_u } && bank.ec.getFloat('lv',it).with{it < max_v && it > min_v } && bank.ec.getFloat('lw',it).with{it < max_w && it > min_w } }
	        
	        //println(' testing dc cuts ')
	        //def hit_pos = (0..<banks.traj.rows()).find{(banks.traj.getInt('pindex',it) == 0 && 
	    //     banks.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() &&
	        //banks.traj.getByte('layer',it) == 12 )}.collect{ ['cx','cy','cz'].collect{ii-> banks.traj.getFloat(ii,it)} }
	        //def dcsect  = (0..<banks.trck.rows()).find{ banks.trck.getInt('pindex',it) == 0 && 
	    //       banks.trck.getByte('detector',it) == DetectorType.DC.getDetectorId() }.collect{banks.trck.getByte('sector',it)//} 
	    //    println(dcsect)
	        //banks.traj.getFloat('cx',it) } //['cx','cy','cz'].collect{ banks.traj.getFloat(it,ii) }}
	        //['cx','cy','cz'].collect{ (banks.traj.getFloat(it,ii)) } }
	        //println(test)
	        //quick test to check hit position  - check only pindex == 0 - mostly electron
	        def hit_pos = (0..<banks.traj.rows()).findResult{(banks.traj.getShort('pindex',it) == 0 &&
								       banks.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() &&
								       banks.traj.getByte('layer',it) == 12 )
								     ? 'xy'.collect{ axis -> banks.traj.getFloat(axis,it) } : null }
	        def sectordc = (0..<banks.trck.rows()).findResult{(banks.trck.getShort('pindex',it) == 0 && 
								         banks.trck.getByte('detector',it) == DetectorType.DC.getDetectorId()) ? banks.trck.getByte('sector',it) : null }

	        if( hit_pos != null && sectordc != null){
		hit_rot = rotateDCHitPosition(hit_pos,sectordc-1)
		histos.computeIfAbsent("dcr1xy",histoBuilders.dcxy).fill(hit_rot.get(0),hit_rot.get(1))
		    }
	    }
    }
    
    out.mkdir("/test/")
    out.cd("/test/")
    histos.values().each{out.addDataSet(it)}
    out.writeFile("dcr1xy.hipo")

    reader.close();
}
