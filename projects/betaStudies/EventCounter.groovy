package projects.betaStudies

class EventCounter{

    private static int counter = 0
    public static synchronized int getGlobalCounter(){
	return counter++
    }
}
    
