class ExclusiveCuts:

    def __init__(self):
        self.epkpkmxe_min = None
        self.epkpkmxe_max = None
        self.ekpkmX_min = None
        self.ekpkmX_max = None
        self.epkpX_min = None
        self.epkpX_max = None
        self.epkmX_min = None
        self.epkmX_max = None
        self.lambda1520_min = None
        self.lambda1520_max = None
        self.lambda1820_min = None
        self.lambda1820_max = None

    def setCuts(self,datatype, tag, sig):
        print(" IN THE NEW CLASS YEAH!")
        self.lambda1520_min = 1.49
        self.lambda1520_max = 1.587
        self.lambda1820_min = 1.75
        self.lambda1820_max = 1.89


        if datatype == 'sim' and 'inb' in tag or 'outb' in tag :
            print('init simulation cuts')
            self.epkpkmxe_min = 0.0222 - 0.0262*sig
            self.epkpkmxe_max = 0.0222 + 0.0262*sig
            
            self.epkpX_min = 0.234 - 0.0348*sig
            self.epkpX_max = 0.234 + 0.0348*sig
            
            self.epkmX_min = 0.228 - 0.0284*sig
            self.epkmX_max = 0.228 + 0.0284*sig
            
            self.ekpkmX_min = 0.891 - 0.0285*sig
            self.ekpkmX_max = 0.891 + 0.0285*sig
            
        if datatype == 'data' and 'inb' in tag:
            self.epkpkmxe_min = 0.08 - 0.0398*sig
            self.epkpkmxe_max = 0.08 + 0.0398*sig
            
            self.epkpX_min = 0.248 - 0.059*sig
            self.epkpX_max = 0.248 + 0.059*sig
            
            self.epkmX_min = 0.248 - 0.055*sig
            self.epkmX_max = 0.248 + 0.055*sig
            
            self.ekpkmX_min = 0.904 - 0.0719*sig
            self.ekpkmX_max = 0.904 + 0.0719*sig
            
        if datatype == 'data' and 'outb' in tag:
            self.epkpkmxe_min = 0.08 - 0.0398*sig
            self.epkpkmxe_max = 0.08 + 0.0398*sig
            
            self.epkpX_min = 0.248 - 0.059*sig
            self.epkpX_max = 0.248 + 0.059*sig
            
            self.epkmX_min = 0.248 - 0.055*sig
            self.epkmX_max = 0.248 + 0.055*sig
            
            self.ekpkmX_min = 0.904 - 0.0719*sig
            self.ekpkmX_max = 0.904 + 0.0719*sig
            
            
    def passEPKPKMX(self,epkpkmX):
        return ( (epkpkmX.E() < self.epkpkmxe_max and epkpkmX.E() > self.epkpkmxe_min ))

    def passEPKPX(self,epkpX):
        return epkpX.M2() < self.epkpX_max and epkpX.M2() > self.epkpX_min

    def passEPKMX(self,epkmX):
        return epkmX.M2() < self.epkmX_max and epkmX.M2() > self.epkmX_min

    def passEKPKMX(self,ekpkmX):
        return ekpkmX.M2() < self.ekpkmX_max and ekpkmX.M2() > self.ekpkmX_min

    def passLambda1520Cut(self, lmda):
        return lmda.M() < self.lambda1520_min or lmda.M() > self.lambda1520_max
        

    def passLambda1820Cut(self,lmda):
        return lmda.M() < self.lambda1820_min or lmda.M() > self.lambda1820_max

        
    
