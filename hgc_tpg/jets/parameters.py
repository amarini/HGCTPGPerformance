from array import array
import numpy as np

class parameters :
    def __init__(self):
        ##
        self.calibration={"do":False,
                #"type":"inversion",
                "type":"none",
                "constants":[
                1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 
                ],
                "cut":[ ## raw, 0 below the cut
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                ],
                ## 
                "pu":"PU_cone_vs_eta_v0.root:prof_cone_C3D_pTgt0",
                "pt-eta":"",
                ## NO PU
                #"pt-eta":"4.4374 -0.948102 *TMath::Log(TMath::Max(x,20.)) + 0.0686934* TMath::Log(TMath::Max(x,20.)) * TMath::Log(TMath::Max(x,20.))", ## nopu
                ## pu200
                #"pt-eta":"0.577536*TMath::Erfc(-1.11941*x+4.30494)",
                
                ## pu140
                #"pt-eta": "(0.447 + 0.08173*TMath::Log(TMath::Max(x,20.)) -0.009685*TMath::Power(TMath::Log(TMath::Max(x,20.)),2)  ) * TMath::Erfc(-1.4256*TMath::Log(TMath::Max(x,20.)) +4.8851)",
                }
        
        ##
        self.matrix_calibration={"file":"matrix_calibration.root"}

        ## cluster parameters
        self.cluster = {"dR": float(0.2),"extraRadius":[0.1]}

        ## tree name
        self.output = {"tree":"jets","file":"output.root"}

        ## to run
        self.torun=['efficiency_pteta','fake_pteta','resolution_pt','turn_on','EnergyProfile','VbfSelection']
        #self.torun=['matrix_calibration'] ## set also calib to false
        ## 
        self.efficiency = {"ptbins":array('d',[20,30,40,50]),"etabins":array('d',[1.7,2.0,2.2,2.5,2.7]),"abseta":True,"dR":float(0.3)}
        self.fake = {"ptbins":array('d',[20,30,40,50]),"etabins":array('d',[1.7,2.0,2.2,2.5,2.7]),"abseta":True,"dR":float(0.3)}
        #self.resolution = {"bins":np.arange(-10.0,10.0, 0.01,dtype='d'),"dR":float(0.3),"ptbins":array('d',[30.,40,50,60,70,80,90,100,120,150,200,300,500,1000,2000,3000])}
        self.resolution = {"bins":np.arange(-2.0,10, .1,dtype='d'),"dR":float(0.3),"ptbins":array('d',[30.,40,50,60,70,80,90,100,120,150,200,300,500,1000,2000,3000])}
        #self.turnon = {"ptbins":array('d',range(30,1000)),"dR":float(0.3),"ptcut":[float(35.),float(80.),float(120),float(500.)]}
        self.turnon = {"ptbins":array('d',range(30,1000)),"dR":float(0.3),"ptcut":[float(50),float(100),float(125),float(150),float(200.)]}
        self.layer_deposits= {'ptcut':[30,40,50,80,100,200,500,1000],
                'dR':0.4
                }
