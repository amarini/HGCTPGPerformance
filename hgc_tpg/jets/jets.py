## producer of jet out of clusters
## ntuples version
## Andrea Carlo Marini -- Thu Jun 29 15:30:22 CEST 2017

import ROOT
import hgc_tpg.utilities.functions as f
import numpy as np
import math
from glob import glob
import sys
from datetime import datetime

## using pyjet
from pyjet import cluster, DTYPE_EP, DTYPE_PTEPM


class BaseRun:
    def __init__(self,outer):
        self.outer=outer

    def run(self):
        return self

    def finalize(self):
        return self

class VbfSelection(BaseRun):
    def __init__(self,outer):
        self.outer=outer
        self.histos={}
        self.histos["CutFlow"] = ROOT . TH1D("CutFlow","AllGen,match02,match01,matchboth",10,0,10)

        self.histos["Mjj_gen"] = ROOT . TH1D("Mjj_gen","Mjj gen level",250,0,1000)
        self.histos["Mjj_endcaps_gen"] = ROOT . TH1D("Mjj_endcaps_gen","Mjj gen level. ",250,0,1000)
        self.histos["Mjj_reco"] = ROOT . TH1D("Mjj_reco","Mjj reco level",250,0,1000)

        self.histos["axis1"] = ROOT . TH1D("axis1","-Log(axis1)",100,0,10.)
        self.histos["axis2"] = ROOT . TH1D("axis2","-Log(axis2)",100,0,10.)
        self.histos["mult"] = ROOT . TH1D("mult","mult",100,0,100)
        self.histos["mult1GeV"] = ROOT . TH1D("mult1GeV","mult",100,0,100)
        self.histos["ptD"] = ROOT . TH1D("ptD","ptD",100,0,1.)
        self.histos["dR"] = ROOT . TH1D("dR","-Log(dR)",100,0,10.)

    def run( self):
        ## find VBF genjets
        j1_gen,j2_gen,dijet_gen=None,None,None
        mjj=0
        dijet_gen_endcaps=None

        for igen, gj1 in enumerate(self.outer.genjets_):
            if gj1.Pt()<30: continue
            isE1=abs(gj1.Eta()) > 3.0 or abs(gj1.Eta()) < 1.5 
            for jgen, gj2 in enumerate(self.outer.genjets_):
                if igen <= jgen: continue
                if gj2.Pt()<30: continue
                isE2=abs(gj2.Eta()) > 3.0 or abs(gj2.Eta()) < 1.5 
                dijet_gen=ROOT.TLorentzVector(0,0,0,0,)
                dijet_gen += gj1
                dijet_gen += gj2
                if dijet_gen.M() > mjj: 
                    j1_gen=gj1
                    j2_gen=gj2
                    mjj=dijet_gen.M()
                if isE1 and isE2 and ( dijet_gen_endcaps== None or dijet_gen.M() > dijet_gen_endcaps.M()):
                    dijet_gen_endcaps=dijet_gen
        ## 
        self.histos["Mjj_gen"].Fill(mjj)
        self.histos["Mjj_endcaps_gen"].Fill(dijet_gen_endcaps.M() if dijet_gen_endcaps != None else 0)
        ## reco loop
        j1_trg,j2_trg,dijet_trg=None,None,None
        j1_idx,j2_idx = -1,-1
        mjj=0
        for ijet, j1 in enumerate(self.outer.jets_):
            if j1.Pt()<30: continue
            if abs(j1.Eta()) > 3.0 or abs(j1.Eta()) < 1.5 : continue
            for jjet, j2 in enumerate(self.outer.genjets_):
                if ijet <= jjet: continue
                if j2.Pt()<30: continue
                if abs(j2.Eta()) > 3.0 or abs(j2.Eta()) < 1.5 : continue
                dijet_trg=ROOT.TLorentzVector(0,0,0,0,)
                dijet_trg += j1
                dijet_trg += j2
                if dijet_trg.M() > mjj: 
                    j1_trg=j1
                    j2_trg=j2
                    mjj=dijet_trg.M()
                    j1_idx=ijet
                    j2_idx=jjet
        self.histos["Mjj_reco"].Fill(mjj)

        if j1_trg != None and j2_trg!=None:
            for idx in [j1_idx, j2_idx]:
                ## fill properties
                self.histos["axis1"]    . Fill(-math.log(self.outer.jetprop_[idx]["axis1"]) if self.outer.jetprop_[idx]["axis1"] >0 else 100)
                self.histos["axis2"]    . Fill(-math.log(self.outer.jetprop_[idx]["axis2"]) if self.outer.jetprop_[idx]["axis2"] >0 else 100)
                self.histos["mult"]     . Fill(self.outer.jetprop_[idx]["mult"])
                self.histos["mult1GeV"] . Fill(self.outer.jetprop_[idx]["mult1GeV"])
                self.histos["ptD"]      . Fill(self.outer.jetprop_[idx]["ptD"])
                self.histos["dR"]       . Fill(self.outer.jetprop_[idx]["dR"])
        ##
        return self

    def finalize(self):
        for h_str in self.histos  :
            h=self.histos[h_str]
            h.Write()
        return self

class EnergyProfile(BaseRun):
    def __init__(self,outer):
        self.outer=outer
        self.histos={}
        self.dR=0.4

        self.histos["EnergyProfile"] = ROOT . TH2D("EnergyProfile","EnergyProfile;E0.1/E0.2;E0.2/ETrue",1000,0,10,1000,0,10)
        self.histos["CutFlow"] = ROOT . TH1D("CutFlow","AllGen,match02,match01,matchboth",10,0,10)

    def run(self):
        ## match true -> E02 -> E01
        for igen, gj in enumerate(self.outer.genjets_):
            if gj.Pt()<30: continue
            if abs(gj.Eta()) > 3.0 or abs(gj.Eta()) < 1.5 : continue
            ET  = gj.Energy()
            E02 = 0.
            E01 = 0.
            match1=False
            match2=False
            #print "Considering Jet PT=",gj.Pt(),",",gj.Eta(),gj.Phi()
            dRmin=self.dR
            for itr, jet in enumerate(self.outer.jetsR_[0.2]):
                #print " RECO 02 Jet PT=",jet.Pt(),":",jet.Eta(),jet.Phi(),"DR=", gj.DeltaR(jet)
                if jet.Pt() < 20: continue
                if gj.DeltaR(jet) < dRmin: 
                    dRmin=gj.DeltaR(jet)
                    E02 = jet.Energy()
                    match1 = True
                    #break
            dRmin=self.dR
            for itr, jet in enumerate(self.outer.jetsR_[0.1]):
                if jet.Pt() < 20: continue
                if gj.DeltaR(jet) < dRmin: 
                    dRmin=gj.DeltaR(jet)
                    E01 = jet.Energy()
                    match2 = True
                    #break
        
            self.histos["CutFlow"].Fill(0)
            if match1: self.histos["CutFlow"].Fill(1)
            if match2: self.histos["CutFlow"].Fill(2)
            if match1 and match2:
                self.histos["CutFlow"].Fill(3)
                #self.histos["EnergyProfile"].Fill(ET/E02,E01/E02)
                self.histos["EnergyProfile"].Fill(E01/E02,E02/ET)


        return self

    def finalize(self):
        for h_str in self.histos  :
            h=self.histos[h_str]
            h.Write()
        return self

class efficiency_pteta(BaseRun):
    def __init__(self,outer):
        self.outer=outer

        ptbins=outer.cfg.efficiency['ptbins']
        etabins=outer.cfg.efficiency['etabins']
        self.eff_num = ROOT.TH2D("eff_num","eff. num",len(ptbins)-1,ptbins,len(etabins)-1,etabins)
        self.eff_den = ROOT.TH2D("eff_den","eff. den",len(ptbins)-1,ptbins,len(etabins)-1,etabins)
        self.useAbsEta = outer.cfg.efficiency["abseta"]
        self.dR= outer.cfg.efficiency["dR"]

    def run(self):
        for igen, gj in enumerate(self.outer.genjets_):
            match=False
            for itr, jet in enumerate(self.outer.jets_):
                if gj.DeltaR(jet) < self.dR: 
                    match = True
                    break
            eta = gj.Eta()
            pt = gj.Pt()
            if self.useAbsEta: eta = abs(eta)
            self.eff_den.Fill(pt,eta)
            if match:  self.eff_num.Fill(pt,eta)
        return self

    def finalize(self):        
        self.eff_h = self.eff_num.Clone("eff")
        self.eff_h.Divide(self.eff_den)
        for h in [ self.eff_h,self.eff_num,self.eff_den] :
            h.Write()
        return self

class fake_pteta(BaseRun):
    def __init__(self,outer):
        self.outer=outer

        ptbins=outer.cfg.fake['ptbins']
        etabins=outer.cfg.fake['etabins']
        self.fake_num = ROOT.TH2D("fake_num","fake. num",len(ptbins)-1,ptbins,len(etabins)-1,etabins)
        self.fake_den = ROOT.TH2D("fake_den","fake. den",len(ptbins)-1,ptbins,len(etabins)-1,etabins)
        self.useAbsEta = outer.cfg.fake["abseta"]
        self.dR= outer.cfg.fake["dR"]

    def run(self):
        for itr, jet in enumerate(self.outer.jets_):
            match=False
            for igen, gj in enumerate(self.outer.genjets_):
                if gj.DeltaR(jet) < self.dR: 
                    match = True
                    break
            eta = jet.Eta()
            pt = jet.Pt()
            if self.useAbsEta: eta = abs(eta)
            self.fake_den.Fill(pt,eta)
            if match:  self.fake_num.Fill(pt,eta)
        return self

    def finalize(self):        
        self.fake_h = self.fake_num.Clone("fake")
        self.fake_h.Divide(self.fake_den)
        for h in [ self.fake_h,self.fake_num,self.fake_den] :
            h.Write()
        return self

class resolution_pt(BaseRun):
    def __init__(self,outer):
        self.outer=outer

        bins=outer.cfg.resolution['bins']
        self.ptbins=outer.cfg.resolution['ptbins']
        self.res_h = {}
        for idx in range(0,len(self.ptbins)-1):
            self.res_h['pt%.0f_%.0f'%(self.ptbins[idx],self.ptbins[idx+1])]=ROOT.TH1D("res_pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1]),"res. pt %.0f %.0f"%(self.ptbins[idx],self.ptbins[idx+1]),len(bins)-1,bins)
            self.res_h['gen_pt%.0f_%.0f'%(self.ptbins[idx],self.ptbins[idx+1])]=ROOT.TH1D("genres_pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1]),"res. pt %.0f %.0f"%(self.ptbins[idx],self.ptbins[idx+1]),len(bins)-1,bins)
            self.res_h['recodist_pt%.0f_%.0f'%(self.ptbins[idx],self.ptbins[idx+1])]=ROOT.TH1D("recodist_pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1]),"res. pt %.0f %.0f"%(self.ptbins[idx],self.ptbins[idx+1]),100,self.ptbins[idx],self.ptbins[idx+1])
        self.dR= outer.cfg.resolution["dR"]

    def run(self):
        for itr, jet in enumerate(self.outer.jets_):
            if jet.Pt()<20 : continue
            #if abs(jet.Eta()) > 3.0 or abs(jet.Eta()) < 1.5 : continue
            if abs(jet.Eta()) > 2.7 or abs(jet.Eta()) < 1.7 : continue  ## RECO stricter 
            match=False
            pt_gen=0
            for igen, gj in enumerate(self.outer.genjets_):
                if abs(gj.Eta()) > 3.0 or abs(gj.Eta()) < 1.5 : continue
                if gj.DeltaR(jet) < self.dR: 
                    match = True
                    pt_gen=gj.Pt()
                    break
            if match and pt_gen>5:  
                res = (jet.Pt() - pt_gen)/pt_gen
                genres = pt_gen/jet.Pt()
                for idx in range(0,len(self.ptbins)-1):
                    if pt_gen >= self.ptbins[idx] and pt_gen< self.ptbins[idx+1]:
                        self.res_h["pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1])].Fill(res)
                    ## to derive correctios: gen/reco and reco mean
                    if jet.Pt() >= self.ptbins[idx] and jet.Pt()< self.ptbins[idx+1]:
                        self.res_h["gen_pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1])].Fill(genres)
                        self.res_h["recodist_pt%.0f_%.0f"%(self.ptbins[idx],self.ptbins[idx+1])].Fill(jet.Pt())
        return self

    def finalize(self):        
        for h_str in self.res_h  :
            h=self.res_h[h_str]
            h.Write()
        return self

class turn_on(BaseRun):
    def __init__(self,outer):
        self.outer=outer

        self.ptcut=outer.cfg.turnon['ptcut']
        ptbins=outer.cfg.turnon['ptbins']
        self.dR= outer.cfg.turnon["dR"]
        
        self.turnon_num = {}
        self.turnon_den = {}
        for pt in self.ptcut:
            self.turnon_num["%.0f"%pt] = ROOT.TH1D("turnon_pt%.0f_num"%pt,"turnon.",len(ptbins)-1,ptbins)
            self.turnon_den["%.0f"%pt] = ROOT.TH1D("turnon_pt%.0f_den"%pt,"turnon.",len(ptbins)-1,ptbins)

    def run(self):
        for igen, gj in enumerate(self.outer.genjets_):
            match=False
            pt_reco=0.
            if abs(gj.Eta()) > 3.0 or abs(gj.Eta()) < 1.5 : continue
            for itr, jet in enumerate(self.outer.jets_):
                if abs(jet.Eta()) > 3.0 or abs(jet.Eta()) < 1.5 : continue
                if gj.DeltaR(jet) < self.dR: 
                    match = True
                    pt_reco=jet.Pt()
                    break
            if match:
                pt_gen=gj.Pt()
                for pt in self.ptcut:
                    self.turnon_den["%.0f"%pt].Fill(pt_gen)
                    if pt_reco> pt: self.turnon_num["%.0f"%pt].Fill(pt_gen)
        return self

    def finalize(self):        
        histos=[]
        for pt in self.ptcut:
            turnon_h = self.turnon_num["%.0f"%pt].Clone("turnon_pt%.0f"%pt)
            turnon_h.Divide(self.turnon_den["%.0f"%pt])
            histos.extend([turnon_h,self.turnon_num["%.0f"%pt],self.turnon_den["%.0f"%pt]])
        for h in histos :
            h.Write()
        return self

class layer_deposits(BaseRun):
    def __init__(self,outer):
        self.outer=outer
        self.ptcuts=outer.cfg.layer_deposits['ptcut']

        self.dR=outer.cfg.layer_deposits['dR']

        self.histos={}
        for ipt in range(0,len(self.ptcuts) -1 ):
         for layer in range(1,29):
            self.histos["raw_layer%d_pt%.0f_%.0f"%(layer,self.ptcuts[ipt],self.ptcuts[ipt+1])] = ROOT.TH1D("raw_layer%d_pt%.0f_%.0f"%(layer,self.ptcuts[ipt],self.ptcuts[ipt+1]),"",1000,0,1000)
            self.histos["cl_layer%d_pt%.0f_%.0f"%(layer,self.ptcuts[ipt],self.ptcuts[ipt+1])] = ROOT.TH1D("cl_layer%d_pt%.0f_%.0f"%(layer,self.ptcuts[ipt],self.ptcuts[ipt+1]),"",1000,0,1000)

    def run(self):

        for igen, gj in enumerate(self.outer.genjets_):
            match=False
            if abs(gj.Eta()) > 3.0 or abs(gj.Eta()) < 1.5 : continue
            pt_gen=gj.Pt()

            pt0 =-1
            pt1 =-1
            for ipt in range( 0,len(self.ptcuts) -1 ):
                if pt_gen >= self.ptcuts[ipt] and pt_gen<self.ptcuts[ipt+1]:
                    pt0 = self.ptcuts[ipt]
                    pt1 = self.ptcuts[ipt+1]
            if pt0 < 0: continue ## not in the pt list

            isReco=False
            for itr, jet in enumerate(self.outer.jets_):
                if jet.Pt()<20 : continue
                if abs(jet.Eta()) > 3.0 or abs(jet.Eta()) < 1.5 : continue
                if gj.DeltaR(jet) < 0.4: isReco=True
            if not isReco: continue

            for icl, di in enumerate(self.outer.cl_raw_):
                ieta =  self.outer.cl_eta_[icl]
                iphi = self.outer.cl_phi_[icl] 
                if f.deltaR(gj.Eta(),ieta, gj.Phi(),iphi ) >= self.dR: continue
                ilayer = self.outer.cl_layer_[icl]
                ipt =  self.outer.cl_pt_[icl]
                iraw = self.outer.cl_raw_[icl]
                self.histos["raw_layer%d_pt%.0f_%.0f"%(ilayer,pt0,pt1)] . Fill (iraw) 
                self.histos["cl_layer%d_pt%.0f_%.0f"%(ilayer,pt0,pt1)] . Fill (ipt) 

        return self

    def finalize(self):        
        for name in self.histos:
            self.histos[name].Write()
        return self

class matrix_calibration(BaseRun):
    def __init__(self,outer):
        self.outer=outer
        self.dR=0.4

        self.histos={}
        self.histos["N"] = ROOT.TH1D("norm","norm",1,0,1)
        self.histos["M"] = ROOT.TH2D("matrix","matrix",50,0,50,50,0,50)
        self.histos["v"] = ROOT.TH1D("vector","vector",50,0,50)
        

    def run(self):
        for igen, gj in enumerate(self.outer.genjets_):
            match=False
            if abs(gj.Eta()) > 3.0 or abs(gj.Eta()) < 1.5 : continue
            pt_gen=gj.Pt()
            if pt_gen < 40 : continue

            isReco=False
            for itr, jet in enumerate(self.outer.jets_):
                if jet.Pt()<20 : continue
                #if jet.Pt()<200 : continue ##FIXME
                if abs(jet.Eta()) > 3.0 or abs(jet.Eta()) < 1.5 : continue
                if gj.DeltaR(jet) < 0.4: isReco=True
            if not isReco: continue

            for icl, di in enumerate(self.outer.cl_raw_):
                ieta =  self.outer.cl_eta_[icl]
                iphi = self.outer.cl_phi_[icl] 
                ilayer = self.outer.cl_layer_[icl]
                if f.deltaR(gj.Eta(),ieta, gj.Phi(),iphi ) >= self.dR: continue
                for jcl, dj in enumerate(self.outer.cl_raw_):
                    jeta =  self.outer.cl_eta_[jcl]
                    jphi = self.outer.cl_phi_[jcl] 
                    jlayer = self.outer.cl_layer_[jcl]
                    if f.deltaR(gj.Eta(),jeta, gj.Phi(),jphi ) >= self.dR: continue
                    match = True
                    self.histos["M"].Fill(ilayer,jlayer,di*dj)
                self.histos["v"].Fill(ilayer,di*pt_gen)

            if match:
                self.histos["N"].Fill(1)

        return self

    def finalize(self):        
        for name in self.histos:
            self.histos[name].Write()
        return self


class jet_clustering:
    ''' Jet Clustering '''
    def __init__(self, input_file, conf):
        self.chain = ROOT.TChain("hgcalTriggerNtuplizer/HGCalTriggerNtuple")
        for template_name in input_file.split(','):
            if '*' in template_name:
                for f in glob(template_name):
                    self.chain.Add(f)
            else: ## regular files, xrootd, ...
                self.chain.Add(template_name)

        self.output = ROOT.TFile(conf.output['file'],"RECREATE")
        self.tree_out = ROOT.TTree(conf.output ['tree'],conf.output['tree'])
        self.cfg = conf

        self.doCalibration=conf.calibration["do"]

        self.extraRadius=conf.cluster["extraRadius"]
        self.clusterR = conf.cluster["dR"]

        if self.doCalibration:
            print "============= USING CALIBRATION SEQUENCE ===========:"
            self.printCalib=0
            self.calibration=conf.calibration['type']

            if conf.calibration['type']=='inversion':
                self.const=conf.calibration["constants"][:]
                self.cut = conf.calibration["cut"][:]
                print "<*> Using inversion method for calibration:"
                print "    coeff:" + ','.join(["%.3f"%x for x in self.const])
                print "     cuts:" + ','.join(["%.3f"%x for x in self.cut])

            self.pu = conf.calibration["pu"]
            if self.pu !="":
                fname,pname = self.pu.split(":")[0], self.pu.split(":")[1]
                self.fPU = ROOT.TFile.Open(fname)
                self.pPU = self.fPU.Get(pname)
                print "<*> Using pu corrections:",fname,pname

            self.pteta = conf.calibration["pt-eta"]
            if self.pteta !="":
                self.calibf = ROOT.TF1("myfunc",self.pteta,0,10000)
                print "<*> Using pt-eta method for calibration:"
                print "   pt-eta:" + self.pteta
            else:
                self.calibf = None

        ## Collections of TLorentzVector
        self.jets_ = []
        self.genjets_ = []
        self.time_=None
        #jets with different radius: "0.1" -> []
        self.jetsR_ = {}

    def clear(self):
        self.jets_=[]
        self.jetprop_=[] # jetprop_[0]["xxx"] #radius, axis1, axis2, ptD, mult
        self.genjets_=[]
        self.jetsR_ = {}

    def do_clusters(self,ptV,etaV,phiV,energyV ):
        ''' run the clustering on the input std vectors'''
        n=len(ptV)
        if len(etaV) != n : raise ValueError('input vectors have different length')
        if len(phiV) != n : raise ValueError('input vectors have different length')
        if len(energyV) != n : raise ValueError('input vectors have different length')
        
        input_particles = []
        for i in range(0,n):
            p = ROOT.TLorentzVector()
            p.SetPtEtaPhiE(ptV[i],etaV[i],phiV[i],energyV[i])
            input_particles.append( (p.Pt(),p.Eta(),p.Phi(),p.M())  )
       
        ip = np.array( input_particles, dtype=DTYPE_PTEPM)
        sequence=cluster( ip, algo="antikt",ep=False,R=self.clusterR)
        jets = sequence.inclusive_jets()

        for i, jet in enumerate(jets):
            p = ROOT.TLorentzVector()
            p.SetPtEtaPhiM(  jet.pt, jet.eta, jet.phi, jet.mass)
            self.jets_ .append(p) 

            ## compute properties
            sum_weight = 0.
            sum_pt = 0.
            sum_deta = 0.
            sum_dphi = 0.
            sum_deta2 = 0.
            sum_detadphi = 0.
            sum_dphi2 = 0.
            sum_dr2 = 0.
            sum_dr  = 0.
            mult=0
            mult1GeV=0
            for c in jet.constituents():
                mult += 1
                if c.pt > 1: mult1GeV+=1
                weight= c.pt * c.pt
                deta=c.eta -jet.eta
                dphi=f.deltaPhi(c.phi,jet.phi)
                sum_weight += weight;
                sum_pt += c.pt;
                sum_deta += deta*weight;
                sum_dphi += dphi*weight;
                sum_deta2 += deta*deta*weight;
                sum_detadphi += deta*dphi*weight;
                sum_dphi2 += dphi*dphi*weight;
                sum_dr2 += (deta*deta + dphi*dphi)*weight
                sum_dr += math.sqrt(deta*deta + dphi*dphi)*weight
            #Calculate axis2 and ptD
            a,b,c = 0.,0.,0.
            ave_deta, ave_dphi,ave_deta2,ave_dphi2  = 0., 0.,  0.,  0.
            ave_dr2 = 0.
            if sum_weight > 0:
                ave_deta = sum_deta/sum_weight;
                ave_dphi = sum_dphi/sum_weight;
                ave_deta2 = sum_deta2/sum_weight;
                ave_dphi2 = sum_dphi2/sum_weight;
                ave_dr2 = sum_dr2/sum_weight
                ave_dr = sum_dr/sum_weight
                a = ave_deta2 - ave_deta*ave_deta;                          
                b = ave_dphi2 - ave_dphi*ave_dphi;                          
                c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
            delta = math.sqrt(abs((a-b)*(a-b)+4*c*c));
            axis2 = math.sqrt(0.5*(a+b-delta)) if a+b-delta > 0 else 0.
            axis1 = math.sqrt(0.5*(a+b+delta)) if a+b+delta > 0 else 0.
            ptD   = math.sqrt(sum_weight)/sum_pt if sum_weight > 0 else 0. 
            radius= math.sqrt(ave_dr2 - ave_dr*ave_dr) if ave_dr2 - ave_dr*ave_dr>0 else 0.
            self.jetprop_.append({
                     "ptD":ptD,
                     "axis1":axis1,
                     "axis2":axis2,
                     "mult":mult,
                     "mult1GeV":mult1GeV,
                     "dR": radius, #-> axis1^2 + axis2^2
                    })

        ## HARD CODED -> To Be matched
        for dR in [0.2]:
            self.jetsR_[dR]=[]
            sequence=cluster( ip, algo="antikt",ep=False,R=dR)
            jets = sequence.inclusive_jets()
            for i, jet in enumerate(jets):
                p = ROOT.TLorentzVector()
                p.SetPtEtaPhiM(  jet.pt, jet.eta, jet.phi, jet.mass)
                self.jetsR_[dR] .append(p) 

        for dR in [0.1]:
            self.jetsR_[dR]=[]
            for j in self.jetsR_[0.2]:
                #compute a E01 jets
                #ptV,etaV,phiV,energyV
                j01 = ROOT.TLorentzVector()
                if j.Pt()<15:  
                    ## fast exit for jets that are useless
                    self.jetsR_[dR].append(j01)
                    continue
                for k in range(0,len(ptV)):
                    p =ROOT.TLorentzVector()
                    p.SetPtEtaPhiE(ptV[k],etaV[k],phiV[k],energyV[k])
                    if p.DeltaR(j)<dR:
                        j01 += p
                self.jetsR_[dR].append(j01)
        return self

    def do_genjets(self,ptV,etaV,phiV,energyV):
        ''' just grab genjets and put in the class collection'''
        for pt, eta, phi, energy in zip( ptV,etaV,phiV,energyV):
            p = ROOT.TLorentzVector()
            p.SetPtEtaPhiE(  pt, eta, phi, energy)
            self.genjets_ .append(p) 
        return self
   

    def print_progress(self,done,tot,n=30,every=200):
        if done % every != 0  and done < tot-1: return
        a=int(float(done)*n/tot )
        l="\r["
        if a==n: l += "="*n
        elif a> 0: l +=  "="*(a-1) + ">" + " "*(n-a)
        else: l += " "*n
        l+="] %.1f%%"%(float(done)*100./tot)
        if self.time_ == None: self.time_ = datetime.now()
        else: 
            new = datetime.now()
            delta=(new-self.time_)
            H=delta.seconds/3600 
            M=delta.seconds/60 -H*60
            S=delta.seconds - M*60-H*3600
            H+= delta.days*24
            if H>0:
                l+= " in %dh:%dm:%ds "%(H,M,S) 
            else:
                l+= " in %dm:%ds "%(M,S) 
            S=int( (delta.seconds+24*3600*delta.days)*float(tot-done)/float(done) )
            M= S/60
            S-=M*60
            H= M/60
            M-=H*60
            l+= " will end in %dh:%dm:%ds                      "%(H,M,S)
        if a==n: l+="\n"
        print l,
        sys.stdout.flush()
        return self

    def loop(self):
        #todo=[efficiency_pteta(self) ] 
        todo=[]
        print "-> Running on",
        for x in self.cfg.torun:
            print x,
            exec("todo.append("+x+"(self))")
        print 
        
        nentries=self.chain.GetEntries()
        for ientry,entry in enumerate(self.chain):
            self.print_progress(ientry,nentries)
            self.clear()
            c3d_pt_ = np.array(entry.cl3d_pt)
            c3d_eta_ = np.array(entry.cl3d_eta)
            c3d_phi_ = np.array(entry.cl3d_phi)
            c3d_energy_ = np.array(entry.cl3d_energy)

            if "matrix_calibration" in self.cfg.torun or self.doCalibration or 'layer_deposits' in self.cfg.torun:
                c3d_clusters_ = np.array(entry.cl3d_clusters)
                c3d_nclu_ = np.array(entry.cl3d_nclu)

                self.cl_pt_ = np.array(entry.cl_pt)
                self.cl_eta_ = np.array(entry.cl_eta)
                self.cl_phi_ = np.array(entry.cl_phi)
                self.cl_energy_ = np.array(entry.cl_energy)
                self.cl_layer_ = np.array(entry.cl_layer)
                self.cl_ncells_ = np.array(entry.cl_ncells)
                self.cl_cells_ = np.array(entry.cl_cells) ### ?

                ## compute cl raw energy
                self.cl_raw_=[]

                if self.printCalib< 10:
                    print " -----------------2D---------------------------"

                for icl in range(0,len(self.cl_pt_)):
                    layer=self.cl_layer_[icl]
                    raw = max(self.cl_pt_[icl]-self.cut[layer-1],0.) if self.calibration =='inversion' else self.cl_pt_[icl]

                    #self.cl_raw_.append(self.cl_energy_[icl])
                    self.cl_raw_.append(raw)
                    if self.printCalib< 10:
                        print icl,"pt=",self.cl_pt_[icl],"raw=",raw,"\t",raw/self.cl_pt_[icl]
            
            if self. doCalibration and self.calibration=='inversion': 
                #fix 3D -> 2D association#

                #compute calibrated pt
                ptCalib=[]
                for iC in range(0,len(c3d_pt_)):
                    pt3d=0.0
                    #print "DEBUG: 3D Cluster",iC,"has associated N=",c3d_nclu_[iC],"2d clusters"
                    for idx in range(0,c3d_nclu_[iC]):
                        icl = c3d_clusters_[iC][idx]
                        raw = self.cl_raw_[icl]
                        layer= self.cl_layer_[icl]
                        pt3d += self.const[layer-1]*raw
                        #pt3d += self.const[layer-1]*max(raw-self.cut[layer-1],0)
                        #print "DEBUG:",idx,") cl=",icl,"pt_len=",len(self.cl_phi_),"raw=",len(self.cl_raw_)

                    ptCalib.append(pt3d)

                ## update pt and energy using ptCalib
                if self.printCalib< 10:
                    print "--------- CALIBRATION 3D----------"

                for iC in range(0,len(c3d_pt_)):
                    if self.printCalib< 10:
                        print iC,c3d_pt_[iC],c3d_eta_[iC],"->\t", ptCalib[iC],"\t|", ptCalib[iC]/c3d_pt_[iC]
                    c3d_pt_[iC] = ptCalib[iC]



            genjets_pt_ = np.array(entry.genjet_pt)
            genjets_eta_ = np.array(entry.genjet_eta)
            genjets_phi_ = np.array(entry.genjet_phi)
            genjets_energy_ = np.array(entry.genjet_energy)

            ## produce trigger jets and gen jets
            self.do_clusters(c3d_pt_,c3d_eta_,c3d_phi_,c3d_energy_).do_genjets( genjets_pt_,genjets_eta_,genjets_phi_,genjets_energy_)

            if ientry< 10:
                print "ENTRY",ientry
                print "JETS", len(self.jets_)
                for dR in self.jetsR_:
                    print "JETS R = ",dR,len(self.jetsR_[dR])
                print "---------------------------------"

            if self.doCalibration and self.pu != "":
                if self.printCalib< 10: print "--------- PU CALIBRATION ----------"
                for ijet, jet in enumerate(self.jets_):
                    eta= abs(jet.Eta())
                    pt = jet.Pt()
                    b = self.pPU.FindBin(eta)
                    pu = self.pPU.GetBinContent(b)
                    corr =  max( (pt - pu) / pt,0.)
                    jet *= corr
                    if self.printCalib< 10: 
                        print ijet,"JetPt",pt,"pu",pu, "->",jet.Pt()

            if self.doCalibration and self.calibf != None:
                if self.printCalib< 10: print "--------- JET CALIBRATION ----------"
                for ijet, jet in enumerate(self.jets_):
                    if self.printCalib< 10:
                        print ijet,"JetPt",jet.Pt(), "->",
                    if self.calibf != None:
                        jet *= self.calibf.Eval(jet.Pt())
                    if self.printCalib< 10:
                        print jet.Pt()

            if self.doCalibration:
                if self.printCalib< 10:
                    self.printCalib += 1
                    print "--------------------------------"

            ## produce plots
            for x in todo: x.run()
        self.output.cd()
        for x in todo: x.finalize()
        return self
