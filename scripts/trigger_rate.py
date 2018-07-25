import sys
oldargv = sys.argv[:]
import ROOT
sys.argv = oldargv
from glob import glob

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
events = Events("file:test.root")

/eos/cms/store/mc/Upg2017Summer15DR/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/GEN-SIM-RECO
