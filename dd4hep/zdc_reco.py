'''
    An example option file to digitize/reconstruct/clustering calorimeter hits
'''
from Gaudi.Configuration import *
import json
import os
import ROOT

from Configurables import ApplicationMgr, EICDataSvc, PodioInput, PodioOutput, GeoSvc
from GaudiKernel.SystemOfUnits import MeV, GeV, mm, cm, mrad

#detector_name = str(os.environ.get("JUGGLER_DETECTOR", "endcapP_insert"))
#detector_path = str(os.environ.get("DETECTOR_PATH", "."))


if "JUGGLER_DETECTOR_NAME" in os.environ:
    detector_name = os.environ["JUGGLER_DETECTOR_NAME"]
elif "staggered_h3" in os.environ["JUGGLER_SIM_FILE"]:
    detector_name = "zdc_staggered_h3"
elif "staggered_h4" in os.environ["JUGGLER_SIM_FILE"]:
    detector_name = "zdc_staggered_h4"
elif "unstaggered" in os.environ["JUGGLER_SIM_FILE"]:
    detector_name = "zdc_unstaggered"
elif "s1" in os.environ["JUGGLER_SIM_FILE"]:
    detector_name = "zdc_s1_inf"
elif "s2" in os.environ["JUGGLER_SIM_FILE"]:
    detector_name = "zdc_s2_inf"
detector_path = "."
compact_path = os.path.join(detector_path, detector_name)

# input arguments from calibration file
#with open(f'{detector_path}/calibrations/emcal_barrel_calibration.json') as f:
    #calib_data = json.load(f)['electron']

#print(calib_data)

#cb_ecal_sf = float(calib_data['sampling_fraction_img'])
#scifi_barrel_sf = float(calib_data['sampling_fraction_scfi'])

# get sampling fractions from system environment variable, 1.0 by default
#ci_ecal_sf = float(os.environ.get("CI_ECAL_SAMP_FRAC", 0.253))
#cb_hcal_sf = float(os.environ.get("CB_HCAL_SAMP_FRAC", 0.038))
#ci_hcal_sf = float(os.environ.get("CI_HCAL_SAMP_FRAC", 0.025))
#ce_hcal_sf = float(os.environ.get("CE_HCAL_SAMP_FRAC", 0.025))

#ci_hcal_sf = "1."
ci_zdc_sf = "1."
#ci_ecal_sf = "0.03"
#ci_ecal_insert_sf = "0.03"

# input and output
input_sims = [f.strip() for f in str.split(os.environ["JUGGLER_SIM_FILE"], ",") if f.strip()]
output_rec = str(os.environ["JUGGLER_REC_FILE"])
n_events = int(os.environ["JUGGLER_N_EVENTS"])

# geometry service
geo_service = GeoSvc("GeoSvc", detectors=["{}.xml".format(compact_path)], OutputLevel=INFO)
# data service
podioevent = EICDataSvc("EventDataSvc", inputs=input_sims)


# juggler components
from Configurables import Jug__Digi__CalorimeterHitDigi as CalHitDigi
from Configurables import Jug__Reco__CalorimeterHitReco as CalHitReco
# from Configurables import Jug__Reco__CalorimeterHitsMerger as CalHitsMerger
# from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster

# from Configurables import Jug__Reco__ImagingPixelReco as ImCalPixelReco
# from Configurables import Jug__Reco__ImagingTopoCluster as ImagingCluster

# from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG
# from Configurables import Jug__Reco__ImagingClusterReco as ImagingClusterReco

from Configurables import Jug__Fast__InclusiveKinematicsTruth as InclusiveKinematicsTruth

# branches needed from simulation root file
sim_coll = [
    "MCParticles",
#    "HcalEndcapPHits",
#    "HcalEndcapPHitsContributions",
    "ZDCHits",
    "ZDCHitsContributions",
#    "EcalEndcapPHits",
#    "EcalEndcapPHitsContributions",
#    "EcalEndcapPInsertHits",
#    "EcalEndcapPInsertHitsContributions"
]

# input and output
podin = PodioInput("PodioReader", collections=sim_coll)
podout = PodioOutput("out", filename=output_rec)
'''
# Hadron endcap HCal
ci_hcal_daq = dict(
         dynamicRangeADC=200.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_hcal_digi = CalHitDigi("ci_hcal_digi",
         inputHitCollection="HcalEndcapPHits",
         outputHitCollection="HcalEndcapHitsDigi",
         **ci_hcal_daq)
ci_hcal_reco = CalHitReco("ci_hcal_reco",
        inputHitCollection=ci_hcal_digi.outputHitCollection,
        outputHitCollection="HcalEndcapPHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_hcal_sf,
        **ci_hcal_daq)
'''
# ZDC
ci_zdc_daq = dict(
         dynamicRangeADC=800.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_zdc_digi = CalHitDigi("ci_zdc_digi",
         inputHitCollection="ZDCHits",
         outputHitCollection="ZDCHitsDigi",
         **ci_zdc_daq)
ci_zdc_reco = CalHitReco("ci_zdc_reco",
        inputHitCollection=ci_zdc_digi.outputHitCollection,
        outputHitCollection="ZDCHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_zdc_sf,
        **ci_zdc_daq)

# Truth level kinematics
truth_incl_kin = InclusiveKinematicsTruth("truth_incl_kin",
        inputMCParticles = "MCParticles",
        outputInclusiveKinematics = "InclusiveKinematicsTruth"
)

# Output
podout.outputCommands = ['drop *',
        'keep MCParticles',
        'keep *Digi',
        'keep *Reco*',
	'keep Inclusive*']

ApplicationMgr(
    TopAlg = [podin,
            #ci_hcal_digi, ci_hcal_reco,
            ci_zdc_digi, ci_zdc_reco,
            #ci_ecal_digi, ci_ecal_reco,
            #ci_ecal_insert_digi, ci_ecal_insert_reco,
	    truth_incl_kin, podout],
    EvtSel = 'NONE',
    EvtMax = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
)
