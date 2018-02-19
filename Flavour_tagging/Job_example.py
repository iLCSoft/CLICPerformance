from DIRAC.Core.Base import Script
Script.parseCommandLine()
from ILCDIRAC.Interfaces.API.DiracILC import DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *
  
dirac = DiracILC( True , "CLICdet_repo.cfg" )
job = UserJob( )
job.setName( "example" )
job.setInputData(["/ilc/user/p/proloff/stdhep/2f_samples_23_04_2013/whizard_SM_bb_500_90deg.stdhep"])
job.setCLICConfig("ILCSoft-2017-12-21")
job.setOutputSandbox( ["*.root","*.slcio","*.log"] )


ddsim = DDSim()
ddsim.setSteeringFile("CLICPerformance/examples/clic_steer.py")
ddsim.setVersion("ILCSoft-2017-12-21_gcc62")
ddsim.setDetectorModel("CLIC_o3_v14")
ddsim.setInputFile("/ilc/user/p/proloff/stdhep/2f_samples_23_04_2013/whizard_SM_bb_500_90deg.stdhep")
ddsim.setOutputFile("ddsim_output.slcio")
ddsim.setStartFrom(0)
ddsim.setNumberOfEvents(10)
res1 = job.append(ddsim)
if not res1['OK']:
   print res1['Message']
   sys.exit(2)


# -------------------- Comment if gg_had overlay not wanted -----------------------------------------------------#
over = OverlayInput()
over.setBXOverlay(30)
over.setGGToHadInt(0.3)
over.setNumberOfSignalEventsPerJob(100)
over.setBackgroundType("gghad")
over.setPathToFiles( "/ilc/prod/clic/500gev/gghad/CLIC_o3_v13/SIM/00008670/" )
res = job.append(over)
if not res['OK']:
   print res['Message']
   exit()


clicReco = Marlin()
clicReco.setVersion("ILCSoft-2017-12-21_gcc62")
clicReco.getInputFromApp(ddsim)
clicReco.setSteeringFile("CLICPerformance/examples/clicReconstruction.xml")
clicReco.setDetectorModel("CLIC_o3_v14")
clicReco.setOutputFile("Output_DST.slcio")
res = job.append(clicReco)
if not res['OK']:
    print res['Message']
    exit()
   
   
makentp = Marlin()
makentp.setVersion("ILCSoft-2017-12-21_gcc62")
makentp.getInputFromApp(clicReco)
makentp.setSteeringFile("vtx_makentp.xml")
makentp.setDetectorModel("CLIC_o3_v14")
makentp.setOutputFile("Output_Ntuple_xx.root")
res = job.append(makentp)
if not res['OK']:
    print res['Message']
    exit()
   
job.dontPromptMe()
   
job.submit(dirac)
   
