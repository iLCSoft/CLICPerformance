"""Jet studies."""

# !/bin/python

import os
import sys
import math
import ROOT
from ROOT import gStyle
from jetResolution.helperFunctions import checkJetPairMatching
from jetResolution.helperFunctions import checkForRequiredFields, getTree, getLorentzVectors
from jetResolution.histPool import histPool, VectorTuple

#  ROOT.gROOT.LoadMacro("CLIC_style/CLICdpStyle/rootstyle/CLICdpStyle.C+")
#  ROOT.CLICdpStyle()


def main(inFileName):
  """Fill plots."""
  gStyle.SetOptStat(0)

  inFile = ROOT.TFile(inFileName, "READ")
  if not inFile.IsOpen():
    return ("Can't find the file: %s" % inFileName)

  res = getTree(inFile)
  if not res['OK']:
    return res['Message']
  inTree = res['Value']

  requiredBranches = ['d1_mcPDGID', 'd1_mcE', 'd1_mcPx', 'd1_mcPy', 'd1_mcPz',
                      'd2_mcPDGID', 'd2_mcE', 'd2_mcPx', 'd2_mcPy', 'd2_mcPz',
                      'E_trueAll', 'Px_trueAll', 'Py_trueAll', 'Pz_trueAll',
                      'E_trueInv', 'Px_trueInv', 'Py_trueInv', 'Pz_trueInv',
                      'E_totPFO', 'Px_totPFO', 'Py_totPFO', 'Pz_totPFO',
                      'genJetE', 'genJetPx', 'genJetPy', 'genJetPz',
                      'recoJetE', 'recoJetPx', 'recoJetPy', 'recoJetPz']
  treeBranches = [x.GetName() for x in inTree.GetListOfBranches()]
  res = checkForRequiredFields(requiredBranches, treeBranches,
                               "Tree doesn't contains all required branches.")
  if not res['OK']:
    return res['Message']

  counter = 0
  nEventsToProcess = sys.maxint

  hPool = histPool()

  for iEv in inTree:
    if counter >= nEventsToProcess:
      break
    counter += 1

    vTuple = VectorTuple(
        quark1=getLorentzVectors(iEv.d1_mcPx, iEv.d1_mcPy, iEv.d1_mcPz, iEv.d1_mcE),
        quark2=getLorentzVectors(iEv.d2_mcPx, iEv.d2_mcPy, iEv.d2_mcPz, iEv.d2_mcE),
        totalTrueVisible=getLorentzVectors(iEv.Px_trueAll, iEv.Py_trueAll, iEv.Pz_trueAll, iEv.E_trueAll),
        totalTrueInvisible=getLorentzVectors(iEv.Px_trueInv, iEv.Py_trueInv, iEv.Pz_trueInv, iEv.E_trueInv),
        totalPFO=getLorentzVectors(iEv.Px_totPFO, iEv.Py_totPFO, iEv.Pz_totPFO, iEv.E_totPFO),
        genJets=getLorentzVectors(iEv.genJetPx, iEv.genJetPy, iEv.genJetPz, iEv.genJetE),
        recoJets=getLorentzVectors(iEv.recoJetPx, iEv.recoJetPy, iEv.recoJetPz, iEv.recoJetE)
        )

    if not checkJetPairMatching(vTuple.genJets, vTuple.recoJets, 10./180.*math.pi):
      hPool.fillAllHists(vTuple, 'failJetMatching')
      continue

    hPool.fillAllHists(vTuple, 'signalSelection')
    
  hPool.makeJERPlot('signalSelection')
  hPool.saveToFile('hists_%s' % os.path.basename(inFileName))
  return("Successful execution.")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        print ("\n[INFO]\tRead file: {0}".format(sys.argv[1]))
        print(main(sys.argv[1]))
    else:
        print ("\n[ERROR]\tProvide input file. Terminating...")
        sys.exit()
