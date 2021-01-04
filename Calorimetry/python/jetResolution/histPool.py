"""Helper class which contains all histograms."""
import ROOT
import math
import array
from helperFunctions import randomString
from helperFunctions import calculateRMS90
from helperFunctions import S_OK, S_ERROR
from collections import namedtuple


VectorTuple = namedtuple('VectorTuple', "quark1 quark2 totalTrueVisible "
                         "totalTrueInvisible totalPFO genJets recoJets")


class histPool(object):
  """Helper class which contains all histograms."""

  cosThetaRange = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925,
                   0.95, 0.975, 1.0]
  histInfo = [['TH1F', 'cosTheta', 'cosTheta; cos(#theta); Entries',
              50, 0, 1.0],
              ['TH1F', 'true_totalE', '; total true E [GeV]; Entries',
              400, 200, 300],
              ['TH1F', 'gen_totalE', '; gen di-jet E [GeV]; Entries',
              250, 0, 250],
              ['TH1F', 'reco_totalE', '; reco di-jet E [GeV]; Entries',
              250, 0, 250],
              ['TH1F', 'gen_InvMass', '; gen inv.mass [GeV]; Entries',
              250, 0, 250],
              ['TH1F', 'reco_InvMass', '; reco inv.mass [GeV]; Entries',
              250, 0, 250],
              ['TH1F', 'JER', 'cosTheta; cos(#theta); JER [%]',
              len(cosThetaRange) - 1, cosThetaRange]
  ]

  def __init__(self):
    """Initialize."""
    self.histDict2D = {}

  def initNewSetOfHists(self, setName, histInformation):
    """Initialize new set of histograms."""
    if setName in self.histDict2D:
      return True
    tmpDict = {}
    for iHistInfo in histInformation:
      if iHistInfo[0] == 'TH1F':
        if len(iHistInfo) == 5:
          tmpArr = array.array('d', iHistInfo[4])
          tmpHist = ROOT.TH1F(randomString(10), iHistInfo[2], iHistInfo[3],
                              tmpArr)
        if len(iHistInfo) == 6:
          tmpHist = ROOT.TH1F(randomString(10), iHistInfo[2], iHistInfo[3],
                              iHistInfo[4], iHistInfo[5])
        tmpDict[iHistInfo[1]] = tmpHist
    self.histDict2D[setName] = tmpDict

  def createJERFitHistInfo(self):
    """???"""
    outHistInfo = []
    histInfoTemplate = ['TH1F', 'cosThetaBin', '; reco + true neutrino energy [GeV]; Entries',
                        10000, 0., 500.]
    for i in range(1, len(self.cosThetaRange)):
      thetaRange = 'cosTheta: ' + str(self.cosThetaRange[i - 1]) + ' - ' + str(self.cosThetaRange[i])
      binNumber = i - 1
      histInfoEntry = list(histInfoTemplate)
      histInfoEntry[1] += str(binNumber)
      histInfoEntry[2] = thetaRange + histInfoEntry[2]
      outHistInfo.append(histInfoEntry)
    return outHistInfo

  def fillHist(self, histName, val, histSetName=''):
    """Fill histogram.

    Return True if val is successfully added to histogrm, otherwise return
    False.
    """
    if histSetName not in self.histDict2D:
      self.initNewSetOfHists(histSetName, self.histInfo)
      jerFitHistInfo = self.createJERFitHistInfo()
      fitHistSetName = histSetName + '/fit' if histSetName else 'fit'
      self.initNewSetOfHists(fitHistSetName, jerFitHistInfo)

    if histName in self.histDict2D[histSetName]:
      self.histDict2D[histSetName][histName].Fill(val)
      return True
    else:
      print("Can't find the hitogram: %s" % histName)
      return False

  def fillAllHists(self, vTuple, histSetName=''):
    """Fill all histograms."""
    cosTheta = abs(vTuple.quark1.Pz() / vTuple.quark1.P())
    self.fillHist('cosTheta', cosTheta, histSetName)

    totalTrueE = vTuple.totalTrueVisible + vTuple.totalTrueInvisible
    self.fillHist('true_totalE', totalTrueE.E(), histSetName)

    genJetSum = vTuple.genJets[0] + vTuple.genJets[1]
    self.fillHist('gen_InvMass', genJetSum.M(), histSetName)
    self.fillHist('gen_totalE', genJetSum.E(), histSetName)

    recoJetSum = vTuple.recoJets[0] + vTuple.recoJets[1]
    self.fillHist('reco_InvMass', recoJetSum.M(), histSetName)
    self.fillHist('reco_totalE', recoJetSum.E(), histSetName)

    for i in range(1, len(self.cosThetaRange)):
      if cosTheta > self.cosThetaRange[i - 1] and cosTheta < self.cosThetaRange[i]:
        binNumber = i - 1
        fillVal = vTuple.recoJets[0].E() + vTuple.totalTrueInvisible.E()
        fitHistSetName = histSetName + '/fit' if histSetName else 'fit'
        self.fillHist('cosThetaBin' + str(binNumber), fillVal, fitHistSetName)

  def makeJERPlot(self, histSetName=''):
    """Create JER plot."""
    if histSetName not in self.histDict2D:
      return S_ERROR("Hist set with following name is not found: %s" % histSetName)

    for i in range(1, len(self.cosThetaRange)):
      binNumber = i - 1
      fitHistSetName = histSetName + '/fit' if histSetName else 'fit'
      histName = 'cosThetaBin' + str(binNumber)
      hist = self.histDict2D[fitHistSetName][histName]
      if not hist:
        return S_ERROR("No hist found. SetName: %s; histName: %s" % (fitHistSetName, histName))
      if hist.GetEntries() == 0:
        return S_ERROR("Hist is empty. SetName: %s; histName: %s" % (fitHistSetName, histName))

      #  print('Calculate RMS90. histName: %s; histTitle: %s, nEntries: %d, mean: %f; RMS: %f' % (hist.GetName(), hist.GetTitle(), hist.GetEntries(), hist.GetMean(), hist.GetRMS()))
      rms90Data = calculateRMS90(hist)
      resolution = rms90Data['rms90']
      resolutionError = resolution / math.sqrt(rms90Data['totalCounts'])

      self.histDict2D[histSetName]['JER'].SetBinContent(i, resolution)
      self.histDict2D[histSetName]['JER'].SetBinError(i, resolutionError)

    return S_OK()

  def saveToFile(self, fileName):
    """Save all histograms to the file."""
    tmpFile = ROOT.TFile(fileName, "RECREATE")
    for histSetName in self.histDict2D:
      histDict = self.histDict2D[histSetName]
      if histSetName != '':
        tmpFile.mkdir(histSetName)
        tmpFile.cd(histSetName)
      for iHistName, iHist in histDict.items():
        if iHist.GetEntries() != 0:
          iHist.SetName(iHistName)
          iHist.Write()
      tmpFile.cd()
    tmpFile.Close()
