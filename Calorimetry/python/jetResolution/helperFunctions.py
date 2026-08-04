"""Set of helper functions."""

import re
import string
import random
import ROOT
import math
import sys
from ROOT import TLorentzVector as lorVec


def S_OK(value=None):
  """Return value on success.

  :param value: value of the 'Value'
  :return: dictionary { 'OK' : True, 'Value' : value }
  """
  return {'OK': True, 'Value': value}


def S_ERROR(messageString=''):
  """Return value on error condition.

  :param string messageString: error description
  """
  return {'OK': False, 'Message': '[ERROR]\t%s' % str(messageString)}


def printSet(inSet):
  """Print set."""
  if isinstance(inSet, set):
    return re.split(r'\(|\)', str(inSet))[1]
  else:
    return ''


def checkForRequiredFields(requiredFields, providedFields, errMsg=""):
  """Check for required fields."""
  sDiff = set(requiredFields) - set(providedFields)
  if len(sDiff) > 0:
    errMsg = ("%s Missing branches: %s"
              % (errMsg, list(sDiff)))
    return S_ERROR(errMsg)
  return S_OK()


def getTree(inFile):
  """
  Get TTree object from TFile.

  If there are no or more than one TTree objects return error.
  """
  keyList = inFile.GetListOfKeys()
  outTree = None
  for iKey in keyList:
    iObj = inFile.Get(iKey.GetName())
    if type(iObj) != ROOT.TTree:
      continue
    if not outTree:
      outTree = iObj
    else:
      return S_ERROR("more than one TTree object is found")
  if outTree is not None:
    return S_OK(outTree)
  else:
    return S_ERROR("no TTree object found in the file: %s" % inFile.GetName())


def getLorentzVectors(px, py, pz, E):
  """Create Lorentx vector."""
  inputArgs = (px, py, pz, E)
  if sum([isinstance(x, float) for x in inputArgs]) == 4:
    outVec = lorVec()
    outVec.SetPxPyPzE(px, py, pz, E)
    return outVec
  if sum([isinstance(x, ROOT.vector("float")) for x in inputArgs]) == 4:
    output = []
    for iCount in range(0, px.size()):
      outVec = lorVec()
      outVec.SetPxPyPzE(px[iCount], py[iCount], pz[iCount], E[iCount])
      output.append(outVec)
    return output


def checkJetPairMatching(jetCol1, jetCol2, matchingAngle):
  """Check jet pair matching.

  Check if jet pair from jetCol1 are matched to pair from jetCol2 within
  provided angle.
  """
  # the second pair will have inverted indices
  possibleJetCombinationsForTheFirstPair = [(0, 0), (0, 1), (1, 0), (1, 1)]
  for iComb in possibleJetCombinationsForTheFirstPair:
    angle1 = jetCol1[iComb[0]].Angle(jetCol2[iComb[1]].Vect())
    angle2 = jetCol1[1 - iComb[0]].Angle(jetCol2[1 - iComb[1]].Vect())
    if angle1 < matchingAngle and angle2 < matchingAngle:
      return True
  return False


def randomString(stringLength=10):
  """Generate a random string of fixed length."""
  letters = string.ascii_lowercase
  return ''.join(random.choice(letters) for i in range(stringLength))


def calculateRMS90(inHist):
  '''???'''

  FLOAT_MAX = sys.float_info.max


  # Calculate raw properties of distribution
  sum = 0.
  total = 0.
  sx = 0.
  sxx = 0.
  nbins = inHist.GetNbinsX()

  for i in range(0, nbins+1):
    binx = inHist.GetBinLowEdge(i) + (0.5 * inHist.GetBinWidth(i))
    yi = inHist.GetBinContent(i)
    sx = sx + yi * binx
    sxx = sxx + yi * binx * binx
    total = total + yi

  rawMean = sx / total
  rawMeanSquared = sxx / total
  rawRms = math.sqrt(rawMeanSquared - rawMean * rawMean)

  sum = 0.
  is0 = 0

  i = 0
  while ((i <= nbins) and (sum < total / 10.)):
    sum = sum + inHist.GetBinContent(i)
    is0 = i
    i = i + 1

  # Calculate truncated properties
  rmsmin = FLOAT_MAX
  sigma = FLOAT_MAX
  sigmasigma = FLOAT_MAX
  frac = FLOAT_MAX
  efrac = FLOAT_MAX
  mean = FLOAT_MAX
  low = FLOAT_MAX
  rms = FLOAT_MAX
  high = 0.

  for istart in range(0, is0 + 1):
    sumn = 0.
    csum = 0.
    sumx = 0.
    sumxx = 0.
    iend = 0

    i = istart
    while ((i <= nbins) and (csum < 0.9 * total)):

      binx = inHist.GetBinLowEdge(i) + (0.5 * inHist.GetBinWidth(i))
      yi = inHist.GetBinContent(i)
      csum += yi

      if (sumn < 0.9 * total):
        sumn += yi
        sumx += yi * binx
        sumxx += yi * binx * binx
        iend = i

      i = i + 1

    localMean = sumx / sumn
    localMeanSquared = sumxx / sumn
    localRms = math.sqrt(localMeanSquared - localMean * localMean)

    if (localRms < rmsmin):
      mean = localMean
      rms = localRms
      low = inHist.GetBinLowEdge(istart)
      high = inHist.GetBinLowEdge(iend)
      rmsmin = localRms

      # resolution:
      # frac = rms / mean * 100.;
      # resolutin uncertainty:
      # efrac = frac / std::sqrt(total);

    rms90 = rmsmin
    mean90 = mean
    totalCounts = total

  print('%s (%d entries), rawrms: %f, rms90: %f (%f-%f), mean90: %f, mean: %f' % (inHist.GetTitle(), inHist.GetEntries(), rawRms, rms90, low, high, mean90, rawMean))
    #  std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries), rawrms: " << rawRms << ", rms90: " << rmsmin
    #                << " (" << low << "-" << high << "), mean90: " << mean << ", mean: " << rawMean;

  outData = {}
  outData['rms90'] = rms90
  outData['mean90'] = mean90
  outData['totalCounts'] = totalCounts

  return outData

    #  for (unsigned int istart = 0; istart <= is0; ++istart)
    #  {
    #
    #      if (localRms < rmsmin)
    #      {
    #          mean = localMean;
    #          rms = localRms;
    #          low = pTH1F->GetBinLowEdge(istart);
    #          high = pTH1F->GetBinLowEdge(iend);
    #          rmsmin = localRms;
    #
    #          // resolution:
    #          // frac = rms / mean * 100.;
    #
    #          //  resolutin uncertainty:
    #          // efrac = frac / std::sqrt(total);
    #      }
    #  }
    #
    #  rms90 = rmsmin;
    #  mean90 = mean;
    #  totalCounts = total;
    #
    #  if (print)
    #  {
    #      std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries), rawrms: " << rawRms << ", rms90: " << rmsmin
    #                << " (" << low << "-" << high << "), mean90: " << mean << ", mean: " << rawMean;
    #  }
