from ROOT import *
import sys
import math
 
########################################

draw_pulls = True
draw_residuals = True
draw_resolution = False

########################################

# list of files, if more than one plots will be superimposed

list_files = [] #arguments: root file name, name of tree including path, name for legend (optional), color (optional)  
list_files.append(["../files_root/tracking/nocombact_refit/aida_clic_o2_v04_mu-_85deg_5GeV_5000evts.root.root","MyRecoMCTruthLinker_Refit/perftree","#theta = 85#circ",kBlue])
list_files.append(["../files_root/tracking/nocombact_refit/aida_clic_o2_v04_mu-_75deg_5GeV_5000evts.root.root","MyRecoMCTruthLinker_Refit/perftree","#theta = 75#circ",kOrange+7])
list_files.append(["../files_root/tracking/nocombact_refit/aida_clic_o2_v04_mu-_65deg_5GeV_5000evts.root.root","MyRecoMCTruthLinker_Refit/perftree","#theta = 65#circ",kGreen+1])

# parameters for pull plots

nbins_pulls = 50
range_pulls = 10. 


# parameters for residual plots

nbins_residuals = 100
range_residuals_omega = 0.000005 
range_residuals_phi = 0.002 
range_residuals_tanl = 0.002 
range_residuals_d0 = 0.05 
range_residuals_z0 = 0.05 


#parameters for resolution plots

logx = True
list_resolution_var = [] # name of var_y, name of var_x, min value of var_x, max value of var_x, npoints in x, log/linear scale for x
list_resolution_var.append(["(recoPt-truePt)/(truePt*truePt)", "truePt", 0.1, 1000., 10, logx])
list_resolution_var.append(["(recoPt-truePt)/(truePt*truePt)", "trueTheta", 0, 90., 10, not logx])


########################################

gStyle.SetOptFit(111111);
gStyle.SetOptStat(111111);
#gROOT.ProcessLine(".x style/CLICdpSettingsScript.C");
#gStyle.SetOptStat(0000);
#gStyle.SetOptFit(0000);
#gStyle.SetOptFit(1);

########################################


def CheckInput(list):
    for file in list:
        if len(file) < 2:
            print "ERROR: file and/or tree name not given"
            sys.exit(1)
        elif len(file)<3:
            file.append("test "+str(list.index(file)+1))
        if len(file)<4:
            file.append(list.index(file)+1)
  


def GetHisto(f, v, cond=None):

    if not isinstance(cond,str):
        cond = ""
    
    file = ROOT.TFile.Open(f[0], "read")
    tree = file.Get(f[1])

    if v[1] > 0 :
        histo = TH1D( "histo", "histo", v[1], v[2], v[3] )
        
    tree.Project("histo",v[0],cond)
    histo = gDirectory.GetList().FindObject("histo")

    histo.SetTitle( v[0] ) 
    histo.SetLineColor(f[3])
    histo.SetMarkerColor(f[3])

    return histo
        


def DoGaussFit(f, v, histo=None):

    if not isinstance(histo,TH1):
        histo = GetHisto(f, v)

    fitsf = TF1("fit","gaus",v[2],v[3])
    histo.Fit(fitsf,"N","",v[2],v[3])

    fitsf.SetLineColor(histo.GetLineColor())

    return histo ,fitsf



def GetRangeRMSfrac(histo, frac=0.90):

    integral_all = histo.Integral()
    integral_frac = integral_all*frac

    # include under/overflows
    first_bin = 0
    last_bin = histo.GetNbinsX()+1

    # maximum values allowed
    bin_min = first_bin
    bin_max = last_bin
    nbins_min = last_bin-first_bin+1
    diff_min = nbins_min

    # compute rms90: smallest bin range with integral >= integral_frac, if there is more than one range with the same number of bins take the range where the mean value is more centered
    for start in range(first_bin, last_bin):
        integral = 0
        for i in range(start,last_bin):
            integral += histo.GetBinContent(i)
            if integral >= integral_frac:
                nbins = i - start + 1
                bin_of_mean = histo.FindBin(histo.GetMean())
                diff = math.fabs(nbins/2.-bin_of_mean)
                if nbins < nbins_min:
                    nbins_min = nbins
                    bin_min = start
                    bin_max = i
                    diff_min = diff
                elif nbins == nbins_min:
                    if diff < diff_min:
                        diff_min = diff
                        nbins_min = nbins
                        bin_min = start
                        bin_max = i

    return bin_min, bin_max
    
    
    
def GetResolution(f, v):

    # v array:
    #    0      1      2      3       4      5
    # [var_y, var_x, min_x, max_x, npoints, logx]

    graph = TGraphErrors()

    # step and cond computed differently according to log / not log x scale in order to always have equispaced bins 

    if v[5]:
        if v[2]==0:
            v[2] = 0.1
            print "ERROR: log scale does not allow 0 as min value - min now set to 0.1"
        step = math.pow(v[3]/v[2], 1./v[4])
    else:
        step = (v[3]-v[2])/v[4]
    
    
    lowedge_step = v[2]

    for i in range(0,v[4]):

        if v[5]:
            upedge_step = lowedge_step*step
            cond = v[1] + " >= " + str(lowedge_step) + " && "  + v[1] + " < " + str(upedge_step)
            lowedge_step = upedge_step
        else:
            cond = v[1] + " >= " + str(v[2]+i*step) + " && "  + v[1] + " < " + str(v[2]+(i+1)*step)
        print "cond = ", cond
        
        histo_x = GetHisto(f, [v[1], -1, 0, 0])
        
        if histo_x.GetEntries()>0: 
            x = histo_x.GetMean()
            x_err = histo_x.GetMeanError()

            # in two steps because very large outlayers can cause to have all the distribution in one or very few bins biasing the gaussian fit
            helper_histo_y = GetHisto(f, [v[0], -1, 0, 0])
            min_bin, max_bin = GetRangeRMSfrac(helper_histo_y)
            y_min = helper_histo_y.GetBinLowEdge(1)+(min_bin-1)*helper_histo_y.GetBinWidth(1)
            y_max = helper_histo_y.GetBinLowEdge(1)+max_bin*helper_histo_y.GetBinWidth(1)
            
            nbins = 100
            histo_y, fit_y = DoGaussFit(f, [v[0], nbins, y_min, y_max], cond)
            y = fit_y.GetParameter(2)
            y_err = fit_y.GetParError(2)
            graph.SetPoint(i,x,y)
            graph.SetPointError(i,x_err,y_err)
           
        else:
            graph.SetPoint(i,0,0)
            graph.SetPointError(i,0,0)

    return graph
        
    
########################################
########################################

CheckInput(list_files)

########################################

if draw_pulls:
    
    list_var = []
    list_var.append(["pullOmega", nbins_pulls, -range_pulls, range_pulls])
    list_var.append(["pullPhi", nbins_pulls, -range_pulls, range_pulls])
    list_var.append(["pullTanLambda", nbins_pulls, -range_pulls, range_pulls])
    list_var.append(["pullD0", nbins_pulls, -range_pulls, range_pulls])
    list_var.append(["pullZ0", nbins_pulls, -range_pulls, range_pulls])
    
    
    c_data_pulls = TCanvas( 'c_data_pulls', 'pulls', 200, 10, 800, 600 )
    npad = math.sqrt(len(list_var))
    c_data_pulls.Divide(int(math.ceil(npad)), int(math.floor(npad)))

    tab = [[None for x in range(len(list_files)+1)] for x in range(len(list_var)+1)] 
    hists = []
    fits = []
    legs = []
    for var in list_var:
        leg = TLegend(0.15, 0.88-0.04*len(list_files), 0.25, 0.88)
        leg.SetBorderSize(0)
        legs.append(leg)
        c_data_pulls.cd(list_var.index(var)+1)
        for file in list_files:
            h, f = DoGaussFit(file, var)
            hists.append(h)
            fits.append(f)
            hists[-1].Draw("same")
            fits[-1].Draw("same")
            legs[-1].AddEntry(hists[-1], file[2], "l")
            tab[0][list_files.index(file)+1] = file[2]
            tab[list_var.index(var)+1][0] = var[0]
            tab[list_var.index(var)+1][list_files.index(file)+1] = str(fits[-1].GetParameter(1))+" +- "+str(fits[-1].GetParError(1)) +"   "+str(fits[-1].GetParameter(2))+" +- "+str(fits[-1].GetParError(2))
        legs[-1].Draw("same")
    c_data_pulls.Print("pulls.eps")
    
    outfile = open("tab_pulls.txt", "w")
    for row in tab:
        outfile.write(str([str(x) for x in row])+"\n")
    outfile.close()
    


 ########################################
   
    
if draw_residuals:
    
    list_var = []
    list_var.append(["resOmega", nbins_residuals, -range_residuals_omega, range_residuals_omega])
    list_var.append(["resPhi", nbins_residuals, -range_residuals_phi, range_residuals_phi])
    list_var.append(["resTanLambda", nbins_residuals, -range_residuals_tanl, range_residuals_tanl])
    list_var.append(["resD0", nbins_residuals, -range_residuals_d0, range_residuals_d0])
    list_var.append(["resZ0", nbins_residuals, -range_residuals_z0, range_residuals_z0])

    c_data_residuals = TCanvas( 'c_data_residuals', 'residuals', 200, 10, 800, 600 )
    npad = math.sqrt(len(list_var))
    c_data_residuals.Divide(int(math.ceil(npad)), int(math.floor(npad)))

    hists = []
    legs = []
    for var in list_var:
        leg = TLegend(0.15, 0.88-0.04*len(list_files), 0.25, 0.88)
        leg.SetBorderSize(0)
        legs.append(leg)
        c_data_residuals.cd(list_var.index(var)+1)
        for file in list_files:            
            hists.append(GetHisto(file, var))
            hists[-1].Draw("same")
            legs[-1].AddEntry(hists[-1], file[2], "l");
        legs[-1].Draw("same")

    c_data_residuals.Print("residuals.eps")


########################################
    
    
if draw_resolution:


    for var in list_resolution_var:

        c_data_resolution = TCanvas( 'c_data_resolution', 'resolution', 200, 10, 800, 600 )

        graphs = []
        mg = TMultiGraph()
        
        leg = TLegend(0.70, 0.88-0.04*len(list_files), 0.88, 0.88)
        leg.SetBorderSize(0)
    
        for file in list_files:

            graphs.append(GetResolution(file, var))
            graphs[-1].SetMarkerSize(1.)
            graphs[-1].SetMarkerStyle(20)
            graphs[-1].SetLineColor(file[3])
            graphs[-1].SetMarkerColor(file[3])
            leg.AddEntry(graphs[-1], file[2], "lp");
            mg.Add(graphs[-1])

        mg.Draw("AP")
        leg.Draw("same")

        title = "resolution_" +str(list_resolution_var.index(var)+1)+".eps"
        c_data_resolution.Print(title)


########################################












