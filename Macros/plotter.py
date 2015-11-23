from ROOT import *
import math

########################################

draw_pulls = True
draw_residuals = True
draw_resolution = True

########################################

# list of files, if more than one plots will be superimposed

list_files = [] #arguments: root file name, name of tree including path, name for legend (optional), color (optional)  
list_files.append(["/Users/simoniel/workspace/files_root/tracking/aida_clic_o2_v03_mu_15deg_5GeV_1000evts.root","MyRecoMCTruthLinker_Refit/perftree","test1",kBlue])
# list_files.append(["/Users/simoniel/workspace/files_root/tracking/aida_clic_o2_v03_mu_15deg_5GeV_1000evts_testmacro.root","MyRecoMCTruthLinker_Refit/perftree","test2",kRed])


# parameters for pull plots

nbins_pulls = 50
range_pulls = 10. 


# parameters for residual plots

nbins_residuals = 100
range_residuals_omega = 0.0002 
range_residuals_phi = 0.1 
range_residuals_tanl = 0.1 
range_residuals_d0 = 5 
range_residuals_z0 = 5 


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
#gStyle.SetOptFit(1);

########################################


def GetHisto(f, v, cond=None):

    if not isinstance(cond,str):
        cond = ""
    
    file = ROOT.TFile.Open(f[0], "read")
    tree = file.Get(f[1])

    if v[1] > 0 :
        histo = TH1D( "histo", "histo", v[1], v[2], v[3] )
        
    tree.Project("histo",v[0],cond)
    histo = gDirectory.GetList().FindObject("histo")

    if len(f) > 3:
        histo.SetLineColor(f[3])
        histo.SetMarkerColor(f[3])

    return histo
        


def DoGaussFit(f, v, histo=None):

    if not isinstance(histo,TH1D):
        histo = GetHisto(f, v)

    histo.Fit( "gaus","","", v[2], v[3] )
    histo.SetTitle( v[0] )    
    
    fitsf = histo.GetFunction("gaus")
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

    hists = []
    for file in list_files:
        for var in list_var:

            c_data_pulls.cd(list_var.index(var)+1)
            hists.append(DoGaussFit(file, var)[0])
            # hists[-1].SetLineColor(list_files.index(file)+1)
            hists[-1].Draw("same")

    c_data_pulls.Print("pulls.eps")


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
    for file in list_files:
        for var in list_var:

            c_data_residuals.cd(list_var.index(var)+1)
            hists.append(GetHisto(file, var))
            # hists[-1].SetLineColor(list_files.index(file)+1)
            hists[-1].Draw("same")

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
            if len(file) > 3:
                graphs[-1].SetLineColor(file[3])
                graphs[-1].SetMarkerColor(file[3])
            else:
                graphs[-1].SetLineColor(list_files.index(file)+1)
                graphs[-1].SetMarkerColor(list_files.index(file)+1)
            if len(file) > 2:
                leg.AddEntry(graphs[-1], file[2], "lp");
            else:
                leg.AddEntry(graphs[-1], "test "+str(list_files.index(file)+1), "lp");
            mg.Add(graphs[-1])

        mg.Draw("AP")
        leg.Draw("same")

        title = "resolution_" +str(list_resolution_var.index(var)+1)+".eps"
        c_data_resolution.Print(title)


########################################












