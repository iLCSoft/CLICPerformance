from __future__ import division
from ROOT import *
from array import array
import sys
import math

########################################

draw_pulls = False
draw_residuals = False
draw_resolution = False
draw_check_var = False
draw_hit_residuals = False
draw_efficiency = True
draw_purity = False

########################################

# list of files, if more than one plots will be superimposed

list_files = [] #arguments: root file name, name of tree including path, name for legend (optional), color (optional)

list_files.append(["histograms500GeV_o3_v07_fixed.root", "OuterEndcapPlanarDigiProcessor/trktree", "CLIC_o3_v07", kRed])
list_files.append(["histograms500GeV_o2_v4_ttbar.root", "OuterEndcapPlanarDigiProcessor/trktree", "CLIC_o2_v04", kBlue])


# list_files.append(["/Users/simoniel/workspace/files_root/tracking/HEAD-2016-05-20/merge/aida_test_2016-05-20_clic_o2_v04_mu-_5GeV.root", "MyRecoMCTruthLinker_Forward/perftree", "p = 5 GeV", kAzure+7])
# list_files.append(["/Users/simoniel/workspace/files_root/tracking/HEAD-2016-05-20/merge/aida_test_2016-05-20_clic_o2_v04_mu-_20GeV.root", "MyRecoMCTruthLinker_Forward/perftree", "p = 20 GeV", kOrange+7])



# parameters for pull plots

nbins_pulls = 50
range_pulls = 10. 


# parameters for residual plots

nbins_residuals = 100
# barrel
# range_residuals_omega = 0.000005 
# range_residuals_phi = 0.002 
# range_residuals_tanl = 0.002 
# range_residuals_d0 = 0.05 
# range_residuals_z0 = 0.05

# forward
range_residuals_omega = 0.00008 
range_residuals_phi = 0.01 
range_residuals_tanl = 0.005 
range_residuals_d0 = 0.1 
range_residuals_z0 = 0.1 

#parameters for resolution plots

logx = True
list_resolution_var = [] # name of var_y, name of var_x, min value of var_x, max value of var_x, npoints in x, log/linear scale for x
#list_resolution_var.append(["(recoPt-truePt)/(truePt*truePt)", "truePt", 0.1, 1000., 10, logx])
#list_resolution_var.append(["(recoPt-truePt)/(truePt*truePt)", "trueTheta", 10, 170., 32, not logx])

list_resolution_var.append(["pullOmega", "trueTheta", 7.5, 90, 16, not logx])
list_resolution_var.append(["pullPhi", "trueTheta", 7.5, 90, 16, not logx])
list_resolution_var.append(["pullTanLambda", "trueTheta", 7.5, 90, 16, not logx])
list_resolution_var.append(["pullD0", "trueTheta", 7.5, 90, 16, not logx])
list_resolution_var.append(["pullZ0", "trueTheta", 7.5, 90, 16, not logx])

# list_resolution_var.append(["pullOmega", "trueP", 2., 200., 6 , logx])
# list_resolution_var.append(["pullPhi", "trueP", 2., 200., 6 , logx])
# list_resolution_var.append(["pullTanLambda", "trueP", 1., 200., 6 , logx])
# list_resolution_var.append(["pullD0", "trueP", 2., 200., 6 , logx])
# list_resolution_var.append(["pullZ0", "trueP", 2., 200., 6 , logx])


#parameters for check plots

list_check_var = [] # name of var, nbins, min value, max value
list_check_var.append(["recoNhits", 20, 3, 23])
# list_check_var.append(["recoChi2OverNDF", 100, 0, 10])
# list_check_var.append(["trueTheta", 90, 0., 180])
# list_check_var.append(["pow(resOmega/pullOmega,2)", 100, 0, 2.e-9])
# list_check_var.append(["pow(resPhi/pullPhi,2)", 100, 0, 5.e-6])
# list_check_var.append(["pow(resTanLambda/pullTanLambda,2)", 100, 0, 5.e-6])
# list_check_var.append(["pow(resD0/pullD0,2)", 100, 0, 2.e-3])
# list_check_var.append(["pow(resZ0/pullZ0,2)", 100, 0, 5.e-3])

#parameter for effieciency plots

apply_theta_cut = True
apply_pt_cut = True
logx = True
list_efficiency_var = [] # name of var, min value of var, max value of var, npoints in x, log/linear scale for x
list_efficiency_var.append(["pt_reconstructable", 0.1, 1000., 25 , logx])
list_efficiency_var.append(["theta_reconstructable", 25, 155, 25 , not logx])

#parameter for fake rate plots

#apply_theta_cut = True
#logx = True
list_purity_var = [] # name of var, min value of var, max value of var, npoints in x, log/linear scale for x
#list_purity_var.append(["mc_p", 0.1, 1000., 25 , logx])
list_purity_var.append(["mc_p", 0.1, 400., 25 , logx])
#list_purity_var.append(["mc_theta", 0.436, 2.70, 25 , not logx])
list_purity_var.append(["mc_theta", 0.436, 1.57, 10 , not logx])


########################################

#gStyle.SetOptStat(111111)
#gStyle.SetOptFit(111111)
####gROOT.ProcessLine(".x ../macros/style/CLICdpSettingsScript.C")
gStyle.SetOptStat(0000)
gStyle.SetOptFit(0000)
#gStyle.SetOptFit(1)

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
  

def WriteTex(t):

    file_tex = open("tab_pulls_tex.tex", "w")
    file_txt = open("tab_pulls.txt", "w")

    numrows = len(t)    
    numcols = len(t[0])
    nrep = (numcols-1)*2
    nrep2 = (numcols-1)
    tab_format = "\\begin{tabular}{l"+ "r" * nrep + "} \n"
    tab_format_2 = ""
    for iname in range(numcols-1):
        title = t[0][iname+1]
        if title.find('#')>-1:
            title = title.replace('#', '\\')
            title = "$" + title + "$"
        tab_format_2 += "&  \multicolumn{2}{c}{"+ title + "} "
    tab_format_2 += "\\\  \n"
    tab_format_3 = "&  \multicolumn{1}{c}{$\mu_{fit}$} &  \multicolumn{1}{c}{$\sigma_{fit}$} " * nrep2  + " \\\  \n"
    
    file_tex.write("\documentclass[12pt]{article}\n")
    file_tex.write("\usepackage{rotating}\n")
    file_tex.write("\\begin{document}\n")
    file_tex.write("\\begin{table}\centering\n")
    file_tex.write("\pagenumbering{gobble}\n")
    file_tex.write("\\begin{turn}{90}\n")
    file_tex.write(tab_format)
    file_tex.write("\hline\hline \n")
    file_tex.write(tab_format_2)
    file_tex.write(tab_format_3)
    file_tex.write("\hline \n")
    for i in range(numrows-1):
        line = t[i+1][0] 
        file_txt.write('------------------------------------------------------ \n')
        file_txt.write('{:22s} {:15s}  {:15s} \n'.format(line, '     mean', '    sigma'))
        for j in range(numcols-1):
            mu_and_sigma = t[i+1][j+1].split('   ', 1 )
            file_txt.write('{:22s} {:15s}  {:15s} \n'.format(t[0][j+1], mu_and_sigma[0], mu_and_sigma[1]))
            mu_and_sigma[0] = mu_and_sigma[0].replace(" +- ", " \pm ")
            mu_and_sigma[1] = mu_and_sigma[1].replace(" +- ", " \pm ")
            line = line + " & $" + mu_and_sigma[0] + "$ & $" + mu_and_sigma[1] +"$"
        line = line + "\\\  \n"
        file_tex.write(line)
    file_tex.write("\hline \hline  \n")
    file_tex.write("\end{tabular}  \n") 
    file_tex.write("\end{turn}  \n") 
    file_tex.write("\end{table}  \n") 
    file_tex.write("\end{document}  \n")

    file_tex.close()
    file_txt.write('------------------------------------------------------ \n')
    file_txt.close()

    file_txt = open("tab_pulls.txt", "r")
    print file_txt.read()
    file_txt.close()

            
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
# 
    histo.SetLineColor(f[3])
    # histo.SetLineWidth(3)
    histo.SetMarkerColor(f[3])

    return histo
        


def DoGaussFit(f, v, cond=None, norm=False):

    if not isinstance(cond,str):
        cond = ""
        
    histo = GetHisto(f, v, cond)

    histo.Sumw2()
    
    if norm:
        histo.Scale(1./histo.Integral(0,v[1]+1))

    fitsf = TF1("fit","gaus",v[2],v[3])
    histo.Fit(fitsf,"N","",v[2],v[3])

    fitsf.SetLineColor(histo.GetLineColor())
    # fitsf.SetLineWidth(3)

    return histo ,fitsf



def GetRangeRMSfrac(histo, frac=0.90, check=None):

    
    # include under/overflows
    first_bin = 0
    last_bin = histo.GetNbinsX()+1

    integral_all = histo.Integral(first_bin, last_bin)
    integral_frac = integral_all*frac
    effective_frac = 0

    print "integral_all integral_frac = ", [integral_all, integral_frac]

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
                effective_frac = integral
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

    if check:
        effective_frac = float(effective_frac)/float(integral_all)
        return bin_min, bin_max, nbins_min, effective_frac
    
    return bin_min, bin_max
    
    
    
def GetResolution(f, v):

    # v array:
    #    0      1      2      3       4      5
    # [var_y, var_x, min_x, max_x, npoints, logx]

    c_resolution = TCanvas( 'c_resolution', 'resultion fits', 200, 10, 800, 600 )
    npad = math.sqrt(v[4])
    c_resolution.Divide(int(math.ceil(npad)), int(math.floor(npad)))
    hvec = []
    fvec = []

    
    graph = TGraphErrors()
    graph_mean = TGraphErrors()

    # step and cond computed differently according to log / not log x scale in order to always have equispaced bins 

    if v[5]:
        if v[2]==0:
            v[2] = 0.1
            print "ERROR: log scale does not allow 0 as min value - min now set to 0.1"
        step = math.pow(v[3]/v[2], 1./v[4])
    else:
        step = (v[3]-v[2])/v[4]
    
    
    lowedge_step = v[2]

    counter = 0
    for i in range(0,v[4]):

        if v[5]:
            upedge_step = lowedge_step*step
            cond = v[1] + " >= " + str(lowedge_step) + " && "  + v[1] + " < " + str(upedge_step)
            x = (lowedge_step + upedge_step)/2.
            x_err = (upedge_step - lowedge_step)/2.
            lowedge_step = upedge_step
        else:
            cond = v[1] + " >= " + str(v[2]+i*step) + " && "  + v[1] + " < " + str(v[2]+(i+1)*step)
            x = (v[2]+i*step + v[2]+(i+1)*step)/2.
            x_err = step/2.
        print "cond = ", cond

        # in two steps because very large outlayers can cause to have all the distribution in one or very few bins biasing the gaussian fit
        effective_frac = 1
        effective_nbins = -1
        y_min = 0
        y_max =0
        helper_histo_y = GetHisto(f, [v[0], effective_nbins, y_min, y_max], cond)
        #if helper_histo_y.GetEntries()>5:
        if helper_histo_y.GetEntries()>5:
            #while effective_nbins<3:
            while effective_nbins<20:
                # frac =  0.85 + ( 1 - effective_frac ) # effective frac not needed, in rms90 also under/overflow are counted
                # frac =  0.99
                
                #TEST
                y_min = helper_histo_y.GetMean()-helper_histo_y.GetRMS()
                y_max = helper_histo_y.GetMean()+helper_histo_y.GetRMS()
                helper_histo_y = GetHisto(f, [v[0], 100, y_min, y_max], cond)
                #effective_nbins = 300

                frac =  0.80
                print "helper_histo_y integral = ", helper_histo_y.Integral(0, helper_histo_y.GetNbinsX()+1)
                #RMS90
                min_bin, max_bin, effective_nbins, effective_frac = GetRangeRMSfrac(helper_histo_y, frac, True)
                #y_min = helper_histo_y.GetBinLowEdge(1)+(min_bin-1)*helper_histo_y.GetBinWidth(1)
                #y_max = helper_histo_y.GetBinLowEdge(1)+max_bin*helper_histo_y.GetBinWidth(1)
                y_min = helper_histo_y.GetBinLowEdge(1)+(min_bin-1)*helper_histo_y.GetBinWidth(1)
                y_max = helper_histo_y.GetBinLowEdge(1)+(max_bin+2)*helper_histo_y.GetBinWidth(1)
                helper_histo_y = GetHisto(f, [v[0], 300, y_min, y_max], cond)
             
                

            # frac =  0.05
            # print "helper_histo_y integral = ", helper_histo_y.Integral(0, helper_histo_y.GetNbinsX()+1)
            # min_bin, max_bin, effective_nbins, effective_frac = GetRangeRMSfrac(helper_histo_y, frac, True)
            # y_min = helper_histo_y.GetBinLowEdge(1)+(min_bin-1)*helper_histo_y.GetBinWidth(1)
            # y_max = helper_histo_y.GetBinLowEdge(1)+max_bin*helper_histo_y.GetBinWidth(1)
            # helper_histo_y = GetHisto(f, [v[0], 100, y_min, y_max], cond)


            if v[0].find('pull')>-1:
                y_min = -10.
                y_max = 10.

            nbins = 100
            histo_y, fit_y = DoGaussFit(f, [v[0], nbins, y_min, y_max], cond)

            c_resolution.cd(i+1)

            #if histo_y.GetEntries()>20:
            if histo_y.GetEntries()>500:
                hvec.append(histo_y)
                fvec.append(fit_y)
                hvec[-1].Draw("same")
                fvec[-1].Draw("same")

                graph.SetPoint(counter,x,fit_y.GetParameter(2))
                graph.SetPointError(counter,x_err,fit_y.GetParError(2))
                graph_mean.SetPoint(counter,x,fit_y.GetParameter(1))
                graph_mean.SetPointError(counter,x_err,fit_y.GetParError(1))
                counter += 1

            # print "y_min = ", y_min
            # print "y_max = ", y_max
           
        # else:
        #     graph.SetPoint(i,0,0)
        #     graph.SetPointError(i,0,0)

    GetResolution.ncalls += 1
    c_resolution.Print("many_plots_"+str(GetResolution.ncalls)+".eps")
    return graph, graph_mean
GetResolution.ncalls = 0


########################################


def getBinEdges(logx, nbins, min, max):

    step=0.
    
    if logx:
        if min==0:
            min=0.1
            print "ERROR: log scale does not allow 0 as min value - min now set to 0.1"
        step = math.pow(max/min, 1./nbins)
    else:
        step = float(max-min)/nbins   
    
    edge = min

    vec_edges = []
    vec_edges.append(edge)
    
    for i in range(0,nbins):

        if logx:
            edge = edge*step
        else:
            edge = edge+step

        vec_edges.append(edge)

    return vec_edges



            
def FillHisto(f, v, c, h):

    
    file = ROOT.TFile.Open(f[0], "read")
    tree = file.Get(f[1])
         
    tree.Project(h.GetName(),v[0],c)
    h = gDirectory.GetList().FindObject(h.GetName())

    print '-- h.GetName() = ', h.GetName()
    print '-- c = ', c
    print '-- h.GetEntries() = ', h.GetEntries()

    h.SetLineColor(f[3])
    h.SetMarkerColor(f[3])

    return h
        

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
        #leg = TLegend(0.15, 0.88-0.04*len(list_files), 0.25, 0.88)
        leg = TLegend(0.20, 0.88-0.04*len(list_files), 0.30, 0.88)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
        legs.append(leg)
        c_data_pulls.cd(list_var.index(var)+1)
        for file in list_files:
            h, f = DoGaussFit(file, var, "", True)            
            hists.append(h)
            fits.append(f)
            hists[-1].Draw("same")
            fits[-1].Draw("same")
            print "----- entries hists[-1] = ", hists[-1].GetEntries()
            legs[-1].AddEntry(hists[-1], file[2], "l")
            tab[0][list_files.index(file)+1] = file[2]
            tab[list_var.index(var)+1][0] = var[0]
            tab[list_var.index(var)+1][list_files.index(file)+1] = "%0.3f" % (fits[-1].GetParameter(1)) +" +- " + "%0.3f" % (fits[-1].GetParError(1)) +"   "+ "%0.3f" % (fits[-1].GetParameter(2)) +" +- " + "%0.3f" % (fits[-1].GetParError(2)) 
        legs[-1].Draw("same")
    c_data_pulls.Print("pulls.eps")

    WriteTex(tab)
    


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
        leg.SetTextSize(0.03)
        legs.append(leg)
        c_data_residuals.cd(list_var.index(var)+1)
        for file in list_files:            
            hists.append(GetHisto(file, var))
            hists[-1].Scale(1./hists[-1].Integral(0,var[1]+1))
            hists[-1].Draw("same")
            print "-- leg = ", file[2]
            print "-- var = ", var[0]
            print "-- rms = ", hists[-1].GetRMS()
            legs[-1].AddEntry(hists[-1], file[2], "l");
        legs[-1].Draw("same")

    c_data_residuals.Print("residuals.eps")


########################################
    
    
if draw_resolution:


    for var in list_resolution_var:

        graphs = []
        mg = TMultiGraph()

        graphs_mean = []
        mg_mean = TMultiGraph()
        
        leg = TLegend(0.70, 0.88-0.04*len(list_files), 0.88, 0.88)
        # leg = TLegend(0.70, 0.15+0.04*len(list_files), 0.88, 0.15)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
    
        for file in list_files:

            gr_sigma, gr_mean = GetResolution(file, var)
            graphs.append(gr_sigma)
            graphs[-1].SetMarkerSize(1.)
            graphs[-1].SetMarkerStyle(20)
            graphs[-1].SetLineColor(file[3])
            graphs[-1].SetLineWidth(3)
            graphs[-1].SetMarkerColor(file[3])
            leg.AddEntry(graphs[-1], file[2], "lp");
            mg.Add(graphs[-1])

            graphs_mean.append(gr_mean)
            graphs_mean[-1].SetMarkerSize(1.)
            graphs_mean[-1].SetMarkerStyle(20)
            graphs_mean[-1].SetLineColor(file[3])
            graphs_mean[-1].SetLineWidth(3)
            graphs_mean[-1].SetMarkerColor(file[3])
            mg_mean.Add(graphs_mean[-1])

        c_data_resolution = TCanvas( 'c_data_resolution', 'resolution', 200, 10, 800, 600 )
        #gPad.SetLogy()
        if var[5]:
            gPad.SetLogx()
        mg.Draw("AP")
        mg.GetXaxis().SetTitle(var[1]);
        mg.GetYaxis().SetTitle("#sigma("+var[0]+")");
        if var[0].find('pull')>-1:
            line = TLine(var[2], 1., var[3], 1.)
            line.SetLineStyle(7)
            line.SetLineWidth(2)
            line.SetLineColor(kBlack)
            line.Draw("same")
        leg.Draw("same")

        title = "resolution_" +str(list_resolution_var.index(var)+1)+".eps"
        c_data_resolution.Print(title)

        c_data_mean = TCanvas( 'c_data_mean', 'mean', 200, 10, 800, 600 )
        if var[5]:
            gPad.SetLogx()
        mg_mean.Draw("AP")
        mg_mean.GetXaxis().SetTitle(var[1]);
        mg_mean.GetYaxis().SetTitle("#mu("+var[0]+")");
        if var[0].find('pull')>-1:
            line_mean = TLine(var[2], 0., var[3], 0.)
            line_mean.SetLineStyle(7)
            line_mean.SetLineWidth(2)
            line_mean.SetLineColor(kBlack)
            line_mean.Draw("same")
        leg.Draw("same")

        title = "mean_" +str(list_resolution_var.index(var)+1)+".eps"
        c_data_mean.Print(title)

########################################

    
if draw_check_var:
    
    for var in list_check_var:
        c_data_check_var = TCanvas( 'c_data_check_var', 'check variables', 200, 10, 800, 600 )
        #gPad.SetLogx()
        hists = []
        legs = []        
        leg = TLegend(0.15, 0.88-0.04*len(list_files), 0.25, 0.88)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
        legs.append(leg)
        for file in list_files:
            hists.append(GetHisto(file, var))
            if hists[-1].Integral(0,var[1]+1)>0:
                hists[-1].Scale(1./hists[-1].Integral(0,var[1]+1)) #including outliers in integral for correct normalisation
            hists[-1].SetLineWidth(2)
            hists[-1].Draw("same")
            legs[-1].AddEntry(hists[-1], file[2], "l");
        legs[-1].Draw("same")
        hists[0].GetXaxis().SetTitle(var[0]);
        hists[0].GetYaxis().SetTitle("Events per u.a.");
        title = var[0]
        title = title.replace("(","_")
        title = title.replace(")","_")
        title = title.replace(",","_")
        title = title.replace(".","_")
        title = title.replace("/","_")
        # c_data_check_var.Print(var[0]+".eps")
        c_data_check_var.Print(title+".eps")


########################################
   
    
if draw_hit_residuals:

    list_subdet_layers = [[1, [0,1,2,3,4,5]] , [3, [0,1]] , [5, [0,1,2]] ]
    ntotlayers = 0
    for det in list_subdet_layers:
        ntotlayers += len(det[1])
        print "ntotlayers = ", ntotlayers
    
    list_var = []
    list_var.append(["resU", 50, -0.1, 0.1])
    list_var.append(["resV", 50, -1, 1])

    hists = []
    legs = []
    for var in list_var:

        c_data_hit_residuals = TCanvas( 'c_data_hit_residuals', 'hit_residuals', 200, 10, 800, 600 )
        npad = math.sqrt(ntotlayers)
        c_data_hit_residuals.Divide(int(math.ceil(npad)), int(math.floor(npad)))

        counter = 0
        for det in list_subdet_layers:
            for layer in det[1]:
                cond = "subdet == " + str(det[0]) + " && layer == " +str(layer)
                print "--- cond = ", cond
                counter += 1
                c_data_hit_residuals.cd(counter)

                leg = TLegend(0.15, 0.88-0.04*len(list_files), 0.25, 0.88)
                leg.SetBorderSize(0)
                leg.SetTextSize(0.03)
                legs.append(leg)
                legs[-1].AddEntry("", cond, "");
                        
                for file in list_files:
                    factor = 1.
                    if det[0] == 1:
                        factor = 100.
                        if var[0] == "resV":
                            factor = 500.
                    elif det[0] == 3:
                        factor = 5.
                    hists.append(GetHisto(file, [var[0], var[1], var[2]/factor, var[3]/factor], cond))
                    hists[-1].Draw("same")
                    legs[-1].AddEntry(hists[-1], file[2], "l");
                legs[-1].Draw("same")

        c_data_hit_residuals.Print("hit_"+var[0]+".eps")


########################################
    
    
if draw_efficiency:

    for var in list_efficiency_var:

        # arguments of getBinEdges: logx, nbins, min, max
        var_edges = getBinEdges(var[4],var[3],var[1],var[2])
        print ' var_edges = ', var_edges

        hists = []
        graphs = []



        leg = TLegend(0.70, 0.22+0.04*len(list_files), 0.88, 0.22)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
            
           
        c_data_efficiency = TCanvas( 'c_data_efficiency', 'efficiency', 200, 10, 800, 600 )


        
        for file in list_files:


            name_hist_reconstructed = 'h_var_reconstructed_'+str(list_efficiency_var.index(var))
            name_hist_reconstructable = 'h_var_reconstructable_'+str(list_efficiency_var.index(var))
            h_var_reconstructed = TH1F(name_hist_reconstructed, name_hist_reconstructed, var[3] , array('d',var_edges))
            h_var_reconstructable = TH1F(name_hist_reconstructable, name_hist_reconstructable, var[3] , array('d',var_edges))
            h_var_reconstructed.SetDirectory(0)
            h_var_reconstructable.SetDirectory(0)


            #loop on tree

            f = ROOT.TFile.Open(file[0], "read")
            mytree = f.Get(file[1])
            
            for event in mytree:
                index = 0
                vector_variable = getattr(mytree, var[0])
                for element in vector_variable:
                    if apply_theta_cut:
                        pass_theta = math.fabs(math.cos(event.theta_reconstructable[index]*math.pi/180.))<0.89
                    else:
                        pass_theta = True
                    if apply_pt_cut:
                        #print '--pt = ', event.pt_reconstructable[index]
                        pass_pt = event.pt_reconstructable[index] >= 0.1   
                    else:
                        pass_pt = True
                    print '--pass_theta = ', pass_theta
                    if (pass_theta and pass_pt):
                        h_var_reconstructable.Fill( element )
                        if event.is_reconstructed[index]:
                            h_var_reconstructed.Fill( element )
                    index=index+1
            

            h_var_reconstructable.SetLineColor(kOrange+7)
            h_var_reconstructed.SetLineColor(kAzure+7)
                    
            hists.append(h_var_reconstructed)
            hists.append(h_var_reconstructable)


            g_eff = TGraphAsymmErrors()
            g_eff.Divide( h_var_reconstructed, h_var_reconstructable, "v" )
            g_eff.SetMarkerColor(file[3])
            g_eff.SetLineColor(file[3])

            graphs.append(g_eff)
            leg.AddEntry(g_eff, file[2], "lp");


            if var[4]:
                gPad.SetLogx()

            if list_files.index(file)==0:                    
                g_eff.GetYaxis().SetTitle("#epsilon_{trk}")
                if var[0].find('theta')>-1:
                    g_eff.GetXaxis().SetTitle("#theta [rad]")
                elif var[0].find('pt')>-1:
                    g_eff.GetXaxis().SetTitle("p_{T} [GeV]")
                else:
                    g_eff.GetXaxis().SetTitle(var[0])
                g_eff.Draw('AP')

            else:
                g_eff.Draw('P')

                
        leg.Draw('same')
        c_data_efficiency.Print("efficiency_"+var[0]+".eps")

            
########################################
    
    
if draw_purity:

    for var in list_purity_var:

        # arguments of getBinEdges: logx, nbins, min, max
        var_edges = getBinEdges(var[4],var[3],var[1],var[2])
        print ' var_edges = ', var_edges

        hists = []
        graphs = []


        leg = TLegend(0.70, 0.22+0.04*len(list_files), 0.88, 0.22)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
           
        c_data_purity = TCanvas( 'c_data_purity', 'fake rate', 200, 10, 800, 600 )

        
        for file in list_files:


            name_hist_fake = 'h_var_fake_'+str(list_purity_var.index(var))
            name_hist_all = 'h_var_all_'+str(list_purity_var.index(var))
            h_var_fake = TH1F(name_hist_fake, name_hist_fake, var[3] , array('d',var_edges))
            h_var_all = TH1F(name_hist_all, name_hist_all, var[3] , array('d',var_edges))


            #loop on tree

            f = ROOT.TFile.Open(file[0], "read")
            mytree = f.Get(file[1])

            den = 0
            num = 0
            
            for event in mytree:
                index = 0
                vector_variable = getattr(mytree, var[0])
                for element in vector_variable:
                    if apply_theta_cut:
                        pass_theta = math.fabs(math.cos(event.mc_theta[index]))<0.89
                    else:
                        pass_theta = True
                    if apply_pt_cut:
                        #print '--pt = ', event.pt_reconstructable[index]
                        pass_pt = ( math.fabs(event.mc_p[index]*math.sin(event.mc_theta[index])) ) > 0.1   
                    else:
                        pass_pt = True
                    #pass_pt = ( math.fabs(event.mc_p[index]*math.sin(event.mc_theta[index])) ) > 1
                    #bad_track = event.trk_nhits_vtx[index]==4 and event.trk_nhits_trk[index]==0
                    #pt = math.fabs(event.mc_p[index]*math.sin(event.mc_theta[index]))
                    if (pass_theta and pass_pt):
                        h_var_all.Fill( element )
                        den = den+1
                        #h_var_all.Fill( element*math.fabs(math.sin(event.mc_theta[index])) )
                        if event.trk_purity[index]<0.9:
                            num = num+1
                            h_var_fake.Fill( element )
                            #h_var_fake.Fill( element*math.fabs(math.sin(event.mc_theta[index])) )
                    index=index+1
            
            hists.append(h_var_fake)
            hists.append(h_var_all)

            print "----- FAKE RATE num = ", num
            print "----- FAKE RATE den = ", den
            print "----- FAKE RATE = ", num/den

            g_fake = TGraphAsymmErrors()
            g_fake.Divide( h_var_fake, h_var_all, "v" )
            g_fake.SetMarkerColor(file[3])
            g_fake.SetLineColor(file[3])


            graphs.append(g_fake)

            leg.AddEntry(g_fake, file[2], "lp");

            gPad.SetLogy()
            if var[4]:
                gPad.SetLogx()

            if list_files.index(file)==0:
                g_fake.GetYaxis().SetTitle("fake rate")
                if var[0].find('theta')>-1:
                    g_fake.GetXaxis().SetTitle("#theta [rad]")
                elif var[0].find('pt')>-1:
                    g_fake.GetXaxis().SetTitle("p_{T} [GeV]")
                else:
                    g_fake.GetXaxis().SetTitle(var[0])
                g_fake.GetYaxis().SetRangeUser(1.e-5,1.)
                g_fake.Draw('AP')
            else:
                g_fake.Draw('P')

        leg.Draw('same')

        c_data_purity.Print("fake_rate_"+var[0]+".eps")

            
########################################
