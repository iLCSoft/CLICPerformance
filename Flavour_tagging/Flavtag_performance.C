#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLegend.h>
#include <TUnixSystem.h>
#include <TCut.h>
#include <TROOT.h>


void graph(int numax, Double_t* x, Double_t* y, Double_t* ex, Double_t* ey, TString style = "ACP", int col=1, TString name = "LCFIPlus", bool beauty = true, bool bckg = false)
{
    
    TGraphErrors* gr = new TGraphErrors(numax,x,y,ex,ey);
    gr->SetName(name);
    gr->SetTitle("");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0);
    gr->SetMarkerColor(col);
    gr->SetLineColor(col);
    gr->SetLineWidth(3);
    if (bckg) gr->SetLineStyle(2);
    
    if (beauty)
    {
        gr->GetXaxis()->SetLimits(0.39,1.005);
        gr->SetMaximum(1.05);
        gr->SetMinimum(8.5e-4);
        gr->GetYaxis()->SetTitle("Misidentification eff.");
        gr->GetYaxis()->SetTitleSize(0.07);
        gr->GetXaxis()->SetTitleSize(0.07);
        gr->GetYaxis()->SetLabelSize(0.05);
        gr->GetXaxis()->SetLabelSize(0.05);
        gr->GetXaxis()->SetTitle("Beauty eff.");

    }
    else
    {
        gr->GetXaxis()->SetLimits(0.,1.005);
        gr->SetMaximum(1.05);
        gr->SetMinimum(2.0e-3);
        gr->GetYaxis()->SetTitle("Misidentification eff.");
        gr->GetYaxis()->SetTitleSize(0.07);
        gr->GetYaxis()->SetLabelSize(0.05);
        gr->GetXaxis()->SetTitleSize(0.07);
        gr->GetXaxis()->SetLabelSize(0.05);
        gr->GetXaxis()->SetTitle("Charm eff.");

    }
    
    gr->Draw(style);
}

void Flavtag_performance( bool train=false , TString flav = "B" ) {

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetMarkerSize(0.7);
    gROOT->ForceStyle();

    bool B = true;
    if (flav == "C") B = false;
    TString dirName[2];
    dirName[0] = "lcfiweights_flavtag_500gev_20to90deg_CT_vlc12_3TeV_tight";
    dirName[1] = "lcfiweights_flavtag_500gev_20to90deg_fastjet_vlc12_08_conformal_new_new_no_overlay";

    TCanvas* c1 = new TCanvas("c1","canvas1",600,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1.0);

    TLegend* leg1 = new TLegend();
    if (B) leg1 = new TLegend(0.5,0.1,0.72,0.8,"", "PL");
    if (!B) leg1 = new TLegend(0.4,0.22,0.67,0.92,"", "PL");

    const int NMAX=50;
   
    double xb[NMAX];
    double xc[NMAX];
    double xo[NMAX];
    double exb[NMAX];
    double exc[NMAX];
    double exo[NMAX];
    
    for (unsigned int k=0; k<2; k++)
    {
        cout << dirName[k] << endl;
        void* dir = gSystem->OpenDirectory( dirName[k].Data() );

        vector<TString> files;
        files.push_back( "flavwgts_2jets_c0.root" );
        files.push_back( "flavwgts_2jets_c1.root" );
        files.push_back( "flavwgts_2jets_c2.root" );
        files.push_back( "flavwgts_2jets_c3.root" );

        cout << "Total of " << files.size() << " files in folder " << k+1 << endl;

        vector< TTree* > trees;
        vector< vector< TH1D* > > histSets;
        
        for (unsigned int i=0; i<files.size(); ++i) {
            TFile* f = new TFile(TString::Format("%s/%s",dirName[k].Data(),files[i].Data()));
            TDirectoryFile *lcfiplus = (TDirectoryFile*)f->Get("lcfiplus_dataset");
            TTree* ntp = (TTree*) lcfiplus->Get( train ? "TrainTree" : "TestTree" );
            trees.push_back(ntp);
        }
        
        double ntot1(0);
        double ntot2(0);
        double ntot3(0);
        for (unsigned int j=0; j<trees.size(); ++j)
        {
            ntot1 += trees[j]->GetEntries("classID==0");
            ntot2 += trees[j]->GetEntries("classID==1");
            ntot3 += trees[j]->GetEntries("classID==2");
        }
        
        
        for (int i=0; i<NMAX; ++i)
        {
            TCut cut = "jet";
            if (B) cut = TString::Format("jetB>%f",(double)i/(double)NMAX).Data();
            else cut = TString::Format("jetC>%f",(double)i/(double)NMAX).Data();
            double n1(0);
            double n2(0);
            double n3(0);
            for (unsigned int j=0; j<trees.size(); ++j)
            {
                n1 += trees[j]->GetEntries(cut+"classID==0");
                n2 += trees[j]->GetEntries(cut+"classID==1");
                n3 += trees[j]->GetEntries(cut+"classID==2");
            }
            xb[i] = 0;
            xc[i] = 0;
            xo[i] = 0;
            xb[i] = n1/ntot1;
            exb[i] = 0;//sqrt(xb[i]*(1-xb[i])/ntot1);
            xc[i] = n2/ntot2;
            exc[i] = 0;//sqrt(xc[i]*(1-xc[i])/ntot2);
            xo[i] = n3/ntot3;
            exo[i] = 0;//sqrt(xo[i]*(1-xo[i])/ntot3);
        }
        
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        c1->cd();
        
        if (B)
        {
 
            if (k==0)
            {
                graph(NMAX,xb,xc,exb,exc,"ACP",58,"c-jets: with background",B,true);
                graph(NMAX,xb,xo,exb,exo,"CP",80,"light jets: with background",B,true);
            }
            else if (k==1)
            {
                graph(NMAX,xb,xc,exb,exc,"CP",58,"c-jets: no background",B,false);
                graph(NMAX,xb,xo,exb,exo,"CP",80,"light jets: no background",B,false);
            }
        }
        else
        {

            if (k==0)
            {
                graph(NMAX,xc,xb,exc,exb,"ACP",2,"b-jets: with background",B,true);
                graph(NMAX,xc,xo,exc,exo,"CP",80,"light jets: with background",B,true);
            }
            else if (k==1)
            {
                graph(NMAX,xc,xb,exc,exb,"CP",2,"b-jets: no background",B,false);
                graph(NMAX,xc,xo,exc,exo,"CP",80,"light jets: no background",B,false);
            }
        }
        
    }
    c1->SetLogy(1);
    c1->BuildLegend(0.2,0.62,0.55,0.85,"", "PL");
    
    TPaveText *pt = new TPaveText(0.4,0.38,0.7,1.12,"br");
    if (!B)
    {
        c1->BuildLegend(0.2,0.62,0.55,0.85,"", "PL");
        pt = new TPaveText(0.05,0.48,0.45,1.12,"br");
    }
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetLineColor(0);
    pt->SetLineWidth(0);
    
    TText *pt_LaTex = pt->AddText("CLICdp work in progress");
    pt->Draw();
    pad1->Modified();
    c1->cd();
    
    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
    c1->ToggleToolBar();
    
    c1->SaveAs("flavtag_eff.pdf");
    
    
}
