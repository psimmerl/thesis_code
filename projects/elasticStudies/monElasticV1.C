#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"
#include "TLatex.h"

// For defining the functions
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "TProfile.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLine.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include <map>

std::vector<double> z, x, y, errorz;

void monElasticV1( const char *fin, int run, int rebinx = 1) //,  int sect)
{
  std::string fname(fin);
  fname+=std::to_string(run);

  std::string can_name = "can";
  double min_x = 0.7; // was 0.9 
  double max_x = 1.3;

  std::map<int, double> acc_beam_charge;
  acc_beam_charge[5306] =  417438.34;
  acc_beam_charge[5418] =  51939.55;
  acc_beam_charge[5419] =41387.21;
  acc_beam_charge[5443]= 56469.41;
  acc_beam_charge[5444] = 118864.08;
  acc_beam_charge[5453] = 236320.25;
  

  TLatex *lab =new  TLatex();
  lab->SetNDC();
  lab->SetTextFont(42);
  lab->SetTextSize(0.05);
  lab->SetTextColor(1);

  TFile *ff = new TFile(fin);

  TCanvas *c1 = new TCanvas(can_name.c_str(), can_name.c_str(), 900, 900);
  c1->Divide(3,2);
  std::string fpdf_in("mon_elastic_" + std::to_string(run) +".pdf[");
  c1->Print(fpdf_in.c_str());
  
  ofstream myfile;
  myfile.open ("mon_incl_elastic_r"+std::to_string(run)+"_per_sector.txt");


  for( int ss = 0; ss < 6; ss++ ){
    std::string hist_name("hw_el_pass_all_"+std::to_string(ss)+"_r"+std::to_string(run));
    TH1F* h_testing = (TH1F*)ff->Get(hist_name.c_str());
    h_testing->RebinX(rebinx);
    std::cout<< " h_testing entries " << h_testing->GetEntries() << std::endl;
    std::cout << " ss " << ss << std::endl;
    c1->cd(ss+1);
    double bin_size = h_testing->GetXaxis()->GetBinCenter(2) -h_testing->GetXaxis()->GetBinCenter(1);
    
    std::cout << " bin size " << bin_size << std::endl;
   
    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange range;
  
    ROOT::Fit::BinData data(opt, range );// &x[0], &y[0] );
    ROOT::Fit::FillData(data,h_testing);
  
    std::string f1_name = "f1_" + std::to_string(ss);
    TF1 *f1 = new TF1(f1_name.c_str(),"gaus(0) + [3] + [4]*x +  [5]*x*x");

    ROOT::Math::WrappedMultiTF1 fitFunction(*f1, f1->GetNdim());
    
    //double par0[6] = { 300, 1, 0.1, -450, 650, -80};
    ROOT::Fit::Fitter fitter;
    fitter.SetFunction( fitFunction, false);
    //fitter.Config().SetParamsSettings(6, par0);
    fitter.Config().ParSettings(0).SetLimits(0,10000);
    fitter.Config().ParSettings(1).SetLimits(0.5, 1.5);
    fitter.Config().ParSettings(2).SetLimits(0,0.7);
    
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2");
    fitter.LikelihoodFit(data);
      //Fit(data);

    const ROOT::Fit::FitResult & res = fitter.Result();
    res.Print(std::cout);
    // copy all fit result info (values, chi2, etc..) in TF1
    f1->SetFitResult(res);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // draw the results
    std::string fit_name = "fit_check"+std::to_string(ss);
    std::string fit_name_bg = "fit_bg_check"+std::to_string(ss);
    std::string fit_name_signal = "fit_signal_check"+std::to_string(ss);
    TF1 *fit_check = new TF1(fit_name.c_str(),"gaus(0) + [3] + [4]*x + [5]*x*x",min_x, max_x);
    TF1 *fit_bg = new TF1(fit_name_bg.c_str(),"[0] + [1]*x + [2]*x*x",min_x, max_x);
    TF1 *fit_signal = new TF1(fit_name_signal.c_str(),"gaus(0)",min_x, max_x);
    fit_check->SetFitResult(res);

    fit_signal->SetParameter(0,fit_check->GetParameter(0));
    fit_signal->SetParameter(1,fit_check->GetParameter(1));
    fit_signal->SetParameter(2,fit_check->GetParameter(2));
    
    fit_bg->SetParameter(0,fit_check->GetParameter(3));
    fit_bg->SetParameter(1,fit_check->GetParameter(4));
    fit_bg->SetParameter(2,fit_check->GetParameter(5));

    fit_check->SetLineColor(kRed);
    fit_bg->SetLineColor(kGreen);
    fit_signal->SetLineColor(kBlue);
    h_testing->SetTitle("");
    h_testing->Draw("colz");
    fit_check->Draw("same");
    fit_bg->Draw("same");
    fit_signal->Draw("same");
    
    std::string title("ep #rightarrow eX, El. (FD). S" + std::to_string(ss));
    std::string xtitle("W (Gev)");
    std::string ytitle("counts");
    lab->DrawLatex(0.1, 0.925, title.c_str());
    lab->DrawLatex(0.5, 0.015, xtitle.c_str());
    lab->SetTextAngle(90);
    lab->DrawLatex(0.04, 0.5, ytitle.c_str());
    lab->SetTextAngle(0);
    
    TLine *l0 = new TLine(0.938, 100 , 0.938, 100);
    l0->SetLineColor(kRed);
    l0->SetLineWidth(2);
    l0->Draw("same");
    
    double mu = fit_check->GetParameter(1);
    double sigma = fit_check->GetParameter(2);
    
    double integral_min = fit_check->GetParameter(1) - 3.0*fit_check->GetParameter(2);
    double integral_max = fit_check->GetParameter(1) + 3.0*fit_check->GetParameter(2);
    double integrated_events = fit_check->Integral(integral_min, integral_max);
    double background_events = fit_bg->Integral(integral_min, integral_max)/bin_size;
    double signal_events = fit_signal->Integral(integral_min, integral_max)/bin_size;
    double sig_bg_ratio  = signal_events / background_events;
    double sig_fc_ratio = signal_events/acc_beam_charge[run];
    
    std::string out_latex("#mu = " + std::to_string(mu) + " #sigma = " + std::to_string(sigma));
    lab->DrawLatex(0.15, 0.84, out_latex.c_str());
    out_latex="Signal " + std::to_string(signal_events);
    lab->DrawLatex(0.15, 0.74, out_latex.c_str());
    out_latex="BG " + std::to_string(background_events);
    lab->DrawLatex(0.15, 0.70, out_latex.c_str());    
    TString text;
    text = Form("S/FC %0.6f", sig_fc_ratio);
    lab->DrawLatex(0.15, 0.65, text.Data());
    
    myfile << ss << " " << sig_fc_ratio << std::endl;
            
    
  }
  myfile.close();
  std::string fpdf_current="mon_elastic_" + std::to_string(run) +".pdf";    
  c1->Print(fpdf_current.c_str());

  std::string fpdf_out = "mon_elastic_" + std::to_string(run) +".pdf]";    
  c1->Print(fpdf_out.c_str());

}
