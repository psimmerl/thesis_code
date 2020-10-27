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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TGraphErrors.h"
#include <map>

std::vector<double> z, x, y, errorz;

void monElastic2D( const char *fin, int run, int rebiny = 1) //,  int sect)
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
  std::string fpdf_in("mon_elastic2d_" + std::to_string(run) +".pdf[");
  c1->Print(fpdf_in.c_str());
  
  for( int ss = 0; ss < 6; ss++ ){

    c1->Clear();
    c1->Divide(3,4);
    std::string hist_name("hwtheta_el_pass_all_"+std::to_string(ss)+"_r"+std::to_string(run));

    TH2F* h2_testing = (TH2F*)ff->Get(hist_name.c_str());
    //    h2_testing->RebinY(rebiny);

    std::vector<double> x_center;
    std::vector<double> y_val;
    std::vector<double> y_err;
    std::vector<double> zeros;
    
    ofstream myfile;
    myfile.open ("mon_incl_elastic_sigfc_vs_theta_r"+std::to_string(run)+"_s"+std::to_string(ss)+".txt");

    for( int bb = 1; bb < h2_testing->GetXaxis()->GetNbins(); bb++ ){
      std::string h_testing_name = hist_name + "_" + std::to_string(bb);
      TH1D* h_testing = h2_testing->ProjectionY(h_testing_name.c_str(), bb, bb+1);
      h_testing->RebinX(rebiny); 
      double bin_size = h_testing->GetXaxis()->GetBinCenter(2) -h_testing->GetXaxis()->GetBinCenter(1);  
      std::cout << " bin size " << bin_size << std::endl;
   
      std::cout << " n bins " << h_testing->GetXaxis()->GetNbins() << std::endl;
      std::cout<< " h_testing entries " << h_testing->GetEntries() << std::endl;
      std::cout << " ss " << ss << std::endl;
      c1->cd(bb);
    
      ROOT::Fit::DataOptions opt;
      ROOT::Fit::DataRange range;
      
      ROOT::Fit::BinData data(opt, range );// &x[0], &y[0] );
      ROOT::Fit::FillData(data,h_testing);
      
      std::string f1_name = "f1_" + std::to_string(ss) + "_bb" + std::to_string(bb);
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
      //fitter.LikelihoodFit(data);
      fitter.Fit(data);
      
      const ROOT::Fit::FitResult & res = fitter.Result();
      res.Print(std::cout);
    // copy all fit result info (values, chi2, etc..) in TF1
      f1->SetFitResult(res);
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // draw the results
      std::string fit_name = "fit_check"+std::to_string(ss) + "_tb" + std::to_string(bb);
      std::string fit_name_bg = "fit_bg_check"+std::to_string(ss) + "_tb" + std::to_string(bb);
      std::string fit_name_signal = "fit_signal_check"+std::to_string(ss) + "_tb" + std::to_string(bb);
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
      
      std::string title("ep #rightarrow eX, El. (FD). S" + std::to_string(ss) + "#Theta Bin " + std::to_string(bb));
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
      double x_bin_center = h2_testing->GetXaxis()->GetBinCenter(bb);

      double integral_min = fit_check->GetParameter(1) - 3.0*fit_check->GetParameter(2);
      double integral_max = fit_check->GetParameter(1) + 3.0*fit_check->GetParameter(2);
      double integrated_events = fit_check->Integral(integral_min, integral_max)/bin_size;
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
      out_latex="S/FC " + std::to_string(signal_events/acc_beam_charge[run]);
      lab->DrawLatex(0.15, 0.65, out_latex.c_str());
          
      std::cout << " Creating Graph with Data " << std::endl;
      std::cout << " sector " << ss << " x " << x_bin_center << " mu " << mu << " sigma " << sigma << " sig/fc " << signal_events/acc_beam_charge[run] << std::endl; 
      if( integrated_events == 0 ) continue;
      x_center.push_back(x_bin_center);
      zeros.push_back(0.0);
      y_val.push_back(signal_events/acc_beam_charge[run]);
      y_err.push_back(0);//sqrt(signal_events)/acc_beam_charge[run]);
      myfile << x_bin_center << " " << sig_fc_ratio << std::endl;


    }
    myfile.close();

    std::string fpdf_current="mon_elastic2d_" + std::to_string(run) +".pdf";    
    c1->Print(fpdf_current.c_str());

    c1->Clear();
    c1->Divide(1,1);
    c1->cd(1);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //draw graphs now
    TGraphErrors *g_rates = new TGraphErrors(x_center.size(), &x_center[0], &y_val[0], &zeros[0], &y_err[0]);
    std::string g_title_name = "g_"+hist_name;
    g_rates->SetTitle(g_title_name.c_str());
    g_rates->SetMarkerColor(kRed);
    g_rates->SetMarkerSize(0.8);
    g_rates->SetMarkerStyle(8);
    g_rates->GetXaxis()->SetLimits(3,15);
    c1->SetLogy();
    g_rates->SetMinimum(0.0);
    g_rates->SetMaximum(0.0006);
    g_rates->Draw("AP");
    c1->Print(fpdf_current.c_str());
    
  }

  //close file
  std::string fpdf_out = "mon_elastic2d_" + std::to_string(run) +".pdf]";    
  c1->Print(fpdf_out.c_str());

}
