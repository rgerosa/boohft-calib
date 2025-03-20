import os
import sys
import glob
import argparse
import subprocess
import shutil
import time
import ROOT
import numpy as np
from array import array 

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-dir', type=str, default='', help='input directory with files');
parser.add_argument('-c', '--cat-id', type=str, default='', choices=['cat0','cat1','cat2','cat3','cat4'], help='category identifier');
parser.add_argument('-m', '--mass-obs', type=str, default='mSD', help='mass observable to fit');
parser.add_argument('-r', '--rebin-factor', type=int, default=10, help='rebin factor for plots');
parser.add_argument('-o', '--output-directory', type=str, default='', help='name of the output directory');
args = parser.parse_args()

ROOT.gInterpreter.ProcessLine('#include "CMS_style.h"')
ROOT.setTDRStyle();
ROOT.gStyle.SetOptStat(0);

ROOT.RooMsgService.instance().setSilentMode(True);
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR) ;

def getChi2(h,d):
    chi2 = 0;
    ndf  = 0;
    for i in range(0,h.GetNbinsX()):
        if h.GetBinContent(i+1) == 0 : continue;
        res = h.GetBinContent(i+1)-d.GetBinContent(i+1);
        chi2 += (res**2)/(h.GetBinError(i+1)**2);
        ndf = ndf + 1;
    return chi2,ndf;

def getMaxAndHWHM(h):
    binMax = h.GetMaximumBin();
    xMax = h.GetXaxis().GetBinCenter(binMax);
    yMax = h.GetBinContent(binMax);    
    halfMax = yMax / 2.0;
    binLeft = binMax;
    while binLeft > 1 and h.GetBinContent(binLeft) > halfMax:
        binLeft = binLeft-1;
    binRight = binMax;
    while binRight > 1 and h.GetBinContent(binRight) > halfMax:
        binRight = binRight+1;

    xLeft=h.GetXaxis().GetBinCenter(binLeft);
    xRight=h.GetXaxis().GetBinCenter(binRight);
    hwhm = (xRight-xLeft)/2.

    return xMax,hwhm

f_zjet = ROOT.TFile(args.input_dir+"/Zto2Q.root","READ");
f_wjet = ROOT.TFile(args.input_dir+"/Wto2Q.root","READ");
f_qcd = ROOT.TFile(args.input_dir+"/QCD.root","READ");
f_data = ROOT.TFile(args.input_dir+"/JetMET_2022EE.root","READ");

h_zjet_pass = f_zjet.Get("hist_ak8_"+args.mass_obs);
h_wjet_pass = f_wjet.Get("hist_ak8_"+args.mass_obs);
h_qcd_pass = f_qcd.Get("hist_ak8_"+args.mass_obs);
h_data_pass = f_data.Get("hist_ak8_"+args.mass_obs);

h_zjet_fail = f_zjet.Get("hist_ak8_"+args.mass_obs+"_f");
h_wjet_fail = f_wjet.Get("hist_ak8_"+args.mass_obs+"_f");
h_qcd_fail = f_qcd.Get("hist_ak8_"+args.mass_obs+"_f");
h_data_fail = f_data.Get("hist_ak8_"+args.mass_obs+"_f");

obs = ROOT.RooRealVar(args.mass_obs,"",h_data_pass.GetXaxis().GetBinCenter(h_data_pass.GetMaximumBin()),h_data_pass.GetXaxis().GetXmin(),h_data_pass.GetXaxis().GetXmax())
obs.setBins(h_data_pass.GetNbinsX());
obs_list = ROOT.RooArgList();
obs_list.add(obs);

w_pass = ROOT.RooWorkspace("w_pass","");
w_fail = ROOT.RooWorkspace("w_fail","");

## data histogram conversion
rh_data_pass = ROOT.RooDataHist("data_pass","",obs_list,h_data_pass)
rh_data_fail = ROOT.RooDataHist("data_fail","",obs_list,h_data_fail)
w_pass.Import(rh_data_pass);
w_fail.Import(rh_data_fail);

######################
## QCD MC modelling ##
######################

rh_qcd_pass = ROOT.RooDataHist("qcd_pass","",obs_list,h_qcd_pass)
rh_qcd_fail = ROOT.RooDataHist("qcd_fail","",obs_list,h_qcd_fail)

## order tune in MC
qcd_coef_pass_1 = ROOT.RooRealVar("qcd_coef_pass_1","qcd_coef_pass_1",0.001,-10,10)
qcd_coef_pass_2 = ROOT.RooRealVar("qcd_coef_pass_2","qcd_coef_pass_2",0.001,-10,10)
qcd_coef_pass_3 = ROOT.RooRealVar("qcd_coef_pass_3","qcd_coef_pass_3",0.001,-10,10)
qcd_coef_pass_4 = ROOT.RooRealVar("qcd_coef_pass_4","qcd_coef_pass_4",0.001,-10,10)
qcd_coef_pass_5 = ROOT.RooRealVar("qcd_coef_pass_5","qcd_coef_pass_5",0.001,-10,10)

qcd_coef_fail_1 = ROOT.RooRealVar("qcd_coef_fail_1","qcd_coef_fail_1",0.001,-10,10)
qcd_coef_fail_2 = ROOT.RooRealVar("qcd_coef_fail_2","qcd_coef_fail_2",0.001,-10,10)
qcd_coef_fail_3 = ROOT.RooRealVar("qcd_coef_fail_3","qcd_coef_fail_3",0.001,-10,10)
qcd_coef_fail_4 = ROOT.RooRealVar("qcd_coef_fail_4","qcd_coef_fail_4",0.001,-10,10)
qcd_coef_fail_5 = ROOT.RooRealVar("qcd_coef_fail_5","qcd_coef_fail_5",0.001,-10,10)

if args.cat_id == "cat0":
    pdf_qcd_pass = ROOT.RooChebychev("pdf_qcd_pass","",obs,[qcd_coef_pass_1,qcd_coef_pass_2,qcd_coef_pass_3]);
    pdf_qcd_fail = ROOT.RooChebychev("pdf_qcd_fail","",obs,[qcd_coef_fail_1,qcd_coef_fail_2,qcd_coef_fail_3,qcd_coef_fail_4,qcd_coef_fail_5]);
elif args.cat_id == "cat1" or args.cat_id == "cat2":
    pdf_qcd_pass = ROOT.RooChebychev("pdf_qcd_pass","",obs,[qcd_coef_pass_1,qcd_coef_pass_2,qcd_coef_pass_3,qcd_coef_pass_4]);
    pdf_qcd_fail = ROOT.RooChebychev("pdf_qcd_fail","",obs,[qcd_coef_fail_1,qcd_coef_fail_2,qcd_coef_fail_3,qcd_coef_fail_4,qcd_coef_fail_5]);
else:
    pdf_qcd_pass = ROOT.RooChebychev("pdf_qcd_pass","",obs,[qcd_coef_pass_1,qcd_coef_pass_2,qcd_coef_pass_3,qcd_coef_pass_4,qcd_coef_pass_5]);
    pdf_qcd_fail = ROOT.RooChebychev("pdf_qcd_fail","",obs,[qcd_coef_fail_1,qcd_coef_fail_2,qcd_coef_fail_3,qcd_coef_fail_4,qcd_coef_fail_5]);
    
fit_qcd_pass_res = pdf_qcd_pass.fitTo(rh_qcd_pass,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));
fit_qcd_fail_res = pdf_qcd_fail.fitTo(rh_qcd_fail,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));

pdf_qcd_pass_norm = ROOT.RooRealVar(pdf_qcd_pass.GetName()+"_norm","",rh_qcd_pass.sumEntries())
pdf_qcd_fail_norm = ROOT.RooRealVar(pdf_qcd_fail.GetName()+"_norm","",rh_qcd_fail.sumEntries())
pdf_qcd_pass_norm.setConstant(False);
pdf_qcd_fail_norm.setConstant(False);

h_fit_qcd_pass = pdf_qcd_pass.createHistogram("h_fit_qcd_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_qcd_fail = pdf_qcd_fail.createHistogram("h_fit_qcd_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_qcd_pass.Scale(pdf_qcd_pass_norm.getVal()*args.rebin_factor);
h_fit_qcd_fail.Scale(pdf_qcd_fail_norm.getVal()*args.rebin_factor);

w_pass.Import(pdf_qcd_pass);
w_fail.Import(pdf_qcd_fail);
w_pass.Import(pdf_qcd_pass_norm);
w_fail.Import(pdf_qcd_fail_norm);

h_fit_qcd_pass_test = pdf_qcd_pass.createHistogram("h_fit_qcd_pass_test",obs);
h_fit_qcd_pass_test.Scale(pdf_qcd_pass_norm.getVal());
chi2_qcd_pass,ndf_qcd_pass = getChi2(h_qcd_pass,h_fit_qcd_pass_test);
chi2_qcd_pass = chi2_qcd_pass/(ndf_qcd_pass-fit_qcd_pass_res.floatParsFinal().getSize());

h_fit_qcd_fail_test = pdf_qcd_fail.createHistogram("h_fit_qcd_fail_test",obs);
h_fit_qcd_fail_test.Scale(pdf_qcd_fail_norm.getVal());
chi2_qcd_fail,ndf_qcd_fail = getChi2(h_qcd_fail,h_fit_qcd_fail_test);
chi2_qcd_fail = chi2_qcd_fail/(ndf_qcd_fail-fit_qcd_fail_res.floatParsFinal().getSize());

label = ROOT.TLatex();
label.SetTextAlign(12);
label.SetNDC();
label.SetTextSize(label.GetTextSize()*0.8);

c = ROOT.TCanvas("c","",600,600);

h_qcd_pass.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_qcd_pass.GetYaxis().SetTitle("Events");
h_qcd_pass.SetMarkerColor(ROOT.kBlack);
h_qcd_pass.SetLineColor(ROOT.kBlack);
h_qcd_pass.SetMarkerSize(0.6);
h_qcd_pass.SetMarkerStyle(20);
h_qcd_pass.Draw("EP");
h_fit_qcd_pass.SetLineColor(ROOT.kRed);
h_fit_qcd_pass.SetLineWidth(2);
h_fit_qcd_pass.Draw("hist same");
h_qcd_pass.Draw("EPsame");

ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_qcd_pass));
#label.DrawLatex(0.65,0.75,"Integral=%.2f"%(pdf_qcd_pass_norm.getVal()));
c.SaveAs(args.output_directory+"/qcd_pass_fit_"+args.cat_id+".png","png");

h_qcd_fail.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_qcd_fail.GetYaxis().SetTitle("Events");
h_qcd_fail.SetMarkerColor(ROOT.kBlack);
h_qcd_fail.SetLineColor(ROOT.kBlack);
h_qcd_fail.SetMarkerSize(0.6);
h_qcd_fail.SetMarkerStyle(20);
h_qcd_fail.Draw("EP");
h_fit_qcd_fail.SetLineColor(ROOT.kRed);
h_fit_qcd_fail.SetLineWidth(2);
h_fit_qcd_fail.Draw("hist same");
h_qcd_fail.Draw("EPsame");
ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_qcd_fail));
#label.DrawLatex(0.65,0.75,"Integral=%.2f"%(pdf_qcd_fail_norm.getVal()));
c.SaveAs(args.output_directory+"/qcd_fail_fit_"+args.cat_id+".png","png");

######################
## Z+jets modelling ##
######################

rh_zjet_pass = ROOT.RooDataHist("zjet_pass","",obs_list,h_zjet_pass)
rh_zjet_fail = ROOT.RooDataHist("zjet_fail","",obs_list,h_zjet_fail)

zjet_mz = ROOT.RooRealVar("zjet_mz","zjet_mz",91.18);
zjet_gammaz = ROOT.RooRealVar("zjet_gammaz","zjet_gammaz",2.49);
zjet_mz.setConstant(True);
zjet_gammaz.setConstant(True);

pdf_zjet_bw = ROOT.RooGenericPdf("pdf_zjet_bw","pdf_zjet_bw","@0/(pow(@0*@0 - @1*@1,2) + @2*@2*@0*@0*@0*@0/(@1*@1))",ROOT.RooArgList(obs,zjet_mz,zjet_gammaz));

zjet_mean_gaus_pass = ROOT.RooRealVar("zjet_mean_gaus_pass","zjet_mean_gaus_pass",0,-10,10);
zjet_sigma_gaus_pass = ROOT.RooRealVar("zjet_sigma_gaus_pass","zjet_sigma_gaus_pass",1,0,10);
zjet_mean_gaus_fail = ROOT.RooRealVar("zjet_mean_gaus_fail","zjet_mean_gaus_fail",0,-10,10);
zjet_sigma_gaus_fail = ROOT.RooRealVar("zjet_sigma_gaus_fail","zjet_sigma_gaus_fail",1,0,10);

pdf_zjet_gaus_pass = ROOT.RooGaussian("pdf_zjet_gaus_pass","pdf_zjet_gaus_pass",obs,zjet_mean_gaus_pass,zjet_sigma_gaus_pass)
pdf_zjet_gaus_fail = ROOT.RooGaussian("pdf_zjet_gaus_fail","pdf_zjet_gaus_fail",obs,zjet_mean_gaus_fail,zjet_sigma_gaus_fail)

pdf_zjet_sig_pass = ROOT.RooFFTConvPdf("pdf_zjet_sig_pass","pdf_zjet_sig_pass",obs,pdf_zjet_bw,pdf_zjet_gaus_pass);
pdf_zjet_sig_fail = ROOT.RooFFTConvPdf("pdf_zjet_sig_fail","pdf_zjet_sig_fail",obs,pdf_zjet_bw,pdf_zjet_gaus_fail);

zjet_coef_pass_1 = ROOT.RooRealVar("zjet_coef_pass_1","zjet_coef_pass_1",0.001,-10,10)
zjet_coef_pass_2 = ROOT.RooRealVar("zjet_coef_pass_2","zjet_coef_pass_2",0.001,-10,10)
zjet_coef_pass_3 = ROOT.RooRealVar("zjet_coef_pass_3","zjet_coef_pass_3",0.001,-10,10)
zjet_coef_pass_4 = ROOT.RooRealVar("zjet_coef_pass_3","zjet_coef_pass_4",0.001,-10,10)
zjet_coef_fail_1 = ROOT.RooRealVar("zjet_coef_fail_1","zjet_coef_fail_1",0.001,-10,10)
zjet_coef_fail_2 = ROOT.RooRealVar("zjet_coef_fail_2","zjet_coef_fail_2",0.001,-10,10)
zjet_coef_fail_3 = ROOT.RooRealVar("zjet_coef_fail_3","zjet_coef_fail_3",0.001,-10,10)
zjet_coef_fail_4 = ROOT.RooRealVar("zjet_coef_fail_4","zjet_coef_fail_4",0.001,-10,10)

if args.cat_id == "cat0":
    pdf_zjet_bkg_pass = ROOT.RooChebychev("pdf_zjet_bkg_pass","",obs,[zjet_coef_pass_1]);
    pdf_zjet_bkg_fail = ROOT.RooChebychev("pdf_zjet_bkg_fail","",obs,[zjet_coef_fail_1,zjet_coef_fail_2,zjet_coef_fail_3,zjet_coef_fail_4]);
else:
    pdf_zjet_bkg_pass = ROOT.RooChebychev("pdf_zjet_bkg_pass","",obs,[zjet_coef_pass_1,zjet_coef_pass_2]);
    pdf_zjet_bkg_fail = ROOT.RooChebychev("pdf_zjet_bkg_fail","",obs,[zjet_coef_fail_1,zjet_coef_fail_2,zjet_coef_fail_3,zjet_coef_fail_4]);

zjet_frac_pass = ROOT.RooRealVar("zjet_frac_pass","zjet_frac_pass",0.1,0.,1.);
zjet_frac_fail = ROOT.RooRealVar("zjet_frac_fail","zjet_frac_fail",0.1,0.,1.);
    
pdf_zjet_pass = ROOT.RooAddPdf("pdf_zjet_pass","pdf_zjet_pass",[pdf_zjet_sig_pass,pdf_zjet_bkg_pass],[zjet_frac_pass],True);
pdf_zjet_fail = ROOT.RooAddPdf("pdf_zjet_fail","pdf_zjet_fail",[pdf_zjet_sig_fail,pdf_zjet_bkg_fail],[zjet_frac_fail],True);

pdf_zjet_pass_norm = ROOT.RooRealVar(pdf_zjet_pass.GetName()+"_norm","",rh_zjet_pass.sumEntries())
pdf_zjet_fail_norm = ROOT.RooRealVar(pdf_zjet_fail.GetName()+"_norm","",rh_zjet_fail.sumEntries())
pdf_zjet_pass_norm.setConstant(True);
pdf_zjet_fail_norm.setConstant(True);

fit_zjet_pass_res = pdf_zjet_pass.fitTo(rh_zjet_pass,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));
fit_zjet_fail_res = pdf_zjet_fail.fitTo(rh_zjet_fail,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));

h_fit_zjet_pass = pdf_zjet_pass.createHistogram("h_fit_zjet_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_fail = pdf_zjet_fail.createHistogram("h_fit_zjet_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_pass.Scale(pdf_zjet_pass_norm.getVal()*args.rebin_factor);
h_fit_zjet_fail.Scale(pdf_zjet_fail_norm.getVal()*args.rebin_factor);

h_fit_zjet_bkg_pass = pdf_zjet_bkg_pass.createHistogram("h_fit_zjet_bkg_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_bkg_fail = pdf_zjet_bkg_fail.createHistogram("h_fit_zjet_bkg_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_bkg_pass.Scale(pdf_zjet_pass_norm.getVal()*args.rebin_factor*(1-zjet_frac_pass.getVal()));
h_fit_zjet_bkg_fail.Scale(pdf_zjet_fail_norm.getVal()*args.rebin_factor*(1-zjet_frac_fail.getVal()));

h_fit_zjet_sig_pass = pdf_zjet_sig_pass.createHistogram("h_fit_zjet_sig_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_sig_fail = pdf_zjet_sig_fail.createHistogram("h_fit_zjet_sig_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_zjet_sig_pass.Scale(pdf_zjet_pass_norm.getVal()*args.rebin_factor*zjet_frac_pass.getVal());
h_fit_zjet_sig_fail.Scale(pdf_zjet_fail_norm.getVal()*args.rebin_factor*zjet_frac_fail.getVal());

h_fit_zjet_pass_test = pdf_zjet_pass.createHistogram("h_fit_zjet_pass_test",obs);
h_fit_zjet_pass_test.Scale(pdf_zjet_pass_norm.getVal());
chi2_zjet_pass,ndf_zjet_pass = getChi2(h_zjet_pass,h_fit_zjet_pass_test);
chi2_zjet_pass = chi2_zjet_pass/(ndf_zjet_pass-fit_zjet_pass_res.floatParsFinal().getSize());
max_zjet_pass, hwhm_zjet_pass = getMaxAndHWHM(h_fit_zjet_sig_pass);

h_fit_zjet_fail_test = pdf_zjet_fail.createHistogram("h_fit_zjet_fail_test",obs);
h_fit_zjet_fail_test.Scale(pdf_zjet_fail_norm.getVal());
chi2_zjet_fail,ndf_zjet_fail = getChi2(h_zjet_fail,h_fit_zjet_fail_test);
chi2_zjet_fail = chi2_zjet_fail/(ndf_zjet_fail-fit_zjet_fail_res.floatParsFinal().getSize());
max_zjet_fail, hwhm_zjet_fail = getMaxAndHWHM(h_fit_zjet_sig_fail);

zjet_coef_pass_1.setConstant(True);
zjet_coef_pass_2.setConstant(True);
zjet_coef_pass_3.setConstant(True);
zjet_coef_pass_4.setConstant(True);
zjet_coef_fail_1.setConstant(True);
zjet_coef_fail_2.setConstant(True);
zjet_coef_fail_3.setConstant(True);
zjet_coef_fail_4.setConstant(True);
zjet_mean_gaus_pass.setConstant(True);
zjet_mean_gaus_fail.setConstant(True);
zjet_sigma_gaus_pass.setConstant(True);
zjet_sigma_gaus_fail.setConstant(True);

w_pass.Import(pdf_zjet_pass);
w_fail.Import(pdf_zjet_fail);
w_pass.Import(pdf_zjet_pass_norm);
w_fail.Import(pdf_zjet_fail_norm);

h_zjet_pass.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_zjet_pass.GetYaxis().SetTitle("Events");
h_zjet_pass.SetMarkerColor(ROOT.kBlack);
h_zjet_pass.SetLineColor(ROOT.kBlack);
h_zjet_pass.SetMarkerSize(0.6);
h_zjet_pass.SetMarkerStyle(20);
h_zjet_pass.Draw("EP");
h_fit_zjet_bkg_pass.SetLineColor(ROOT.kBlue);
h_fit_zjet_bkg_pass.SetLineWidth(2);
h_fit_zjet_bkg_pass.Draw("hist same");
h_fit_zjet_pass.SetLineColor(ROOT.kRed);
h_fit_zjet_pass.SetLineWidth(2);
h_fit_zjet_pass.Draw("hist same");
h_zjet_pass.Draw("EPsame");
ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_zjet_pass));
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_zjet_pass))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_zjet_pass));
#label.DrawLatex(0.65,0.65,"Integral=%.2f"%(pdf_zjet_pass_norm.getVal()));
c.SaveAs(args.output_directory+"/zjet_pass_fit_"+args.cat_id+".png","png");

h_zjet_fail.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_zjet_fail.GetYaxis().SetTitle("Events");
h_zjet_fail.SetMarkerColor(ROOT.kBlack);
h_zjet_fail.SetLineColor(ROOT.kBlack);
h_zjet_fail.SetMarkerSize(0.6);
h_zjet_fail.SetMarkerStyle(20);
h_zjet_fail.Draw("EP");
h_fit_zjet_bkg_fail.SetLineColor(ROOT.kBlue);
h_fit_zjet_bkg_fail.SetLineWidth(2);
h_fit_zjet_bkg_fail.Draw("hist same");
h_fit_zjet_fail.SetLineColor(ROOT.kRed);
h_fit_zjet_fail.SetLineWidth(2);
h_fit_zjet_fail.Draw("hist same");
h_zjet_fail.Draw("EPsame");
ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_zjet_fail));
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_zjet_fail))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_zjet_fail));
#label.DrawLatex(0.65,0.65,"Integral=%.2f"%(pdf_zjet_fail_norm.getVal()));
c.SaveAs(args.output_directory+"/zjet_fail_fit_"+args.cat_id+".png","png");

######################
## W+jets modelling ##
######################

rh_wjet_pass = ROOT.RooDataHist("wjet_pass","",obs_list,h_wjet_pass)
rh_wjet_fail = ROOT.RooDataHist("wjet_fail","",obs_list,h_wjet_fail)

wjet_mw = ROOT.RooRealVar("wjet_mw","wjet_mw",80.37);
wjet_gammaw = ROOT.RooRealVar("wjet_gammaw","wjet_gammaw",2.08);
wjet_mw.setConstant(True);
wjet_gammaw.setConstant(True);

pdf_wjet_bw = ROOT.RooGenericPdf("pdf_wjet_bw","pdf_wjet_bw","@0/(pow(@0*@0 - @1*@1,2) + @2*@2*@0*@0*@0*@0/(@1*@1))",ROOT.RooArgList(obs,wjet_mw,wjet_gammaw));

wjet_mean_gaus_pass = ROOT.RooRealVar("wjet_mean_gaus_pass","wjet_mean_gaus_pass",0,-10,10);
wjet_sigma_gaus_pass = ROOT.RooRealVar("wjet_sigma_gaus_pass","wjet_sigma_gaus_pass",1,0,10);
wjet_mean_gaus_fail = ROOT.RooRealVar("wjet_mean_gaus_fail","wjet_mean_gaus_fail",0,-10,10);
wjet_sigma_gaus_fail = ROOT.RooRealVar("wjet_sigma_gaus_fail","wjet_sigma_gaus_fail",1,0,10);

pdf_wjet_gaus_pass = ROOT.RooGaussian("pdf_wjet_gaus_pass","pdf_wjet_gaus_pass",obs,wjet_mean_gaus_pass,wjet_sigma_gaus_pass)
pdf_wjet_gaus_fail = ROOT.RooGaussian("pdf_wjet_gaus_fail","pdf_wjet_gaus_fail",obs,wjet_mean_gaus_fail,wjet_sigma_gaus_fail)

pdf_wjet_sig_pass = ROOT.RooFFTConvPdf("pdf_wjet_sig_pass","pdf_wjet_sig_pass",obs,pdf_wjet_bw,pdf_wjet_gaus_pass);
pdf_wjet_sig_fail = ROOT.RooFFTConvPdf("pdf_wjet_sig_fail","pdf_wjet_sig_fail",obs,pdf_wjet_bw,pdf_wjet_gaus_fail);

wjet_coef_pass_1 = ROOT.RooRealVar("wjet_coef_pass_1","wjet_coef_pass_1",0.001,-10,10)
wjet_coef_pass_2 = ROOT.RooRealVar("wjet_coef_pass_2","wjet_coef_pass_2",0.001,-10,10)
wjet_coef_pass_3 = ROOT.RooRealVar("wjet_coef_pass_3","wjet_coef_pass_3",0.001,-10,10)
wjet_coef_pass_4 = ROOT.RooRealVar("wjet_coef_pass_3","wjet_coef_pass_4",0.001,-10,10)
wjet_coef_fail_1 = ROOT.RooRealVar("wjet_coef_fail_1","wjet_coef_fail_1",0.001,-10,10)
wjet_coef_fail_2 = ROOT.RooRealVar("wjet_coef_fail_2","wjet_coef_fail_2",0.001,-10,10)
wjet_coef_fail_3 = ROOT.RooRealVar("wjet_coef_fail_3","wjet_coef_fail_3",0.001,-10,10)
wjet_coef_fail_4 = ROOT.RooRealVar("wjet_coef_fail_4","wjet_coef_fail_4",0.001,-10,10)

if args.cat_id == "cat0":
    pdf_wjet_bkg_pass = ROOT.RooChebychev("pdf_wjet_bkg_pass","",obs,[wjet_coef_pass_1]);
    pdf_wjet_bkg_fail = ROOT.RooChebychev("pdf_wjet_bkg_fail","",obs,[wjet_coef_fail_1,wjet_coef_fail_2,wjet_coef_fail_3,wjet_coef_fail_4]);
else:
    pdf_wjet_bkg_pass = ROOT.RooChebychev("pdf_wjet_bkg_pass","",obs,[wjet_coef_pass_1,wjet_coef_pass_2]);
    pdf_wjet_bkg_fail = ROOT.RooChebychev("pdf_wjet_bkg_fail","",obs,[wjet_coef_fail_1,wjet_coef_fail_2,wjet_coef_fail_3,wjet_coef_fail_4]);

wjet_frac_pass = ROOT.RooRealVar("wjet_frac_pass","wjet_frac_pass",0.1,0.,1.);
wjet_frac_fail = ROOT.RooRealVar("wjet_frac_fail","wjet_frac_fail",0.1,0.,1.);
    
pdf_wjet_pass = ROOT.RooAddPdf("pdf_wjet_pass","pdf_wjet_pass",[pdf_wjet_sig_pass,pdf_wjet_bkg_pass],[wjet_frac_pass],True);
pdf_wjet_fail = ROOT.RooAddPdf("pdf_wjet_fail","pdf_wjet_fail",[pdf_wjet_sig_fail,pdf_wjet_bkg_fail],[wjet_frac_fail],True);

pdf_wjet_pass_norm = ROOT.RooRealVar(pdf_wjet_pass.GetName()+"_norm","",rh_wjet_pass.sumEntries())
pdf_wjet_fail_norm = ROOT.RooRealVar(pdf_wjet_fail.GetName()+"_norm","",rh_wjet_fail.sumEntries())
pdf_wjet_pass_norm.setConstant(True);
pdf_wjet_fail_norm.setConstant(True);

fit_wjet_pass_res = pdf_wjet_pass.fitTo(rh_wjet_pass,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));
fit_wjet_fail_res = pdf_wjet_fail.fitTo(rh_wjet_fail,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"));

h_fit_wjet_pass = pdf_wjet_pass.createHistogram("h_fit_wjet_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_fail = pdf_wjet_fail.createHistogram("h_fit_wjet_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_pass.Scale(pdf_wjet_pass_norm.getVal()*args.rebin_factor);
h_fit_wjet_fail.Scale(pdf_wjet_fail_norm.getVal()*args.rebin_factor);

h_fit_wjet_bkg_pass = pdf_wjet_bkg_pass.createHistogram("h_fit_wjet_bkg_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_bkg_fail = pdf_wjet_bkg_fail.createHistogram("h_fit_wjet_bkg_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_bkg_pass.Scale(pdf_wjet_pass_norm.getVal()*args.rebin_factor*(1-wjet_frac_pass.getVal()));
h_fit_wjet_bkg_fail.Scale(pdf_wjet_fail_norm.getVal()*args.rebin_factor*(1-wjet_frac_fail.getVal()));

h_fit_wjet_sig_pass = pdf_wjet_sig_pass.createHistogram("h_fit_wjet_sig_pass",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_sig_fail = pdf_wjet_sig_fail.createHistogram("h_fit_wjet_sig_fail",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_factor));
h_fit_wjet_sig_pass.Scale(pdf_wjet_pass_norm.getVal()*args.rebin_factor*wjet_frac_pass.getVal());
h_fit_wjet_sig_fail.Scale(pdf_wjet_fail_norm.getVal()*args.rebin_factor*wjet_frac_fail.getVal());

h_fit_wjet_pass_test = pdf_wjet_pass.createHistogram("h_fit_wjet_pass_test",obs);
h_fit_wjet_pass_test.Scale(pdf_wjet_pass_norm.getVal());
chi2_wjet_pass,ndf_wjet_pass = getChi2(h_wjet_pass,h_fit_wjet_pass_test);
chi2_wjet_pass = chi2_wjet_pass/(ndf_wjet_pass-fit_wjet_pass_res.floatParsFinal().getSize());
max_wjet_pass, hwhm_wjet_pass = getMaxAndHWHM(h_fit_wjet_sig_pass);

h_fit_wjet_fail_test = pdf_wjet_fail.createHistogram("h_fit_wjet_fail_test",obs);
h_fit_wjet_fail_test.Scale(pdf_wjet_fail_norm.getVal());
chi2_wjet_fail,ndf_wjet_fail = getChi2(h_wjet_fail,h_fit_wjet_fail_test);
chi2_wjet_fail = chi2_wjet_fail/(ndf_wjet_fail-fit_wjet_fail_res.floatParsFinal().getSize());
max_wjet_fail, hwhm_wjet_fail = getMaxAndHWHM(h_fit_wjet_sig_fail);

wjet_coef_pass_1.setConstant(True);
wjet_coef_pass_2.setConstant(True);
wjet_coef_pass_3.setConstant(True);
wjet_coef_pass_4.setConstant(True);
wjet_coef_fail_1.setConstant(True);
wjet_coef_fail_2.setConstant(True);
wjet_coef_fail_3.setConstant(True);
wjet_coef_fail_4.setConstant(True);
wjet_mean_gaus_pass.setConstant(True);
wjet_mean_gaus_fail.setConstant(True);
wjet_sigma_gaus_pass.setConstant(True);
wjet_sigma_gaus_fail.setConstant(True);

w_pass.Import(pdf_wjet_pass);
w_fail.Import(pdf_wjet_fail);
w_pass.Import(pdf_wjet_pass_norm);
w_fail.Import(pdf_wjet_fail_norm);

h_wjet_pass.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_wjet_pass.GetYaxis().SetTitle("Events");
h_wjet_pass.SetMarkerColor(ROOT.kBlack);
h_wjet_pass.SetLineColor(ROOT.kBlack);
h_wjet_pass.SetMarkerSize(0.6);
h_wjet_pass.SetMarkerStyle(20);
h_wjet_pass.Draw("EP");
h_fit_wjet_bkg_pass.SetLineColor(ROOT.kBlue);
h_fit_wjet_bkg_pass.SetLineWidth(2);
h_fit_wjet_bkg_pass.Draw("hist same");
h_fit_wjet_pass.SetLineColor(ROOT.kRed);
h_fit_wjet_pass.SetLineWidth(2);
h_fit_wjet_pass.Draw("hist same");
h_wjet_pass.Draw("EPsame");
ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_wjet_pass));
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_wjet_pass))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_wjet_pass));
#label.DrawLatex(0.65,0.65,"Integral=%.2f"%(pdf_wjet_pass_norm.getVal()));
c.SaveAs(args.output_directory+"/wjet_pass_fit_"+args.cat_id+".png","png");

h_wjet_fail.GetXaxis().SetTitle(args.mass_obs+" (GeV)");
h_wjet_fail.GetYaxis().SetTitle("Events");
h_wjet_fail.SetMarkerColor(ROOT.kBlack);
h_wjet_fail.SetLineColor(ROOT.kBlack);
h_wjet_fail.SetMarkerSize(0.6);
h_wjet_fail.SetMarkerStyle(20);
h_wjet_fail.Draw("EP");
h_fit_wjet_bkg_fail.SetLineColor(ROOT.kBlue);
h_fit_wjet_bkg_fail.SetLineWidth(2);
h_fit_wjet_bkg_fail.Draw("hist same");
h_fit_wjet_fail.SetLineColor(ROOT.kRed);
h_fit_wjet_fail.SetLineWidth(2);
h_fit_wjet_fail.Draw("hist same");
h_wjet_fail.Draw("EPsame");
ROOT.CMS_lumi(c,"",False,False,True);
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_wjet_fail));
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_wjet_fail))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_wjet_fail));
#label.DrawLatex(0.65,0.65,"Integral=%.2f"%(pdf_wjet_fail_norm.getVal()));
c.SaveAs(args.output_directory+"/wjet_fail_fit_"+args.cat_id+".png","png");

os.system("mkdir -p "+args.output_directory);

f_out_pass = ROOT.TFile(args.output_directory+"/workspace_"+args.cat_id+"_pass.root","RECREATE");
f_out_pass.cd()
w_pass.Write("w");
h_qcd_pass.Write();
h_fit_qcd_pass.Write();
f_out_pass.Close();

f_out_fail = ROOT.TFile(args.output_directory+"/workspace_"+args.cat_id+"_fail.root","RECREATE");
f_out_fail.cd()
w_fail.Write("w");
h_qcd_fail.Write();
h_fit_qcd_fail.Write();
f_out_fail.Close();
