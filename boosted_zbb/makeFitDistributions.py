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
parser.add_argument('-i', '--input-file', type=str, default='', help='input directory with files');
parser.add_argument('-c', '--cat-id', type=str, default='', choices=['cat0','cat1','cat2','cat3','cat4'], help='category identifier');
parser.add_argument('-p', '--cut-id', type=str, default='', choices=['pass','fail'], help='identify if is pass or fail');
parser.add_argument('-f', '--pre-fit', action='store_true', help='plot pre-fit instead of post-fit');
parser.add_argument('-o', '--output-directory', type=str, default='', help='name of the output directory');
args = parser.parse_args()

ROOT.gInterpreter.ProcessLine('#include "CMS_style.h"')
ROOT.setTDRStyle();
ROOT.gStyle.SetOptStat(0);

f_input = ROOT.TFile(args.input_file,"READ");

if args.pre_fit:
    h_zjet = f_input.Get("shapes_prefit/"+args.cut_id+"/zjet");
    h_wjet = f_input.Get("shapes_prefit/"+args.cut_id+"/wjet");
    h_qcd = f_input.Get("shapes_prefit/"+args.cut_id+"/qcd");
    h_data = f_input.Get("shapes_prefit/"+args.cut_id+"/data");
    h_tot = f_input.Get("shapes_prefit/"+args.cut_id+"/total");
else:    
    h_zjet = f_input.Get("shapes_fit_s/"+args.cut_id+"/zjet");
    h_wjet = f_input.Get("shapes_fit_s/"+args.cut_id+"/wjet");
    h_qcd = f_input.Get("shapes_fit_s/"+args.cut_id+"/qcd");
    h_data = f_input.Get("shapes_fit_s/"+args.cut_id+"/data");
    h_tot = f_input.Get("shapes_fit_s/"+args.cut_id+"/total");

c = ROOT.TCanvas("c","",600,650);
c.cd();

pad = ROOT.TPad("pad","",0,0,1,1);
pad.SetFillColor(0);
pad.SetFillStyle(0);
pad.SetTickx(1);
pad.SetTicky(1);
pad.SetBottomMargin(0.30);
pad.SetRightMargin(0.06);
pad.Draw();
pad.cd();

hs = ROOT.THStack("hs","");
h_qcd.SetFillColorAlpha(ROOT.kAzure+10,0.5);
h_wjet.SetFillColorAlpha(ROOT.kBlue,0.5);
h_zjet.SetFillColorAlpha(ROOT.kOrange+1,0.5);
h_qcd.SetLineColor(1);
h_wjet.SetLineColor(1);
h_zjet.SetLineColor(1);
hs.Add(h_qcd);
hs.Add(h_wjet);
hs.Add(h_zjet);

h_data.SetMarkerColor(ROOT.kBlack);
h_data.SetLineColor(ROOT.kBlack);
h_data.SetMarkerSize(0.8);
h_data.SetMarkerStyle(20);

h_temp = h_tot.Clone("h_temp");
h_temp.GetYaxis().SetTitle("Events");
h_temp.GetXaxis().SetTitleSize(0);
h_temp.GetXaxis().SetLabelSize(0);
if args.cut_id == "fail":
    h_temp.GetYaxis().SetRangeUser(0.,h_tot.GetMaximum()*1.2);
else:
    h_temp.GetYaxis().SetRangeUser(0.,h_tot.GetMaximum()*1.5);
h_temp.Draw("hist");

hs.Draw("hist same");
h_data.Draw("EP0same");
h_data.Draw("EP1same");
ROOT.CMS_lumi(c,"",False,False,True);
pad.RedrawAxis("sameaxis");
pad.RedrawAxis("g");

leg = ROOT.TLegend(0.6,0.75,0.90,0.90);    
leg.SetFillColor(0);
leg.SetFillStyle(0);
leg.SetBorderSize(0);
leg.AddEntry(h_data,"Data","PE1");
leg.AddEntry(h_zjet,"Z+jets","F");
leg.AddEntry(h_wjet,"W+jets","F");
leg.AddEntry(h_qcd,"QCD","F");
leg.Draw("same");

pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.96);
pad2.SetFillColor(0);
pad2.SetGridy(1);
pad2.SetFillStyle(0);
pad2.SetTickx(1);
pad2.SetTicky(1);
pad2.SetTopMargin(0.71);
pad2.SetBottomMargin(0.10);
pad2.SetRightMargin(0.06);
pad2.Draw();
pad2.cd();

h_ratio = h_data.Clone("h_ratio");
h_tot_clone = h_tot.Clone("h_tot_clone");
for i in range(1,h_tot_clone.GetNbinsX()+1): h_tot_clone.SetBinError(i, 0);
for i in range(0,h_ratio.GetN()):
    h_ratio.SetPoint(i,h_ratio.GetX()[i],h_ratio.GetY()[i]/h_tot_clone.GetBinContent(i+1));
    h_ratio.SetPointError(i,h_data.GetErrorXlow(i),h_data.GetErrorXhigh(i),h_data.GetErrorYlow(i)/h_tot_clone.GetBinContent(i+1),h_data.GetErrorYhigh(i)/h_tot_clone.GetBinContent(i+1));
    
h_ratio.SetLineColor(ROOT.kBlack);
h_ratio.SetMarkerColor(ROOT.kBlack);
h_ratio.SetMarkerSize(0.8);

h_tot_ratio = h_tot.Clone("h_tot_ratio");
for i in range(1,h_tot_ratio.GetNbinsX()+1):
    h_tot_ratio.SetBinContent(i, 1);
    h_tot_ratio.SetBinError(i, 0);

h_tot_ratio.SetMarkerSize(0);
h_tot_ratio.SetLineWidth(2);
h_tot_ratio.SetLineColor(ROOT.kRed);
h_tot_ratio.SetFillColor(0);
h_tot_ratio.GetXaxis().SetLabelSize(0.);
h_tot_ratio.GetYaxis().CenterTitle();
h_tot_ratio.GetYaxis().SetTitle("Data/Bkg.");
h_tot_ratio.GetYaxis().SetTickLength(0.08);
h_tot_ratio.GetXaxis().SetTitle("mass (GeV)")
if args.cut_id == "fail":
    h_tot_ratio.GetYaxis().SetRangeUser(0.8,1.2);
else:
    h_tot_ratio.GetYaxis().SetRangeUser(0.5,1.5);

h_tot_ratio_err = h_tot.Clone("h_tot_ratio_err");
h_tot_ratio_err.Divide(h_tot_clone);
h_tot_ratio_err.SetLineColor(0);
h_tot_ratio_err.SetLineWidth(0);
h_tot_ratio_err.SetFillColor(ROOT.kGray);

h_tot_ratio.GetYaxis().SetNdivisions(506);

h_tot_ratio.Draw("hist");
if not args.pre_fit:
    h_tot_ratio_err.Draw("E2same");
h_tot_ratio.Draw("hist same");
h_ratio.Draw("PE0 same");
h_ratio.Draw("PE1 same");
pad2.RedrawAxis("sameaxis");
pad2.RedrawAxis("g");

os.system("mkdir -p "+args.output_directory);
if args.pre_fit:
    c.SaveAs(args.output_directory+"/distribution_"+args.cat_id+"_"+args.cut_id+"_prefit.png","png");
    c.SaveAs(args.output_directory+"/distribution_"+args.cat_id+"_"+args.cut_id+"_prefit.pdf","pdf");
else:
    c.SaveAs(args.output_directory+"/distribution_"+args.cat_id+"_"+args.cut_id+"_postfit.png","png");
    c.SaveAs(args.output_directory+"/distribution_"+args.cat_id+"_"+args.cut_id+"_postfit.pdf","pdf");
