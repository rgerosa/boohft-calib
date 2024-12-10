import argparse
import copy
import yaml
import json
import ast
import os
import ROOT

ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetPadTopMargin(0.07);
ROOT.gStyle.SetPadBottomMargin(0.13);
ROOT.gStyle.SetPadLeftMargin(0.12);
ROOT.gStyle.SetPadRightMargin(0.05);
ROOT.gStyle.SetTitleFont(42);
ROOT.gStyle.SetTitleColor(1);
ROOT.gStyle.SetTitleTextColor(1);
ROOT.gStyle.SetTitleFillColor(10);
ROOT.gStyle.SetTitleFontSize(0.05);
ROOT.gStyle.SetTitleColor(1, "XYZ");
ROOT.gStyle.SetTitleFont(42, "XYZ");
ROOT.gStyle.SetTitleSize(0.05, "XYZ");
ROOT.gStyle.SetTitleXOffset(0.9);
ROOT.gStyle.SetTitleYOffset(1.05);
ROOT.gStyle.SetLabelColor(1, "XYZ");
ROOT.gStyle.SetLabelFont(42, "XYZ");
ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
ROOT.gStyle.SetLabelSize(0.04, "XYZ");
ROOT.gStyle.SetAxisColor(1, "XYZ");
ROOT.gStyle.SetStripDecimals(1); 
ROOT.gStyle.SetTickLength(0.025, "XYZ");
ROOT.gStyle.SetNdivisions(510, "XYZ");
ROOT.gStyle.SetPadTickX(1); 
ROOT.gStyle.SetPadTickY(1);

parser = argparse.ArgumentParser()
parser.add_argument('--era', '-e', type=str, default="none", help='era to be considered',
                    choices=['2022preEE','2022postEE','2023preBPIX','2023postBPIX'])
parser.add_argument("--xmin",type=float, default=0.8, help="min x-val i.e. PNet Xbb score")
parser.add_argument("--xmax",type=float, default=1.0, help="max x-val i.e. PNet Xbb score")
parser.add_argument("--xerr",type=float, default=0.01, help="error on x-points")
parser.add_argument('--use-fit-band', action='store_true', help='use uncertainty band from the fit');
parser.add_argument('--no-fit', action='store_true', help='no fit up/dw variations');

args = parser.parse_args()   

sf_values = {}

if args.era == "2022preEE":
    sf_values["sf_values_pT_250_350"] = [(0.9925,0.858,-0.106,0.117),(0.9675,1.086,-0.054,0.057),(0.925,0.935,-0.064,0.060),(0.85,0.897,-0.056,0.059)]
    sf_values["sf_values_pT_350_450"] = [(0.9925,0.879,-0.093,0.062),(0.9675,0.976,-0.037,0.038),(0.925,0.950,-0.049,0.041),(0.85,1.004,-0.059,0.050)]
    sf_values["sf_values_pT_450_550"] = [(0.9925,0.952,-0.047,0.049),(0.9675,0.984,-0.039,0.038),(0.925,0.943,-0.049,0.047),(0.85,1.060,-0.057,0.056)]
    sf_values["sf_values_pT_550_700"] = [(0.9925,1.020,-0.051,0.047),(0.9675,0.987,-0.045,0.045),(0.925,1.059,-0.046,0.039),(0.85,0.999,-0.049,0.050)]
    sf_values["sf_values_pT_700"] = [(0.9925,0.943,-0.070,0.093),(0.9675,0.868,-0.069,0.070),(0.925,1.025,-0.080,0.079),(0.85,1.128,-0.073,0.076)]
elif args.era == "2022postEE":
    sf_values["sf_values_pT_250_350"] = [(0.9925,1.130,-0.083,0.065),(0.9675,1.081,-0.033,0.035),(0.925,0.931,-0.030,0.031),(0.85,0.926,-0.033,0.034)]
    sf_values["sf_values_pT_350_450"] = [(0.9925,1.040,-0.033,0.042),(0.9675,1.072,-0.021,0.020),(0.925,0.993,-0.024,0.024),(0.85,0.935,-0.031,0.031)]
    sf_values["sf_values_pT_450_550"] = [(0.9925,1.072,-0.031,0.034),(0.9675,1.005,-0.019,0.019),(0.925,1.004,-0.020,0.020),(0.85,0.962,-0.027,0.027)]
    sf_values["sf_values_pT_550_700"] = [(0.9925,1.077,-0.032,0.033),(0.9675,0.991,-0.026,0.024),(0.925,0.999,-0.026,0.024),(0.85,0.985,-0.035,0.034)]
    sf_values["sf_values_pT_700"] = [(0.9925,1.114,-0.060,0.052),(0.9675,0.965,-0.045,0.046),(0.925,0.951,-0.052,0.054),(0.85,0.945,-0.065,0.063)]
elif args.era == "2023preBPIX":
    sf_values["sf_values_pT_250_350"] = [(0.9925,1.085,-0.073,0.077),(0.9675,0.948,-0.038,0.038),(0.925,0.953,-0.037,0.038),(0.85,0.917,-0.041,0.043)]
    sf_values["sf_values_pT_350_450"] = [(0.9925,0.944,-0.033,0.035),(0.9675,0.945,-0.024,0.025),(0.925,0.899,-0.024,0.024),(0.85,0.916,-0.035,0.035)]
    sf_values["sf_values_pT_450_550"] = [(0.9925,0.884,-0.031,0.032),(0.9675,0.881,-0.024,0.023),(0.925,0.973,-0.023,0.021),(0.85,0.947,-0.030,0.030)]
    sf_values["sf_values_pT_550_700"] = [(0.9925,0.885,-0.040,0.042),(0.9675,0.924,-0.027,0.027),(0.925,0.950,-0.027,0.027),(0.85,0.888,-0.038,0.038)]
    sf_values["sf_values_pT_700"] = [(0.9925,0.965,-0.066,0.070),(0.9675,0.927,-0.061,0.060),(0.925,0.941,-0.084,0.081),(0.85,0.829,-0.098,0.098)]
elif args.era == "2023postBPIX":
    sf_values["sf_values_pT_250_350"] = [(0.9925,1.112,-0.115,0.125),(0.9675,0.990,-0.055,0.049),(0.925,0.889,-0.049,0.053),(0.85,0.926,-0.060,0.062)]
    sf_values["sf_values_pT_350_450"] = [(0.9925,0.902,-0.062,0.065),(0.9675,0.940,-0.028,0.028),(0.925,0.932,-0.033,0.035),(0.85,0.935,-0.049,0.045)]
    sf_values["sf_values_pT_450_550"] = [(0.9925,0.898,-0.059,0.057),(0.9675,0.863,-0.031,0.026),(0.925,0.893,-0.031,0.032),(0.85,0.962,-0.058,0.052)]
    sf_values["sf_values_pT_550_700"] = [(0.9925,0.966,-0.088,0.081),(0.9675,0.836,-0.034,0.036),(0.925,0.954,-0.042,0.044),(0.85,0.985,-0.059,0.060)]
    sf_values["sf_values_pT_700"] = [(0.9925,1.009,-0.095,0.105),(0.9675,1.017,-0.081,0.087),(0.925,0.884,-0.165,0.138),(0.85,0.945,-0.132,0.140)]
    

## buil a graph
sf_graph = [];
sf_graph_up = [];
sf_graph_dw = [];
for key,val in sf_values.items():
    graph = ROOT.TGraphAsymmErrors();
    graph.SetName(key)
    graph_up = ROOT.TGraphAsymmErrors();
    graph_up.SetName(key+"_up")
    graph_dw = ROOT.TGraphAsymmErrors();
    graph_dw.SetName(key+"_dw")
    for ipoint,entry in enumerate(reversed(val)):
        graph.SetPoint(ipoint,entry[0],entry[1]);
        graph_dw.SetPoint(ipoint,entry[0],entry[1]+entry[2]);
        graph_up.SetPoint(ipoint,entry[0],entry[1]+entry[3]);
        err_x_min = args.xerr;
        err_x_max = args.xerr;
        if entry[0]-args.xerr < args.xmin:
            err_x_min = abs(entry[0]-args.xmin)
        if entry[0]+args.xerr > args.xmax:
            err_x_max = abs(entry[0]-args.xmax)            
        graph.SetPointError(ipoint,err_x_min,err_x_max,abs(entry[2]),abs(entry[3]));
        graph_dw.SetPointError(ipoint,err_x_min,err_x_max,abs(entry[2]),abs(entry[3]));
        graph_up.SetPointError(ipoint,err_x_min,err_x_max,abs(entry[2]),abs(entry[3]));
    sf_graph.append(graph);
    sf_graph_up.append(graph_up);
    sf_graph_dw.append(graph_dw);

## make a fit
func = ROOT.TF1("f", "[0]+[1]*x",args.xmin,args.xmax)
func_up = ROOT.TF1("f_up", "[0]+[1]*x",args.xmin,args.xmax)
func_dw = ROOT.TF1("f_dw", "[0]+[1]*x",args.xmin,args.xmax)
canvas = ROOT.TCanvas("canvas","",600,600);

ROOT.gStyle.SetOptFit(1111);

for igraph,graph in enumerate(sf_graph):
    canvas.cd();
    frame  = canvas.DrawFrame(args.xmin,0.7,args.xmax,1.2)
    frame.GetXaxis().SetLabelSize(0.85*frame.GetXaxis().GetLabelSize());
    frame.GetYaxis().SetLabelSize(0.85*frame.GetYaxis().GetLabelSize());
    frame.GetXaxis().SetTitleSize(0.85*frame.GetXaxis().GetTitleSize());
    frame.GetYaxis().SetTitleSize(0.85*frame.GetYaxis().GetTitleSize());
    frame.GetXaxis().SetTitleOffset(1.20);
    frame.GetYaxis().SetTitleOffset(1.20);
    frame.GetXaxis().SetTitle("PNET X_{bb} score");
    frame.GetYaxis().SetTitle("Scale factor bb");
    frame.GetXaxis().SetLimits(graph.GetX()[0]-graph.GetEXlow()[0],args.xmax);
    frame.Draw();
    
    if args.use_fit_band:
        result = graph.Fit(func,"MRNS");    
        interval = ROOT.TGraphErrors()
        step = abs(func.GetXmax()-func.GetXmin())/100;
        for i in range(0,101):
            x = func.GetXmin()+i*step;
            interval.SetPoint(i,x,func.Eval(x));
        
        print("Chi2 ",result.Chi2()," ndf ",result.Ndf()," reduced chi2 ",result.Chi2()/result.Ndf()," p-value ",result.Prob())

        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(interval,0.68);
        interval.SetFillColorAlpha(ROOT.kRed,0.3);
        interval.Draw("E3same");
        func.SetLineColor(ROOT.kRed);
        func.SetLineWidth(2);
        func.SetFillColor(ROOT.kRed);
        func.Draw("Lsame");
        graph.SetMarkerStyle(20);
        graph.SetMarkerSize(0.8);
        graph.SetLineColor(ROOT.kBlack);
        graph.SetMarkerColor(ROOT.kBlack);
        graph.Draw("PEsame");
        frame.Draw("sameaxis");
        canvas.SaveAs(graph.GetName()+"_"+args.era+".png","png");
        canvas.SaveAs(graph.GetName()+"_"+args.era+".pdf","pdf");
    else:
        result = graph.Fit(func,"MRNS");
        print("Central fit: Chi2 ",result.Chi2()," ndf ",result.Ndf()," reduced chi2 ",result.Chi2()/result.Ndf()," p-value ",result.Prob())
        if args.no_fit:
            print("Skip fit");
        else:
            result_up = sf_graph_up[igraph].Fit(func_up,"MRNS");
            print("Up fit: Chi2 ",result_up.Chi2()," ndf ",result_up.Ndf()," reduced chi2 ",result_up.Chi2()/result_up.Ndf()," p-value ",result_up.Prob())
            result_dw = sf_graph_dw[igraph].Fit(func_dw,"MRNS");
            print("Dw fit: Chi2 ",result_dw.Chi2()," ndf ",result_dw.Ndf()," reduced chi2 ",result_dw.Chi2()/result_dw.Ndf()," p-value ",result_dw.Prob())
        sf_band = ROOT.TGraphAsymmErrors()
        step = abs(func.GetXmax()-func.GetXmin())/100;
        for i in range(0,101):
            x = func.GetXmin()+i*step;
            if args.no_fit:
                sf_band.SetPoint(i,x,func.Eval(x));
                sf_band.SetPointError(i,0,0,abs(graph.Eval(x,0,"S")-sf_graph_dw[igraph].Eval(x,0,"S")),abs(graph.Eval(x,0,"S")-sf_graph_dw[igraph].Eval(x,0,"S")));
            else:
                sf_band.SetPoint(i,x,func.Eval(x));
                sf_band.SetPointError(i,0,0,abs(func_dw.Eval(x)-func.Eval(x)),abs(func.Eval(x)-func_up.Eval(x)));
        sf_band.SetFillColorAlpha(ROOT.kRed,0.3);
        sf_band.Draw("E3same");
        func.SetLineColor(ROOT.kRed);
        func.SetLineWidth(2);
        func.SetFillColor(ROOT.kRed);
        func.Draw("Lsame");
        graph.SetMarkerStyle(20);
        graph.SetMarkerSize(0.8);
        graph.SetLineColor(ROOT.kBlack);
        graph.SetMarkerColor(ROOT.kBlack);
        graph.Draw("PEsame");
        frame.Draw("sameaxis");
        canvas.SaveAs(graph.GetName()+"_"+args.era+".png","png");
        canvas.SaveAs(graph.GetName()+"_"+args.era+".pdf","pdf");
