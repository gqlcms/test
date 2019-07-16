#! /usr/bin/env python
import os
import glob
import math
import datetime
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse   import OptionParser
from time       import gmtime, strftime
from array import array
print '\n';
from ROOT import gROOT, TPaveLabel, TPie, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph, TGraphErrors,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, TVectorD, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite

parser = OptionParser()
parser.add_option('--channel',    action="store",type="string",dest="channel"    ,default="mu" )
parser.add_option('--name',       action="store",type="string",dest="name"       ,default="M" )
parser.add_option('--ST_low' ,    action="store",type="int"   ,dest="ST_low"     ,default=1000 )
parser.add_option('--ST_high',    action="store",type="int"   ,dest="ST_high"    ,default=13000)
parser.add_option('--mass_low' ,  action="store",type="int"   ,dest="mass_low"   ,default=0    )
parser.add_option('--mass_high',  action="store",type="int"   ,dest="mass_high"  ,default=13000)
parser.add_option('--Dimass_low' ,action="store",type="int"   ,dest="Dimass_low" ,default=1000 )
parser.add_option('--Dimass_high',action="store",type="int"   ,dest="Dimass_high",default=13000)
parser.add_option('--REGION',     action="store",type="string",dest="REGION"     ,default="PS0")
parser.add_option('--scale1',     action="store",type="int"   ,dest="scale1"     ,default=1    ) #High Mass
parser.add_option('--scale2',     action="store",type="int"   ,dest="scale2"     ,default=1    ) #Low Mass
parser.add_option('--piechart',   action="store",type="int"   ,dest="piechart"   ,default=0    )
parser.add_option('--loop',       action="store",type="int"   ,dest="loop"       ,default=0    )
parser.add_option('--tau',        action="store",type="float" ,dest="tau"        ,default=0.4  )
parser.add_option('--wtag',       action="store",type="string",dest="wtag"      ,default="0.99"  )
parser.add_option('--rtag',       action="store",type="string",dest="rtag"      ,default="0.99"  )
parser.add_option('--Var_N3',     action="store",type="string",dest="Var_N3"    ,default="MJJJ"  )
parser.add_option('--Var_N2',     action="store",type="string",dest="Var_N2"    ,default="MJJc"  )
(options, args) = parser.parse_args()

def UnderOverFlow1D(h):
    Bins=h.GetNbinsX();
    h.SetBinContent( 1,  h.GetBinContent(1)+h.GetBinContent(0) );
    h.SetBinError(   1,  math.sqrt( h.GetBinError(1)*h.GetBinError(1) + h.GetBinError(0)*h.GetBinError(0)) );
    h.SetBinContent( Bins,  h.GetBinContent(Bins)+h.GetBinContent(Bins+1) );
    h.SetBinError(   Bins,  math.sqrt( h.GetBinError(Bins)*h.GetBinError(Bins) + h.GetBinError(Bins+1)*h.GetBinError(Bins+1)) );
    return h;

def Integerization(h):
    Bins=h.GetNbinsX();
    for i in range(1,Bins+1):
        ###half=h.GetBinContent(i)/2.; h.SetBinContent(i,half); #scale psudo-data yields by 0.5
        #print h.GetBinContent(i)
        if (h.GetBinContent(i)-int(h.GetBinContent(i)))>0.5 : value=int(h.GetBinContent(i))+1;
        else : value=int(h.GetBinContent(i));
        h.SetBinContent( i, value          );
        h.SetBinError(   i, math.sqrt(value) );
        #print h.GetBinContent(i)
    return h;

class ANALYSIS:
    def __init__(self, channel , fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):
        self.setTDRStyle();
        self.channel     = channel;
        self.color_palet = {'data':1, 'QCD':2, 'VV':62, 'STop':8, 'TTbar':80, 'ZJets':6, 'WJets':90, 'Signal':1, 'Uncertainty':1, }
        self.vpt_cut   = 200;
        self.Nb_cut    = "==0"
        self.Nak4_cut  = "<5"
        self.tau_cut   = "<"
        self.and_or    = "&&"
        self.MET_cut   = 40; self.lpt_cut= 55;
        if self.channel=="el":
            self.MET_cut= 80;self.lpt_cut= 45;
        self.deltaPhi_METj_cut=2.0;
        if channel=="had":#================ SETININGS FOR 0lep-SEARCH ========================
            if options.REGION[1:3] in ["R1","R2","R3","S2"]     : self.nak8jet2 = 2; self.nak8jet3 = 2;
            if options.REGION[1:3] in ["R4","R5","R6","R7","S3"]: self.nak8jet2 = 3; self.nak8jet3 = 3;
            if options.REGION      in ["PS0","SR"]              : self.nak8jet2 = 2; self.nak8jet3 = 3;
            options.scale1=100000; options.scale2=3000; options.ST_low=1350;
        if channel=="mu" :#================ SETININGS FOR 1lep-SEARCH ========================
            if options.REGION[1:3] in ["R1","R2","R3","S1"]     : self.nak8jet1 = 1; self.nak8jet2 = 1;
            if options.REGION[1:3] in ["R4","R5","R6","R7","S2"]: self.nak8jet1 = 2; self.nak8jet2 = 2;
            if options.REGION      in ["PS0","SR"]              : self.nak8jet1 = 1; self.nak8jet2 = 2;   #use only "SR" for nj = 1 and 2 
            if "tt" in options.REGION                           : self.Nb_cut =">0"; self.Nak4_cut="<7";   #in ["CRtt",   "CR1tt","CR2tt","CR3tt",   "CR4tt","CR5tt","CR6tt","CR7tt" ]:
            if "W"  in options.REGION                           : 
                self.tau_cut= ">";
                if options.REGION in ["CR4W","CR5W","CR6W","CR7W" ]:self.and_or="||";#in ["CRW",   "CR1W","CR2W","CR3W",   "CR4W","CR5W","CR6W","CR7W" ]:
            options.scale1=1000;  options.scale2=30;

    #================ SETININGS FOR Canvas/pads/histos and more ==================
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");        self.tdrStyle.SetCanvasBorderMode(0);        self.tdrStyle.SetCanvasColor(kWhite);        self.tdrStyle.SetCanvasDefH(700);        self.tdrStyle.SetCanvasDefW(700);        self.tdrStyle.SetCanvasDefX(0);          self.tdrStyle.SetCanvasDefY(0);
        #For the Pad:
        self.tdrStyle.SetPadBorderMode(0);        self.tdrStyle.SetPadColor(kWhite);        self.tdrStyle.SetPadGridX(False);        self.tdrStyle.SetPadGridY(False);        self.tdrStyle.SetGridColor(0);        self.tdrStyle.SetGridStyle(3);        self.tdrStyle.SetGridWidth(1);      
        #For the frame:
        self.tdrStyle.SetFrameBorderMode(0);        self.tdrStyle.SetFrameBorderSize(1);        self.tdrStyle.SetFrameFillColor(0);        self.tdrStyle.SetFrameFillStyle(0);        self.tdrStyle.SetFrameLineColor(1);        self.tdrStyle.SetFrameLineStyle(1);        self.tdrStyle.SetFrameLineWidth(1);
        #For the histo:
        self.tdrStyle.SetHistLineColor(1);        self.tdrStyle.SetHistLineStyle(0);        self.tdrStyle.SetHistLineWidth(1);        self.tdrStyle.SetEndErrorSize(2);              self.tdrStyle.SetMarkerStyle(20);      self.tdrStyle.SetErrorX(0.);
        #For the fit/function:
        self.tdrStyle.SetOptFit(1);        self.tdrStyle.SetFitFormat("5.4g");        self.tdrStyle.SetFuncColor(2);        self.tdrStyle.SetFuncStyle(1);        self.tdrStyle.SetFuncWidth(1);      
        #For the date:
        self.tdrStyle.SetOptDate(0);      
        #For the statistics box:
        self.tdrStyle.SetOptFile(0); self.tdrStyle.SetOptStat(0); self.tdrStyle.SetStatColor(kWhite); self.tdrStyle.SetStatFont(42); self.tdrStyle.SetStatFontSize(0.025); self.tdrStyle.SetStatTextColor(1); self.tdrStyle.SetStatFormat("6.4g"); self.tdrStyle.SetStatBorderSize(1); self.tdrStyle.SetStatH(0.1); self.tdrStyle.SetStatW(0.15);
        #Margins:
        self.tdrStyle.SetPadTopMargin(0.05);        self.tdrStyle.SetPadBottomMargin(0.13);        self.tdrStyle.SetPadLeftMargin(0.18);        self.tdrStyle.SetPadRightMargin(0.06);      
        #For the Global title:
        self.tdrStyle.SetOptTitle(0);        self.tdrStyle.SetTitleFont(42);        self.tdrStyle.SetTitleColor(1);        self.tdrStyle.SetTitleTextColor(1);        self.tdrStyle.SetTitleFillColor(10);        self.tdrStyle.SetTitleFontSize(0.05);
        if options.REGION[0:2] is "SR": self.SetTitleSize(0.04, "XYZ");
        #For the axis titles:
        self.tdrStyle.SetTitleColor(1, "XYZ");        self.tdrStyle.SetTitleFont(42, "XYZ");        self.tdrStyle.SetTitleSize(0.06, "XYZ");  
        if options.REGION[0:2] is "SR": self.tdrStyle.SetTitleSize(0.04, "XYZ");
        self.tdrStyle.SetTitleXOffset(0.8);        self.tdrStyle.SetTitleYOffset(0.8);      
        #For the axis labels:
        self.tdrStyle.SetLabelColor(1, "XYZ");        self.tdrStyle.SetLabelFont(42, "XYZ");        self.tdrStyle.SetLabelOffset(0.007, "XYZ");        self.tdrStyle.SetLabelSize(0.04, "XYZ");
        if options.REGION[0:2] is "SR": self.tdrStyle.SetLabelSize(0.02, "XYZ");
        #For the axis:
        self.tdrStyle.SetAxisColor(1, "XYZ");        self.tdrStyle.SetStripDecimals(kTRUE);        self.tdrStyle.SetTickLength(0.03, "XYZ");        self.tdrStyle.SetNdivisions(510, "XYZ");        self.tdrStyle.SetPadTickX(1);       self.tdrStyle.SetPadTickY(1);      
        #Change for log plots:
        self.tdrStyle.SetOptLogx(0); self.tdrStyle.SetOptLogy(0); self.tdrStyle.SetOptLogz(0);
        #Postscript options:
        self.tdrStyle.SetPaperSize(20.,20.); self.tdrStyle.cd();

    
    def DefineSelection_0lep(self):#==========[ 0lep REGIONS & SELECTION DEFINITION ]===========================================
        REGION=options.REGION; Mjmax="Mj_maxc"; Mjmid="Mj_mid"; Mjmin="Mj_minc";  MJJc="MJJc"; MJJJ="MJJJ"; t1=str(options.tau);t2="0.4"; DM_max_min="("+Mjmax+"-"+Mjmin+")"; DM_max_mid="("+Mjmax+"-"+Mjmid+")"; DM_mid_min="("+Mjmid+"-"+Mjmin+")";
        t41max = "tau41_max"; t41mid = "tau41_mid"; t41min = "tau41_min";
        t31max = "tau31_max"; t31mid = "tau31_mid"; t31min = "tau31_min";
        t21max = "tau21_max"; t21mid = "tau21_mid"; t21min = "tau21_min";  
        t42max = "tau42_max"; t42mid = "tau42_mid"; t42min = "tau42_min";  
        DPhijj = "Phij_12>0.7 && Phij_13>0.7 && Phij_23>0.7  ";
        EtaJCut= " abs(Etaj)<2. && abs(Etaj_2)<2. ";
        Filters= "  1  &&  ";
        c="0.4"; g="0.25"; g2=options.wtag; g3=options.rtag; 
     #   Filters= "  passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0  &&  ";

        common    = "  %s<(ST)&&(ST)<%s"%(options.ST_low,options.ST_high)+"  &&  40<Mj_max  &&  ((Nj8==4&&MJJJ>1000)||( Nj8==2 && 1000<"+MJJc+" && 40<"+Mjmin+") || ( Nj8==3 && 1000<"+MJJJ+" && 40<"+Mjmid+" )) && "; #Mj_maxc and Mj_minc
      #  common     = " %s<(ST)&&(ST)<%s"%(options.ST_low,options.ST_high)+" &&  40<Mj_max  &&  (( Nj8==2 && 1000<MJJ && 40<Mj_min) || ( Nj8>=3 && 1000<MJJJ && 40<Mj_mid )) &&  ";
        common1    = "  %s<(ST)&&(ST)<%s"%(options.ST_low,options.ST_high)+"  &&  40<Mj_max  &&  ((Nj8==4&&MJJJ>1000)||( Nj8==2 && 1000<"+MJJc+" && 40<Mj_min) || ( Nj8==3 && 1000<"+MJJJ+" && 40<"+Mjmid+" )) && "
        Subjtns2  = " (("+Mjmax+"<100&&"+t41max+"<0.25) || ("+Mjmax+">100&&"+t42max+"<0.4))   &&   (("+Mjmin+"<100&&"+t41min+"<0.25) || ("+Mjmin+">100&&"+t42min+"<0.4))   &&   "; 
        Subjtns3  = " ("+t21max+"<"+t1+")   &&   ("+t21mid+"<"+t1+")  &&   ("+t21min+"<"+t1+") && "; 
    #   Subjtns3  = " (("+Mjmax+"<130&&"+t21max+"<"+t1+") || ("+Mjmax+">130&&"+t42max+"<0.4))   &&   (("+Mjmid+"<130&&"+t21mid+"<"+t1+") || ("+Mjmid+">130&&"+t42mid+"<0.4))   &&   (((40<"+Mjmin+"&&"+Mjmin+"<130&&"+t21min+"<"+t1+") || ("+Mjmin+">130&&"+t42min+"<0.4)) || ("+Mjmin+"<40&&"+t21max+"<1.45&&"+t21mid+"<1.4)) && ";
#        Subjtns3  = " (("+Mjmax+"<100&&"+t41max+"<0.25) || ("+Mjmax+">100&&"+t42max+"<0.4))   &&   (("+Mjmid+"<100&&"+t41mid+"<0.25) || ("+Mjmid+">100&&"+t42mid+"<0.4))   &&   (((40<"+Mjmin+"&&"+Mjmin+"<100&&"+t41min+"<0.25) || ("+Mjmin+">100&&"+t42min+"<0.4)) || ("+Mjmin+"<40&&"+t41max+"<1.45&&"+t41mid+"<1.4)) && ";
        Constrains= "";
#        Constrains= "&&PTj>500&&PTj_2>500&&HT>1200 "#  %s<"%(options.mass_low)+Mjmin+"&&"+Mjmin+"<%s    &&  (( Nj8==2 && %s<"%(options.mass_high,options.Dimass_low)+MJJc+"&&"+MJJc+"<%s) || ( Nj8==3 && %s<"%(options.Dimass_high,options.Dimass_low)+MJJJ+"&&"+MJJJ+"<%s)) && "%(options.Dimass_high);
        #-------------------------- PRESELECTION -----------------------------------------------------
        if REGION == "PS0" : selection = Filters + "" +            Constrains+"(1>0)";
        if REGION == "PS2" : selection = Filters + MJJc+">-1000" +            Constrains +" &&Nj8==2";
        if REGION == "PS3" : selection = Filters + MJJJ+">1000" +            Constrains +" &&Nj8==3";
        #------------------------------ SIGNAL REGIONS -----------------------------------------------
        ISR1 = common1 + " Nj8==3  &&  num_ak4jetsex<=2  &&   70<Mj_max2&&Mj_max2<100  &&  70<Mj_min2&&Mj_min2<100    &&  PTj_3/PTj<0.2                 &&  tau42_max2<"+c+"  &&  tau21_min2<"+c+"  &&  tau41_max2<"+g+" && tau41_min2<"+g+"  &&  nbtag_deep==0  ";
        ISR2 = common1 + " Nj8==3  &&  num_ak4jetsex<=2  &&  100<Mj_max2&&Mj_max2<200  &&  70<Mj_min2&&Mj_min2<100    &&  PTj_3/PTj<0.2                  &&  tau42_max2<"+c+"  &&  tau21_min2<"+c+"  &&  tau41_max2<"+g+" && tau41_min2<"+g+"  &&  nbtag_deep==0  ";
        ISR3 = common1 + " Nj8==3  &&  num_ak4jetsex<=2  &&  200<Mj_max2              &&  70<Mj_min2&&Mj_min2<100    &&  PTj_3/PTj<0.2                  &&  tau42_max2<"+c+"  &&  tau21_min2<"+c+"  &&  tau41_max2<"+g+" && tau41_min2<"+g+"  &&  nbtag_deep==0  ";
        ISR4 = common1 + " Nj8==4   &&  num_ak4jetsex<=2 && PTj_4/PTj<0.2   &&   70<Mj_mid&&Mj_max<100  &&  70<Mj_min                                                        &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  nbtag_deep==0  ";
        ISR5 = common1 + " Nj8==4   &&  num_ak4jetsex<=2 && PTj_4/PTj<0.2   &&    70<Mj_mid&&Mj_max<100  &&     Mj_min<70                                                     &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  nbtag_deep==0  ";
        ISR42 = common + " Nj8==3   &&  num_ak4jetsex<=2 && PTj_min/PTj_mid<0.2   &&   70<Mj_mid&&Mj_max<100  &&  60<Mj_min                                                        &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";
        ISR52 = common + " Nj8==3   &&  num_ak4jetsex<=2 && PTj_min/PTj_mid<0.2   &&    70<Mj_mid&&Mj_max<100  &&     Mj_min<60                                                     &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";

        ICR1 = common1 + " Nj8==3  &&  num_ak4jetsex<=2 &&  PTj_3/(ST-PTj_3)<0.1  && 70<Mj_max&&Mj_max<100  &&  40<Mj_min&&Mj_min<70                                             &&  tau21_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  "
        ICR2 = common1 + " Nj8==3  &&  num_ak4jetsex<=2  && PTj_3/(ST-PTj_3)<0.1  &&    ((40<Mj_min&&Mj_min<70)||(100<Mj_min&&Mj_min<130))  &&  100<Mj_max&&Mj_max<200                &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";
 
        ICR3 = common1 + " Nj8==3  &&  num_ak4jetsex<=2  && PTj_3/(ST-PTj_3)<0.1    && ((40<Mj_min&&Mj_min<70)||(100<Mj_min&&Mj_min<130))  &&  200<Mj_max                            &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";

        ICR4 = common1 + " Nj8==4   &&  num_ak4jetsex<=2 && PTj_4/(ST-PTj_4)<0.1   && ((100<Mj_max&&(70<Mj_mid&&Mj_mid<100))||(Mj_mid<70&&(70<Mj_max&&Mj_max<100))) && 60<Mj_min    &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";

        ICR5 = common1 + " Nj8==4   &&  num_ak4jetsex<=2 && PTj_4/(ST-PTj_4)<0.1  && ((100<Mj_max&&(70<Mj_mid&&Mj_mid<100))||(Mj_mid<70&&(70<Mj_max&&Mj_max<100))) &&    Mj_min<60 &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";   

    #    if REGION == "SR1" : selection = Filters + common + Subjtns2 + Constrains +" Nj8==2 && num_ak4jetsex<=2 &&  70<"+Mjmax+"&&"+Mjmax+"<100  &&  70<"+Mjmin+"&&"+Mjmin+"<100 ";
        if REGION == "SR1" :
           #       selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  70<Mj_maxc&&Mj_maxc<100  &&  70<Mj_minc&&Mj_minc<100                                            &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+" &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";
                  selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  70<Mj_maxc&&Mj_maxc<100  &&  70<Mj_minc&&Mj_minc<100                                            &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_min>"+options.wtag+"  &&  nbtag_deep==0  ";
        #          selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  70<Mj_maxc&&Mj_maxc<100  &&  70<Mj_minc&&Mj_minc<100                                            &&  tau21_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";
                  selection="("+selection1+") || ("+ISR1+")"; 

        if REGION == "SR2" :
                  selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  100<Mj_maxc&&Mj_maxc<200  &&  70<Mj_minc&&Mj_minc<100                                            &&  (jetAK8puppi_dnnDecorrh4q_max + jetAK8puppi_dnnDecorrw_max)/(jetAK8puppi_dnnDecorrh4q_max + jetAK8puppi_dnnDecorrw_max + jetAK8puppi_dnnDecorrqcd_max)>"+g3+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+"  &&  nbtag_deep==0  "; 
      #            selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  100<Mj_maxc&&Mj_maxc<200  &&  70<Mj_minc&&Mj_minc<100                                            &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";
                  selection="("+selection1+") || ("+ISR2+")";
        if REGION == "SR3" : 
                  selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  200<Mj_maxc              &&  70<Mj_minc&&Mj_minc<100                                            &&  (jetAK8puppi_dnnDecorrh4q_max + jetAK8puppi_dnnDecorrw_max)/(jetAK8puppi_dnnDecorrh4q_max + jetAK8puppi_dnnDecorrw_max + jetAK8puppi_dnnDecorrqcd_max)>"+g3+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+"  &&  nbtag_deep==0  "; 
#                  selection1 = Filters + common + " Nj8==2  &&  num_ak4jetsex<=2  &&  200<Mj_maxc              &&  70<Mj_minc&&Mj_minc<100                                            &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";

                  selection="("+selection1+") || ("+ISR3+")";
        if REGION == "SR4" :
                  selection1 = Filters + common + " Nj8==3   &&  num_ak4jetsex<=2  &&   70<Mj_mid&&Mj_max<100  &&  70<Mj_min                                                        &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_mid>"+g2+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+"   &&  nbtag_deep==0 "; 
 #                 selection1 = Filters + common + " Nj8==3   &&  num_ak4jetsex<=2  &&   70<Mj_mid&&Mj_max<100  &&  60<Mj_min                                                        &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";

                  selection="("+selection1+") || ("+ISR4+")";
        if REGION == "SR5" : 
                  selection1 = Filters + common + " Nj8==3   &&  num_ak4jetsex<=2  &&   70<Mj_mid&&Mj_max<100  &&     Mj_min<70                                                     &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_mid>"+g2+"   &&  nbtag_deep==0  &&  PTj_3/PTj>=0.2"; 
#                  selection1 = Filters + common + " Nj8==3   &&  num_ak4jetsex<=2  &&   70<Mj_mid&&Mj_max<100  &&     Mj_min<60                                                     &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";
                  selection="("+selection1+") || ("+ISR5+")";
             #     selection="(("+selection1+")||("+ISR5+"))&& !("+ISR52+")";
        #------------------------- CONTROL REGIONS for QCD -------------------------------------------
        if REGION == "CR1q": selection1 = common + " Nj8==2  &&  num_ak4jetsex<=2  &&   70<Mj_max&&Mj_max<100  &&  40<Mj_min&&Mj_min<70                                             &&  tau21_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";selection="("+selection1+") || ("+ICR1+")";
        if REGION == "CR2q": selection1 = common + " Nj8==2  &&  num_ak4jetsex<=2  && ((40<Mj_min&&Mj_min<70)||(100<Mj_min&&Mj_min<130))  &&  100<Mj_max&&Mj_max<200                &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";selection="("+selection1+") || ("+ICR2+")";
        if REGION == "CR3q": selection1 = common + " Nj8==2  &&  num_ak4jetsex<=2  && ((40<Mj_min&&Mj_min<70)||(100<Mj_min&&Mj_min<130))  &&  200<Mj_max                            &&  tau42_max<"+c+"  &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight==0  ";selection="("+selection1+") || ("+ICR3+")";


     #   if REGION == "CR4q": selection1 = common + " Nj8==3  &&  num_ak4jetsex<=2  && ((100<Mj_max&&60<Mj_mid&&Mj_mid<100))||(Mj_mid<60&&(60<Mj_max&&Mj_max<100))) && 60<Mj_min    &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_mid>"+g2+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+"  &&  num_bJet_tight==0  ";     selection="(("+selection1+") || ("+ICR4+"))&&(!("+ICR1+"))&&(!("+ICR2+"))&&(!("+ICR3+"))";
        #if REGION == "CR4q": selection1 = common + " Nj8==3  &&  num_ak4jetsex<=2  && ((100<Mj_max&&60<Mj_mid&&Mj_mid<100))||(Mj_mid<60&&(60<Mj_max&&Mj_max<100))) && 60<Mj_min    &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";     selection="(("+selection1+") || ("+ICR4+"))&&(!("+ICR1+"))&&(!("+ICR2+"))&&(!("+ICR3+"))";

        if REGION == "CR4q": selection1 = common + " Nj8==3  &&  num_ak4jetsex<=2  && ((100<Mj_max&&(60<Mj_mid&&Mj_mid<100))||(Mj_mid<60&&(60<Mj_max&&Mj_max<100))) && 60<Mj_min    &&  jetAK8puppi_dnnDecorrW_max>"+g2+"  &&  jetAK8puppi_dnnDecorrW_mid>"+g2+"  &&  jetAK8puppi_dnnDecorrW_min>"+g2+"  &&  num_bJet_tight==0  ";     selection="(("+selection1+") || ("+ICR4+"))&&(!("+ICR1+"))&&(!("+ICR2+"))&&(!("+ICR3+"))";

#        if REGION == "CR4q": selection = Filters + common + " Nj8==3   &&  num_ak4jetsex<=2  &&   60<Mj_mid&&Mj_max<100  &&  60<Mj_min                                                        &&  jetAK8puppi_dnnDecorrW_max<"+g2+"  &&  jetAK8puppi_dnnDecorrW_mid<"+g2+"  &&  jetAK8puppi_dnnDecorrW_min<"+g2+"   &&  nbtag_deep==0 ";


        if REGION == "CR5q": selection1 = common + " Nj8==3  &&  num_ak4jetsex<=2  && ((100<Mj_max&&(70<Mj_mid&&Mj_mid<100))||(Mj_mid<70&&(70<Mj_max&&Mj_max<100))) &&    Mj_min<60 &&  tau21_max<"+c+"  &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight==0  ";     selection="(("+selection1+") || ("+ICR5+"))&&(!("+ICR1+"))&&(!("+ICR2+"))&&(!("+ICR3+"))";
        if REGION == "CR2t": selection = common + " Nj8==2  &&  num_ak4jetsex<=4  && ((140<Mj_max&&Mj_max<200&&tau42_max<"+c+")||(60<Mj_max&&Mj_max<100&&tau21_max<"+c+"))  &&  60<Mj_min&&Mj_min<100   &&  tau21_min<"+c+"  &&  tau41_max<"+g+" && tau41_min<"+g+"  &&  num_bJet_tight>=1  ";
        if REGION == "CR5t": selection = common + " Nj8==3  &&  num_ak4jetsex<=4  && ((140<Mj_max&&Mj_max<200&&tau42_max<"+c+")||(60<Mj_max&&Mj_max<100&&tau21_max<"+c+"))  &&  60<Mj_mid&&Mj_mid<100   &&  tau21_mid<"+c+"  &&  tau41_max<"+g+" && tau41_mid<"+g+"  &&  num_bJet_tight>=1  ";
   #     if REGION == "CR1q": selection = Filters + common + Subjtns2 + Constrains +" Nj8==2 && Nj4<=5 &&                     "+DM_max_min+"< 30 && ((40<"+Mjmin+"&&"+Mjmin+"<70)||(100<"+Mjmin+"&&"+Mjmin+"<130))  ";
    #    if REGION == "CR2q": selection = Filters + common + Subjtns2 + Constrains +" Nj8==2 && Nj4<=5 &&  30<"+DM_max_min+"&&"+DM_max_min+"<130 && ((40<"+Mjmin+"&&"+Mjmin+"<70)||(100<"+Mjmin+"&&"+Mjmin+"<130))  &&  (nbtag==0||200<"+Mjmax+"||"+Mjmax+"<150)  ";
     #   if REGION == "CR3q": selection = Filters + common + Subjtns2 + Constrains +" Nj8==2 && Nj4<=5 && 130<"+DM_max_min+"                     && ((40<"+Mjmin+"&&"+Mjmin+"<70)||(100<"+Mjmin+"&&"+Mjmin+"<130))  ";
      #  if REGION == "CR4q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 &&  50<"+Mjmin+" && "+Mjmax+"<120"+"&&"+Mjmid+">50"+ " &&  !(70<"+Mjmid+"&&"+Mjmax+"<100)  ";
       # if REGION == "CR5q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 &&  50>"+Mjmin+" && "+Mjmax+"<120"+"&&"+Mjmid+">50"+ " &&  !(70<"+Mjmid+"&&"+Mjmax+"<100)  ";
#        if REGION == "CR4q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 && (("+Mjmid+" >70&& "+Mjmax+"<130  &&"+Mjmax+">100)||("+Mjmid+" <70&& "+Mjmid+" >40&& ("+Mjmax+"-"+Mjmid+")<60))&&"+Mjmin+"<40";
 #       if REGION == "CR5q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 && (("+Mjmid+" >70&& "+Mjmax+"<130  &&"+Mjmax+">100)||("+Mjmid+" <70&& "+Mjmid+" >40&& ("+Mjmax+"-"+Mjmid+")<60))&&"+Mjmin+">40";
        #if REGION == "CR6q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 &&  30<"+DM_max_mid+"&&"+DM_max_mid+"<130 && ((40<"+Mjmin+"&&"+Mjmin+"<70)||(100<"+Mjmin+"&&"+Mjmin+"<130)) && (nbtag==0||(100<"+Mjmax+"&&"+Mjmax+"<150))";
        #if REGION == "CR7q": selection = Filters + common + Subjtns3 + Constrains +" Nj8==3 && Nj4<=6 && 130<"+DM_max_mid+"                     && ((40<"+Mjmin+"&&"+Mjmin+"<70)||(100<"+Mjmin+"&&"+Mjmin+"<130)) ";
        #-------------------------- CONTROL REGION SR2-top -------------------------------------------
      #  if REGION == "CR2t": selection = Filters + common + Subjtns2 + Constrains +" Nj8==2 && Nj4<=5  && 150<"+Mjmax+"&&"+Mjmax+"<200  &&  70<"+Mjmin+"&&"+Mjmin+"<100 && "+t21min+"<0.4 && nbtag>=1  ";
        #---------------------------------------------------------------------------------------------
        self.Make_Controlplots_for_0lep( selection , "" , self.nak8jet2 ); #FOR eatch SRi  #options.numak8jets );  #self.Make_Controlplots_for_0lep(selection,"preselection",self.nak8jet2);

    def Make_Controlplots_for_0lep(self,selection,tag,Nj,CR=0): #===============================================================
        REGION = options.REGION;
        if REGION[1:3] in ["R1","R2","R3",    "S2","S0"]: Nj=2;
        if REGION[1:3] in ["R4","R5","S3","S0"]: Nj=3;
        if REGION!="PS0" and REGION!="PS2" and REGION!="PS3":  
           if REGION[1:3] in ["R1","R2","R3",    "S2","S0"]:  
             #     self.construct_plot(Nj,"Nj4",selection,"",tag,10,-0.5,9.5, "Number of AK4Jets","Events",0,CR);
              #    self.construct_plot(Nj,"num_ak4jetsex",selection,"",tag,10,-0.5,9.5, "Number of exclusive AK4Jets","Events",0,CR);
                  self.construct_plot(Nj,"(MJJc*(Nj8==2)+MJJ*(Nj8==3))",selection,"",tag,30,1000,4000, "M(JJc) (GeV)","Events/(100 GeV)",0,CR);
             #     self.construct_plot(Nj,"Mj_maxc",selection,"",tag,30,70,200, "M(JJc) (GeV)","Events/(100 GeV)",0,CR);
               #   self.construct_plot(Nj,"MJJ",selection,"",tag,30,1000,4000, "M(JJ) (GeV)","Events/(100 GeV)",0,CR);
           if REGION[1:3] in ["R4","R5","S3","S0"]:       
                #  self.construct_plot(Nj,"Nj4",selection,"",tag,10,-0.5,9.5, "Number of AK4Jets","Events",0,CR);
                 # self.construct_plot(Nj,"num_ak4jetsex",selection,"",tag,10,-0.5,9.5, "Number of exclusive AK4Jets","Events",0,CR);
                  self.construct_plot(Nj,"MJJJ",selection,"",tag,30,1000,4000, "M(JJJc) (GeV)","Events/(200 GeV)",0,CR);
                  #self.construct_plot(Nj,options.Var_N3,selection,"",tag,20,0,100, "M(JJJc) (GeV)","Events/(5 GeV)",0,CR);
        #          self.construct_plot(Nj,"Mj_min",selection,"",tag,20,0,100, "Mj_min (GeV)","Events/(5 GeV)",0,CR);

             #     self.construct_plot(Nj,"PTj_4/ST"    ,selection,"",tag,20,0,0.2,"PTj_4/ST","Events",0,CR);
        else :
           if REGION=="PS0":
              Nj=23; #useless  
           if REGION=="PS2":
              Nj=2;
              Phiptmax="jetAK8puppi_phi";
              Phiptmin="jetAK8puppi_phi_2";
           if REGION=="PS3":
              Nj=3;
              Phiptmax="jetAK8puppi_phi";
              Phiptmin="jetAK8puppi_phi_3";
           if REGION=="PS2" or REGION=="PS3":
              deltaPhiMETmax="pi-fabs(pi-fabs(MET_phi-"+Phiptmax+"))";
              deltaPhiMETmin="pi-fabs(pi-fabs(MET_phi-"+Phiptmin+"))";
     #      self.construct_plot(Nj,"nbtag"                  ,selection,"",tag,5,   -0.5, 4.5,"num of bjets(csv_v2)"     ,"Events/bin"    ,0 , CR);
           self.construct_plot(Nj,"nbtag_deep"                  ,selection,"",tag,5,   -0.5, 4.5,"num of bjets(csv_deep)"     ,"Events/bin"    ,0 , CR);
#           self.construct_plot(Nj,"Mj"                  ,selection,"",tag,30,   0, 300,"Mj (GeV)"     ,"Events/(10 GeV)"    ,0 , CR);
#           self.construct_plot(Nj,"Mj_2"                  ,selection,"",tag,30,   0, 300,"Mj_2 (GeV)"     ,"Events/(10 GeV)"    ,0 , CR);
         #  self.construct_plot(Nj,"MET_et*cos("+deltaPhiMETmax+")"                  ,selection,"",tag,30,   0, 300,"MET parallel to Jet^{ptmax} (GeV)"     ,"Events/(10 GeV)"    ,0 , CR); 
         #  self.construct_plot(Nj,"MET_et*sin("+deltaPhiMETmax+")"                  ,selection,"",tag,30,   0, 300,"MET perpendicular to Jet^{ptmax} (GeV)"     ,"Events/(10 GeV)"    ,0 , CR); 
         #  self.construct_plot(Nj,"MET_et*cos("+deltaPhiMETmin+")"                  ,selection,"",tag,30,   0, 300,"MET parallel to Jet^{ptmin} (GeV)"     ,"Events/(10 GeV)"    ,0 , CR); 
        #   self.construct_plot(Nj,"MET_et*sin("+deltaPhiMETmin+")"                  ,selection,"",tag,30,   0, 300,"MET perpendicular to Jet^{ptmin} (GeV)"     ,"Events/(10 GeV)"    ,0 , CR);

#self.construct_plot(Nj,"(jet_tau2tau1_puppi_mmax+jet_tau3tau1_puppi_mmax+jet_tau4tau1_puppi_mmax)/3" ,selection,"",tag,20,   0,   1,"(#tau_{21}+#tau_{31}+#tau_{41})^{mmax} / 3"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_max"      ,selection,"",tag,36,  40, 220,"M_{j}^{max} (GeV)","Events/(5 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_mid"      ,selection,"",tag,30,   0, 150,"M_{j}^{mid} (GeV)","Events/(5 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_min"      ,selection,"",tag,30,   0, 150,"M_{j}^{min} (GeV)","Events/(5 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_mean"     ,selection,"",tag,40,   0, 100,"<m_{j}> (GeV)"    ,"Events/(10 GeV)",0 , CR);
        #if (options.PlotNum==3 or options.PlotNum==0):
        #self.construct_plot(Nj,"ST+MET_et"               ,selection,"",tag,25, 500,3000,"S_{T} (GeV)"      ,"Events/(100 GeV)",0 , CR );
        #self.construct_plot(Nj,"jet_tau2tau1_puppi_mmax" ,selection,"",tag,20,  0.15,   0.8,"#tau_{21_mmax}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi_mmax" ,selection,"",tag,20,   0.1,   0.6,"#tau_{31_mmax}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi_mmax" ,selection,"",tag,20,   0.1,   0.5,"#tau_{41_mmax}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi_mmax" ,selection,"",tag,20,   0.1,   1,"#tau_{42_mmax}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi_mmid" ,selection,"",tag,20,   0.15,   0.8,"#tau_{21_mmid}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi_mmid" ,selection,"",tag,20,   0.1,   0.6,"#tau_{31_mmid}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi_mmid" ,selection,"",tag,20,   0.1,   0.5,"#tau_{41_mmid}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi_mmid" ,selection,"",tag,20,   0.1,   1,"#tau_{42_mmid}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi_mmin" ,selection,"",tag,20,   0.15,   0.8,"#tau_{21_mmin}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi_mmin" ,selection,"",tag,20,   0.1,   0.6,"#tau_{31_mmin}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi_mmin" ,selection,"",tag,20,   0.1,   0.5,"#tau_{41_mmin}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi_mmin" ,selection,"",tag,20,   0.1,   1,"#tau_{42_mmin}"   ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi"          ,selection,"",tag,22,   0, 220,"m_{j1} (GeV)"     ,"Events/(10 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_2"        ,selection,"",tag,22,   0, 220,"m_{j2} (GeV)"     ,"Events/(10 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_mass_puppi_3"        ,selection,"",tag,22,   0, 220,"m_{j3} (GeV)"     ,"Events/(10 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_pt_puppi"            ,selection,"",tag,30,   0,1500,"PT_{j1} (GeV)"    ,"Events/(50 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_pt_puppi_2"          ,selection,"",tag,30,   0,1500,"PT_{j2} (GeV)"    ,"Events/(50 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jet_pt_puppi_3"          ,selection,"",tag,30,   0,1500,"PT_{j3} (GeV)"    ,"Events/(50 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"jetAK8puppi_eta"         ,selection,"",tag,25,   0, 2.5,"#eta_{j1}"        ,"Events(0.1)"     ,0 , CR);        
        #self.construct_plot(Nj,"jetAK8puppi_eta_2"       ,selection,"",tag,25,   0, 2.5,"#eta_{j2}"        ,"Events(0.1)"     ,0 , CR);
        #self.construct_plot(Nj,"jetAK8puppi_eta_3"       ,selection,"",tag,25,   0, 2.5,"#eta_{j3}"        ,"Events(0.1)"     ,0 , CR);
        #self.construct_plot(Nj,"jet_pt_puppi_23"         ,selection,"",tag,30,   0,1500,"PT_{j2+j3} (GeV)" ,"Events/(50 GeV)" ,0 , CR);
        #self.construct_plot(Nj,"MassVV[0]"               ,selection,"",tag,40,   0,4000,"m_{j1j2} (GeV)"   ,"Events/(100 GeV)",0 , CR);
        #self.construct_plot(Nj,"MassVV[1]"               ,selection,"",tag,40,   0,4000,"m_{j3j1} (GeV)"   ,"Events/(100 GeV)",0 , CR);
        #self.construct_plot(Nj,"MassVV[2]"               ,selection,"",tag,40,   0,4000,"m_{j2j3} (GeV)"   ,"Events/(100 GeV)",0 , CR);
        #self.construct_plot(Nj,"Nj4"                     ,selection,"",tag,10,-0.5, 9.5,"Number of jetAK4" ,"Events"          ,0, CR );
        #self.construct_plot(Nj,"Nj8"                     ,selection,"",tag, 5,-0.5, 4.5,"Number of jetAK8" ,"Events"          ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_d12"     ,selection,"",tag,32,   0, 3.2,"jet_phi12"        ,"Events/(0.1)"    ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_d13"     ,selection,"",tag,32,   0, 3.2,"jet_phi13"        ,"Events/(0.1)"    ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_d23"     ,selection,"",tag,32,   0, 3.2,"jet_phi23"        ,"Events/(0.1)"    ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_dij_max" ,selection,"",tag,32,   0, 3.2,"jet_phi_ij_max"   ,"Events/(0.1)"    ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_dij_min" ,selection,"",tag,32,   0, 3.2,"jet_phi_ij_min"   ,"Events/(0.1)"    ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_dij_min2",selection,"",tag, 2,   0,   6,"jet_phi_ij_min2"  ,"Events/(3)"      ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_phi_dij_min3",selection,"",tag, 4,   0,   4,"jet_phi_ij_min3"  ,"Events/(1)"      ,0, CR );
        #self.construct_plot(Nj,"Pt2dPt1"                 ,selection,"",tag,20,   0,   1,"Pt2/Pt1"          ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"Pt3dPt1"                 ,selection,"",tag,20,   0,   1,"Pt3/Pt1"          ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"Pt2dPt1-Pt3dPt1"         ,selection,"",tag,20,   0,   1,"(Pt2-Pt3)/Pt1"    ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jetAK8puppi_dR12"        ,selection,"",tag,20,   0,   5,"jet_dR12"         ,"Events/(0.25)"   ,0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_dR13"        ,selection,"",tag,20,   0,   5,"jet_dR13"         ,"Events/(0.25)"   ,0, CR );
        #self.construct_plot(Nj,"nbtag"                   ,selection,"",tag, 5,-0.5, 4.5,"number of b-jets" ,"Events"          ,0 , CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi"      ,selection,"",tag,20,   0,   1,"#tau_{21}"        ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi_2"    ,selection,"",tag,20,   0,   1,"#tau_{21_2}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi_3"    ,selection,"",tag,20,   0,   1,"#tau_{21_3}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi"      ,selection,"",tag,20,   0,   1,"#tau_{31}"        ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi_2"    ,selection,"",tag,20,   0,   1,"#tau_{31_2}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau3tau1_puppi_3"    ,selection,"",tag,20,   0,   1,"#tau_{31_3}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi"      ,selection,"",tag,20,   0,   1,"#tau_{41}"        ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi_2"    ,selection,"",tag,20,   0,   1,"#tau_{41_2}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau1_puppi_3"    ,selection,"",tag,20,   0,   1,"#tau_{41_3}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi"      ,selection,"",tag,20,   0,   1,"#tau_{42}"        ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi_2"    ,selection,"",tag,20,   0,   1,"#tau_{42_2}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tau4tau2_puppi_3"    ,selection,"",tag,20,   0,   1,"#tau_{42_3}"      ,"Events/(0.05)"   ,0 , CR);
        #self.construct_plot(Nj,"jet_tauitau1_puppi"      ,selection,"",tag,20,   0,   1,"(#tau_{21(j1)}#tau_{31(j1)}#tau_{41(j1)})^{1/3}","Events/(0.05)",0 , CR);
        #self.construct_plot(Nj,"jet_tauitau1_puppi_2"    ,selection,"",tag,20,   0,   1,"(#tau_{21(j2)}#tau_{31(j2)}#tau_{41(j2)})^{1/3}","Events/(0.05)",0 , CR);
        #self.construct_plot(Nj,"jet_tauitau1_puppi_3"    ,selection,"",tag,20,   0,   1,"(#tau_{21(j3)}#tau_{31(j3)}#tau_{41(j3)})^{1/3}","Events/(0.05)",0 , CR);

    def DefineSelection_1lep(self):#==========[ 1lep REGIONS & SELECTION DEFINITION ]==========================================
        REGION=options.REGION; Nb_cut=self.Nb_cut; ST_low=options.ST_low; ST_high=options.ST_high; Nak4_cut=self.Nak4_cut; lpt_cut=self.lpt_cut; MET_cut=self.MET_cut; tau_cut=self.tau_cut; and_or=self.and_or;  Nj8="num_ak8jets"; Nj4="num_ak4jets";
        mass_low=options.mass_low; mass_high=options.mass_high; Dimass_low=options.Dimass_low; Dimass_high=options.Dimass_high; Dimass_low=options.Dimass_low; Dimass_high=options.Dimass_high;
        Filters ="  passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0  &&  "; ###+ "mtVlepnew<80 && tau41_Mj_minc<0.25 && tau21_Mj_minc<0.35 && delPhilepmet<0.25 && num_ak4jetsex<4&&" # cust for CR(VV)
        common ="  Mj_maxc>40  &&  abs(jetAK8puppi_eta)<2.4  &&  jet_pt_puppi>200  && abs(l_eta)<2.4  &&  l_pt>%s  &&  MET_et>%s  &&  W_pt>200  &&  (num_ak8jets==1||(abs(jetAK8puppi_eta_2)<2.4 && jet_pt_puppi_2>100)) &&  %s<ST&&ST<%s &&  %s<Mj_maxc&&Mj_maxc<%s &&"%(lpt_cut,MET_cut,  ST_low,ST_high,  mass_low,mass_high);
        SR1  = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv &&  70<Mj_maxc&&Mj_maxc<=100                 && tau21_Mj_maxc<0.5 ";
        SR2  = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc<0.5 ";
        SR3  = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv && 200<Mj_maxc                              && tau42_Mj_maxc<0.5 ";
        SR4  = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 && 60<Mj_minc    && tau21_Mj_maxc<0.5 && tau21_Mj_minc<0.5 ";
        SR5  = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 &&    Mj_minc<60 && tau21_Mj_maxc<0.5 ";
        SR6  = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc<0.5 ";
        SR7  = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj && 200<Mj_maxc                              && tau42_Mj_maxc<0.5 ";
        #-----------------------------------------------------------------------------------------------------------------------------------------------------
        CR1W = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv &&  70<Mj_maxc&&Mj_maxc<=100                 && tau21_Mj_maxc>0.5 ";
        CR2W = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc>0.5 ";
        CR3W = Filters + common + Nj8 +"==1 && nbtag==0 && "+Nj4+"<=4 && 1000<m_jlv && 200<Mj_maxc                              && tau42_Mj_maxc>0.5 ";
        CR4W = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 && 60<Mj_minc    && (tau21_Mj_maxc>0.5||tau21_Mj_minc>0.5) ";
        CR5W = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 &&    Mj_minc<60 && tau21_Mj_maxc>0.5 ";
        CR6W = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc>0.5 ";
        CR7W = Filters + common + Nj8 +"==2 && nbtag==0 && "+Nj4+"<=4 && 1000<m_lvj && 200<Mj_maxc                              && tau42_Mj_maxc>0.5 ";
        CR1t = Filters + common + Nj8 +"==1 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_jlv &&  70<Mj_maxc&&Mj_maxc<=100                 && tau21_Mj_maxc<0.5 ";
        CR2t = Filters + common + Nj8 +"==1 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_jlv && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc<0.5 ";
        CR3t = Filters + common + Nj8 +"==1 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_jlv && 200<Mj_maxc                              && tau42_Mj_maxc<0.5 ";
        CR4t = Filters + common + Nj8 +"==2 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 && 60<Mj_minc    && tau21_Mj_maxc<0.5 && tau21_Mj_minc<0.5 ";
        CR5t = Filters + common + Nj8 +"==2 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_lvj &&  60<Mj_maxc&&Mj_maxc<=100 &&    Mj_minc<60 && tau21_Mj_maxc<0.5 ";
        CR6t = Filters + common + Nj8 +"==2 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_lvj && 100<Mj_maxc&&Mj_maxc<=200                 && tau42_Mj_maxc<0.5 ";
        CR7t = Filters + common + Nj8 +"==2 && nbtag>=1 && "+Nj4+"<=6 && 1000<m_lvj && 200<Mj_maxc                              && tau42_Mj_maxc<0.5 ";
        if options.channel=="mu":
            #self.Make_Controlplots_for_1lep(eval(REGION),"",self.nak8jet1);                                   
            self.tdrStyle.SetErrorX(0.5); 
            self.Make_Comparison(eval("S"+REGION[1:3]),eval("C"+REGION[1:4]),"CR"+REGION[3:4],self.nak8jet1); 


    #____________[ FUNCTIONS FOR COMPARION ]_______________________
    def Make_Comparison(self,selection_SR,selection_CR,tag,numak8jets):
        if options.REGION[1:3] in ["R1","R2","R3"]: self.construct_plot( numak8jets,"m_jlv",selection_SR,selection_CR,tag,10,1000,3000,"M(JW) (GeV)" ,"Events/(200 GeV)",0 );
        if options.REGION[1:3] in ["R4","R5"]:      self.construct_plot( numak8jets,"m_lvj",selection_SR,selection_CR,tag,10,1000,3000,"M(JJW) (GeV)","Events/(200 GeV)",0 );
        if options.REGION[1:3] in ["R6","R7"]:      self.construct_plot( numak8jets,"m_jlv",selection_SR,selection_CR,tag,12, 600,3000,"M(JW) (GeV)" ,"Events/(200 GeV)",0 );
        self.construct_plot( numak8jets,"W_pt" ,selection_SR,selection_CR,tag,8, 300,1100,"PT(W) (GeV)" ,"Events",0 );
        self.construct_plot( numak8jets,"W_eta",selection_SR,selection_CR,tag,12, 0,2.4,"eta(W) (GeV)" ,"Events",0 );
        self.construct_plot( numak8jets,"Mj_maxc"        ,selection_SR,selection_CR,tag,24, 0, 120,"Mjmax (GeV)" ,"Events",0 );
        self.construct_plot( numak8jets,"Mj_minc"        ,selection_SR,selection_CR,tag,20, 0, 100,"Mjmin (GeV)" ,"Events",0 );
        self.construct_plot( numak8jets,"jet_pt_puppi"  ,selection_SR,selection_CR,tag,7,400,1100,"PT(j1) (GeV)" ,"Events/(100 GeV)",0 );
        self.construct_plot( numak8jets,"jet_pt_puppi_2",selection_SR,selection_CR,tag,7,400,1100,"PT(j2) (GeV)" ,"Events/(100 GeV)",0 );
        self.construct_plot( numak8jets,"jetAK8puppi_eta"  ,selection_SR,selection_CR,tag,10, 0., 2.4,"|eta|j1"  ,"Events",0 );
        self.construct_plot( numak8jets,"jetAK8puppi_eta_2",selection_SR,selection_CR,tag,10, 0., 2.4,"|eta|j2","Events",0 );
        self.construct_plot( numak8jets,"ST"               ,selection_SR,selection_CR,tag,12, 0.,2400,"ST (GeV)","Events",0 );


    def Make_Controlplots_for_1lep( self , selection , tag , Nj , CR=0 ):   # numak8jets,CR=0) :
        REGION= options.REGION;
        if REGION[1:3] in ["R1","R2","R3",     "S1"     ]:Nj=1;
        if REGION[1:3] in ["R4","R5","R6","R7","S2","S0"]:Nj=2;
        bins=10;low=1000; high=2000;
        if  "PS" in REGION : bins=24; low=1000; high=3400;
        if (options.ST_high<=1000) and (REGION[1:3] in ["R6"])                     : bins=13;  low= 500;  high=2000;
        if (options.ST_high<=1000) and (REGION[1:3] in ["R7"])                     : bins= 3;  low= 500;  high=1500;
        if (options.ST_high<=1000) and (REGION[1:3] in ["R3"])                     : bins= 4;  low=1000;  high=2000;
        if (options.ST_high<=1000) and (REGION[1:3] in ["R4"])                     : bins= 8;  low=1000;  high=2000;
        if (options.ST_high<=1000) and (REGION[1:3] in ["R1","R2","R5"])           : bins=10;  low=1000;  high=2000;
        if (options.ST_high> 1000) and (REGION[1:3] in ["R6","R7"])                : bins=10;  low= 600;  high=2600;
        if (options.ST_high> 1000) and (REGION[1:3] in ["R1","R2","R3","R4","R5"]) : bins= 8;  low=1000;  high=2600;
#        if REGION=="PS0" or REGION=="PS2" or REGION=="PS3": 
 #           self.construct_plot(MET_et,"MET",selection,"",tag,bins,low,high,"M(j,lv) (GeV)"  ,"Events / ("+str( (high-low)/bins )+" GeV)",0, CR );
  #      else:
        if REGION[1:3] in ["R1","R2","R3","R6","R7","S1"] : self.construct_plot(Nj,"m_jlv",selection,"",tag,bins,low,high,"M(j,lv) (GeV)"  ,"Events / ("+str( (high-low)/bins )+" GeV)",0, CR );
        if REGION[1:3] in ["R4","R5"               ,"S2"] : self.construct_plot(Nj,"m_lvj",selection,"",tag,bins,low,high,"M(j,j,lv) (GeV)","Events / ("+str( (high-low)/bins )+" GeV)",0, CR );
        
        #if REGION[0:] is "SR": ##use only "SR" for nj = 1 and 2 # self.nak8jet1 =1; # self.nak8jet2 =2;
        #self.construct_plot(Nj,"jetAK8puppi_eta"     ,selection,"",tag,12, 0,2.4,"|#eta(j1)|","Events/(0.2)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_eta_2"   ,selection,"",tag,12, 0,2.4,"|#eta(j2)|","Events/(0.2)", 0, CR);
        #self.construct_plot(Nj,"Mj_minc"              ,selection,"",tag,12, 0,120,"M_{j}^{min}","Events/(10 GeV)", 0, CR);
        #self.construct_plot(Nj,"Mj_maxc"              ,selection,"",tag,21,40,250,"M_{j}^{max}","Events/(10 GeV)", 0, CR);
        #self.construct_plot(Nj,"m_jlv",selection,"",tag,25,1000,3500,"M(j,lv) (GeV)"   ,"Events / 100 GeV)",0, CR );
        #self.construct_plot(Nj,"jet_mass_puppi_corr" ,selection,"",tag,44, 0,220,"M_{j1(ak8)} (GeV)","Events/(2.5 GeV)", 0, CR );
        #self.construct_plot(Nj,"jetAK8puppi_sdcorr_2",selection,"",tag,40, 0,100,"M_{j2(aK8)} (GeV)","Events/(2.5 GeV)", 0, CR );
        #self.construct_plot(Nj,"num_ak8jets"         ,selection,"",tag,4 ,-0.5, 3.5,"N_{ak8-jets}"   ,"Events" ,0, CR);
        #self.construct_plot(Nj,"MassVV[0]"           ,selection,"",tag,35,0,3500,"M_{W_{1}W_{2}} (GeV)","Events/(100 GeV)",0, CR );
        #self.construct_plot(Nj,"tau42_Mj_maxc",selection,"",tag,20,0,1,"#tau42_Mj_maxc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau31_Mj_maxc",selection,"",tag,20,0,1,"#tau31_Mj_maxc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau21_Mj_maxc",selection,"",tag,20,0,1,"#tau21_Mj_maxc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau41_Mj_maxc",selection,"",tag,20,0,1,"#tau41_Mj_maxc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau42_Mj_minc",selection,"",tag,20,0,1,"#tau42_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau21_Mj_minc",selection,"",tag,20,0,1,"#tau21_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau31_Mj_minc",selection,"",tag,20,0,1,"#tau31_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau41_Mj_minc",selection,"",tag,20,0,1,"#tau41_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau43_Mj_maxc",selection,"",tag,20,0,1,"#tau43_Mj_maxc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau32_Mj_maxc",selection,"",tag,20,0,1,"#tau32_Mj_maxc","Events/(0.05)", 0, CR);        
        #self.construct_plot(Nj,"tau43_Mj_minc",selection,"",tag,20,0,1,"#tau43_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"tau32_Mj_minc",selection,"",tag,20,0,1,"#tau32_Mj_minc","Events/(0.05)", 0, CR);    
        #self.construct_plot(Nj,"tau31_Mj_minc",selection,"",tag,20,0,1,"#tau31_Mj_minc","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"PTw2divPTw1"   ,selection,"",tag,20,0,1,"PTw2/PTw1","Events", 0, CR);
        #self.construct_plot(Nj,"PTw3divPTw1"   ,selection,"",tag,20,0,1,"PTw3/PTw1","Events", 0, CR);
        #self.construct_plot(Nj,"PTw3divPTw2"   ,selection,"",tag,20,0,1,"PTw3/PTw2","Events", 0, CR);
        #self.construct_plot(Nj,"MassVV[1]"     ,selection,"",tag,35,0,3500,"M_{W_{1}W_{3}} (GeV)","Events/(100 GeV)",0, CR );
        #self.construct_plot(Nj,"MassVV[2]"     ,selection,"",tag,35,0,3500,"M_{W_{2}W_{3}} (GeV)","Events/(100 GeV)",0, CR );        
        #self.construct_plot(Nj,"W_pt"          ,selection,"",tag,10,0, 1000,"PT(l,MET) (GeV)","Events/(100 GeV)", 0, CR);
        #self.construct_plot(Nj,"jet_pt_puppi"  ,selection,"",tag,12,0, 1200,"PT_{j1(AK8)} (GeV)","Events/(100 GeV)", 0, CR);
        #self.construct_plot(Nj,"jet_pt_puppi_2",selection,"",tag,10,0, 1000,"PT_{j2(AK8)} (GeV)","Events/(100 GeV)", 0, CR);
        #self.construct_plot(Nj,"ST"            ,selection,"",tag,12,0, 2400,"S_{T} (GeV)","Events/(200 GeV)", 0, CR);
        #self.construct_plot(Nj,"num_ak4jets"   ,selection,"",tag,9,-0.5, 8.5,"N_{AK4jets}","Events", 0, CR);
        #self.construct_plot(Nj,"num_ak4jetsex" ,selection,"",tag,9,-0.5, 8.5,"N_{AK4jets_exc}","Events", 0, CR);
        #self.construct_plot(Nj,"mtVlepnew"     ,selection,"",tag,25,0, 500 ,"MT(l,MET) (GeV)","Events/(20 GeV)", 0, CR);
        #self.construct_plot(Nj,"l_pt"          ,selection,"",tag,8,0, 800,"PT_{l} (GeV)","Events/(100 GeV)", 0, CR);
        #self.construct_plot(Nj,"trackIso"      ,selection,"",tag,20,0,2.5,"trackIso(l) (GeV)","Events/0.125", 0, CR);
        #self.construct_plot(Nj,"muisolation"   ,selection,"",tag,20,0,0.1,"muisolation","Events/0.05", 0, CR);
        #self.construct_plot(Nj,"deltaRlepjet"  ,selection,"",tag,20,3,4      ,"#DeltaR(l,j1)"        ,"Events/(0.25)" , 0, CR);
        #self.construct_plot(Nj,"deltaRlepjet_2",selection,"",tag,20,0,5      ,"#DeltaR(l,j2)"        ,"Events/(0.25)" , 0, CR);
        #self.construct_plot(Nj,"deltaRjet1jet2",selection,"",tag,12,0,6,"#DeltaR(j1,j2)","Events/(0.5)", 0, CR); 
        #self.construct_plot(Nj,"delPhijetmet"  ,selection,"",tag,5,2,3.14159,"|#Delta#phi(j1,MET)|" ,"Events", 0, CR);
        #self.construct_plot(Nj,"delPhijetmet_2",selection,"",tag,10,1,3.14159,"|#Delta#phi(j2,MET)|" ,"Events/(2x0.107)", 0, CR);
        #self.construct_plot(Nj,"delPhijetlep"  ,selection,"",tag,10,1.2,3.14159,"|#Delta#phi(j1,l)|"   ,"Events", 0, CR);
        #self.construct_plot(Nj,"delPhijetlep_2",selection,"",tag,10,1,3.14159,"|#Delta#phi(j2,l)|"   ,"Events/(2x0.107)", 0, CR);
        #self.construct_plot(Nj,"delPhilepmet"    ,selection,"",tag,5,0,1.5,"|#Delta#phi(l,MET)|","Events/(2x0.157)", 0, CR);
        #self.construct_plot(Nj,"MET_et"        ,selection,"",tag,20,0,800,"MET (GeV)","Events/(40 GeV)", 0, CR);
        #self.construct_plot(Nj,"nPV"           ,selection,"",tag,35,5,40,"nPV","Events", 0, CR);
        #self.construct_plot(Nj,"nbtag"         ,selection,"",tag,5,-0.5,4.5,"number of b-jets","Events", 0, CR);
        #self.construct_plot(Nj,"jet_tau2tau1_puppi",selection,"",tag,20,0,1,"#tau_{21_j1}","Events/(0.05)", 0, CR);
        ##self.construct_plot(Nj,"jet_tau2tau1_puppi_2",selection,"",tag,20,0,1,"#tau_{21_j2}","Events/(0.05)", 0, CR);
        ##self.construct_plot(Nj,"jet_tau4tau2_puppi",selection,"",tag,20,0,1,"#tau_{42_j1}","Events/(0.05)", 0, CR);
        ##self.construct_plot(Nj,"jet_tau4tau2_puppi_2",selection,"",tag,20,0,1,"#tau_{42_j2}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau31",selection,"",tag,20,0,1,"#tau_{31_j1}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau32",selection,"",tag,20,0,1,"#tau_{32_j1}","Events/(0.05)", 0, CR);
        ##self.construct_plot(Nj,"jetAK8puppi_tau41",selection,"",tag,20,0,1,"#tau_{41_j1}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau43",selection,"",tag,20,0,1,"#tau_{43_j1}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau31_2",selection,"",tag,20,0,1,"#tau_{31_j2}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau32_2",selection,"",tag,20,0,1,"#tau_{32_j2}","Events/(0.05)", 0, CR);
        ##self.construct_plot(Nj,"jetAK8puppi_tau41_2",selection,"",tag,20,0,1,"#tau_{41_j2}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"jetAK8puppi_tau43_2",selection,"",tag,20,0,1,"#tau_{43_j2}","Events/(0.05)", 0, CR);
        #self.construct_plot(Nj,"l_eta"           ,selection,"",tag,12,0,2.4,"|#eta(l)|"          ,"Events/(0.2)", 0, CR);
        #self.construct_plot(Nj,"deltaetajet1jet2",selection,"",tag,10,0,3  ,"|#Delta#eta(j1,j2)|","Events"   , 0, CR);
        #self.construct_plot(Nj,"deltaetajet1lep" ,selection,"",tag,10,0,3  ,"|#Delta#eta(j1,l)|" ,"Events"   , 0, CR);
        #self.construct_plot(Nj,"deltaetajet2lep" ,selection,"",tag,10,0,3  ,"|#Delta#eta(j2,l)|" ,"Events"   , 0, CR);
        #self.construct_plot_productof3tau(Nj,"jet_tau2tau1_puppi","jetAK8puppi_tau31","jetAK8puppi_tau41",selection,"",tag,20,0,1,"#(tau_{21_j1}*#tau_{31_j1}*#tau_{41_j1})^{1/3}","Events/(0.05)", 0, CR);
        #self.construct_plot_averageof3tau(Nj,"jet_tau2tau1_puppi","jetAK8puppi_tau31","jetAK8puppi_tau41",selection,"",tag,20,0,1,"(#tau_{21_j1}+#tau_{31_j1}+#tau_{41_j1})/3","Events/(0.05)", 0, CR);
        ##self.construct_plot_productof3tau(Nj,"W_pt","jet_pt_puppi","jet_pt_puppi_2",selection,"",tag,30,0,1500,"(PT_{W_{l}}*PT_{j1}*PT_{j2})^{1/3} [GeV]","Events/(50 GeV)", 0, CR);
        ##self.construct_plot_productof2Wpt(Nj,"W_pt","jet_pt_puppi",selection,"",tag,30,0,1800,"(PT_{W_{l}}*PT_{j1})^{1/2} [GeV]","Events/(60 GeV)", 0, CR);
        ##self.construct_plot_averageof3tau(Nj,"MassVV[0]","MassVV[1]","MassVV[2]",selection,"",tag,35,0,3500,"(M_{W_{1}W_{2}}+M_{W_{1}W_{3}}+M_{W_{2}W_{3}})/3","Events/(100 GeV)", 0, CR);
        ##self.construct_plot_averageof2MassVV(Nj,"MassVV[1]","MassVV[2]",selection,"",tag,35,0,3500,"(M_{W_{1}W_{3}}+M_{W_{2}W_{3}})/2","Events/(100 GeV)", 0, CR);

    #def Make_Controlplots_consider_both_jets(self,selection,tag,Nj,CR=0) :
    #    self.construct_plot(Nj,"tau21_Mj_maxc",selection,"",tag,20,0,1,"#tau_{21} for all jets","Events/(0.05)", 0, CR);         # numak8jets,CR=0) :
        #if tag[10:16] in ["tau_40","tau_70"]:
        #if tag[10:17] in ["tau_100","tau_200"]:
        #self.make_controlplot(Nj,"tau42_Mj_maxc",selection,"",tag,20,0,1,"#tau_{42}","Events/(0.05)", 0, CR);




#=========================================================================================
    def construct_plot(self,Nj,variable,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="",logy=1,CR=0):
        channel = options.channel; REGION  = options.REGION; print "--> Selection ==> {"+cut+"}";
        if channel=="had": path = "/eos/cms/store/user/chench/NEW/deepv2/"      +"mu";         lumi = 35.52;
        if channel=="mu" : path = "/eos/cms/store/user/xulyu/data2016/2016/"+ self.channel;lumi = 35.9;   #"/eos/cms/store/user/xulyu/data2016/sj/tree16/"  for new ntuples with ak4_exc 
       #====== DEFINE CANVAS Ccommon for 0lep and 1lep ==========================
        canvas_controlplot = TCanvas(self.channel+"_"+REGION+"_"+variable,self.channel+"_"+REGION+"_"+variable, 700,700);
        if channel=="had": STcut=1250;
        if channel=="mu" : STcut=1000;
        if (REGION[0:2] in ["",""]) and (options.ST_high>STcut): # "SR", "PS"  was -->"",""
            fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
            fPads1.SetTicky(1); fPads1.SetTickx(1)  ##prefit_SYS->GetYaxis()->SetTickLength(0.02);
            fPads1.SetBottomMargin(0.12);  fPads1.SetLeftMargin(0.11 );  fPads1.SetRightMargin(0.03); fPads1.Draw(); fPads1.cd();
        if (REGION[0:2] in ["CR","VR","PS","SR"]) or (options.ST_high<=STcut) or tag[10:14] == "tau_":  #  "" was --> "PS"
            fPads1 = TPad("pad1", "Run2", 0.0, 0.29, 1.00, 1.00);
            fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.29);
            fPads1.SetBottomMargin( 0.007);  fPads1.SetLeftMargin(   0.10);  fPads1.SetRightMargin(  0.03);
            fPads2.SetLeftMargin(   0.10 );  fPads2.SetRightMargin(  0.03);  fPads2.SetBottomMargin( 0.25);
            fPads1.Draw(); fPads2.Draw();  fPads1.cd();
        #====================== DEFINE TREES AND HISTOS ======================================
        if channel=="had":
            t_data    = TChain("PKUTree");
            for DataPeriod in ["B","C","D","E","F","G","H2","H3"]:t_data.Add(path + "_out_"+DataPeriod+".root"); 
            h_data=TH1D("h_data","h_DATA"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_data.Sumw2();  c_data   = h_data.Clone(   "c_data"   );
            if REGION in ["SR4","SR5","CR4q","CR5q","CR4t","CR5t"]:
               t_Signal1 = TChain("PKUTree");  t_Signal1.Add(path+"_out_M3000-R0-4.root"  );  h_Signal1=TH1D("h_Signal1","h_Signal1"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal1.Sumw2();
               t_Signal2 = TChain("PKUTree");  t_Signal2.Add(path+"_out_M3000-R0-5.root"  );  h_Signal2=TH1D("h_Signal2","h_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal2.Sumw2();
               t_Signal3 = TChain("PKUTree");  t_Signal3.Add(path+"_out_M3000-R0-6.root"  );  h_Signal3=TH1D("h_Signal3","h_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal3.Sumw2();
               t_Signal4 = TChain("PKUTree");  t_Signal4.Add(path+"_out_M3500-R0-4.root"  );  h_Signal4=TH1D("h_Signal4","h_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal4.Sumw2();
               t_Signal5 = TChain("PKUTree");  t_Signal5.Add(path+"_out_M3500-R0-5.root"  );  h_Signal5=TH1D("h_Signal5","h_Signal5"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal5.Sumw2();
               t_Signal6 = TChain("PKUTree");  t_Signal6.Add(path+"_out_M3500-R0-6.root"  );  h_Signal6=TH1D("h_Signal6","h_Signal6"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal6.Sumw2();
            else:
               t_Signal1 = TChain("PKUTree");  t_Signal1.Add(path+"_out_M3000-R0-06.root"  );  h_Signal1=TH1D("h_Signal1","h_Signal1"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal1.Sumw2();
               t_Signal2 = TChain("PKUTree");  t_Signal2.Add(path+"_out_M3000-R0-08.root"  );  h_Signal2=TH1D("h_Signal2","h_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal2.Sumw2();
               t_Signal3 = TChain("PKUTree");  t_Signal3.Add(path+"_out_M3000-R0-1.root"  );  h_Signal3=TH1D("h_Signal3","h_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal3.Sumw2();
               t_Signal4 = TChain("PKUTree");  t_Signal4.Add(path+"_out_M3500-R0-06.root"  );  h_Signal4=TH1D("h_Signal4","h_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal4.Sumw2();
               t_Signal5 = TChain("PKUTree");  t_Signal5.Add(path+"_out_M3500-R0-08.root"  );  h_Signal5=TH1D("h_Signal5","h_Signal5"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal5.Sumw2();
               t_Signal6 = TChain("PKUTree");  t_Signal6.Add(path+"_out_M3500-R0-1.root"  );  h_Signal6=TH1D("h_Signal6","h_Signal6"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal6.Sumw2();

            t_QCD     = TChain("PKUTree");  t_QCD.Add(    path+"_PKUTree_QCD.root"  );  h_QCD    =TH1D("h_QCD"    ,"h_QCD"    +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_QCD.Sumw2();   
            t_WJets   = TChain("PKUTree");  t_WJets.Add(  path+"_PKUTree_WJets.root");  h_WJets  =TH1D("h_WJets"  ,"h_WJets"  +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_WJets.Sumw2(); 
            t_ZJets   = TChain("PKUTree");  t_ZJets.Add(  path+"_PKUTree_ZJets.root");  h_ZJets  =TH1D("h_ZJets"  ,"h_ZJets"  +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_ZJets.Sumw2(); 
            t_TTbar   = TChain("PKUTree");  t_TTbar.Add(  path+"_PKUTree_TT.root"   );  h_TTbar  =TH1D("h_TTbar"  ,"h_TTbar"  +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_TTbar.Sumw2(); 
            t_STop    = TChain("PKUTree");  t_STop.Add(   path+"_PKUTree_ST.root"   );  h_STop   =TH1D("h_STop"   ,"h_STop"   +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_STop.Sumw2();  
            t_VV      = TChain("PKUTree");  t_VV.Add(     path+"_PKUTree_VV.root"   );  h_VV     =TH1D("h_VV"     ,"h_VV"     +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_VV.Sumw2();     
            h_TotalMC     =TH1D("h_TotalMC"     ,"h_TotalMC"     +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_TotalMC.Sumw2();
           
            h_TotalMC.Add(h_QCD); h_TotalMC.Add(h_WJets);  h_TotalMC.Add(h_ZJets);  h_TotalMC.Add(h_TTbar); h_TotalMC.Add(h_STop); h_TotalMC.Add(h_VV);

            if tag=="CRq" or REGION=="CR1q" or REGION=="CR2q" or REGION=="CR3q" or REGION=="CR4q" or REGION=="CR5q" or REGION=="CR2t" or REGION=="SR1" or REGION=="SR2" or REGION=="SR3" or REGION=="SR4" or REGION=="SR5":
                c_data=h_data.Clone("c_data");  c_Signal1=h_Signal1.Clone("c_Signal1");  c_Signal2=h_Signal2.Clone("c_Signal2");  c_Signal3=h_Signal3.Clone("c_Signal3");  c_Signal4=h_Signal4.Clone("c_Signal4");  
                c_QCD=h_QCD.Clone("c_QCD");  c_WJets=h_WJets.Clone("c_WJets");  c_ZJets=h_ZJets.Clone("c_ZJets");  c_TTbar=h_TTbar.Clone("c_TTbar");  c_STop=h_STop.Clone("c_STop");  c_VV=h_VV.Clone("c_VV");c_TotalMC=h_TotalMC.Clone("c_TotalMC");
          
        if channel=="mu":
            t_data    = TChain("PKUTree");  t_data.Add(   path+"_PKUTree_2017vvv.root"        );  h_data   =TH1D("h_data"   ,"h_DATA"   +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_data.Sumw2();
            t_Signal1 = TChain("PKUTree");  t_Signal1.Add(path+"_out_case1_off_3000.root"     );  h_Signal1=TH1D("h_Signal1","h_Signal1"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Signal1.Sumw2();
            t_Signal2 = TChain("PKUTree");  t_Signal2.Add(path+"_out_case1_off_1500.root"     );  h_Signal2=TH1D("h_Signal2","h_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Signal2.Sumw2();
            t_Signal3 = TChain("PKUTree");  t_Signal3.Add(path+"_out_case2_off_3000.root"     );  h_Signal3=TH1D("h_Signal3","h_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Signal3.Sumw2();
            t_Signal4 = TChain("PKUTree");  t_Signal4.Add(path+"_out_case2_off_1500.root"     );  h_Signal4=TH1D("h_Signal4","h_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Signal4.Sumw2();
            t_WJets   = TChain("PKUTree");  t_WJets.Add(  path+"_PKUTree_WJetsPt180_xww.root" );  h_WJets  =TH1D("h_WJets"  ,"h_WJets"  +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_WJets.Sumw2();
            t_TTbar   = TChain("PKUTree");  t_TTbar.Add(  path+"_PKUTree_TTBARpowheg_xww.root");  h_TTbar  =TH1D("h_TTbar"  ,"h_TTbar"  +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TTbar.Sumw2();
            t_STop    = TChain("PKUTree");  t_STop.Add (  path+"_PKUTree_SingleTop_xww.root"  );  h_STop   =TH1D("h_STop"   ,"h_STop"   +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_STop.Sumw2();
            t_VV      = TChain("PKUTree");  t_VV.Add   (  path+"_PKUTree_VV_xww.root"         );  h_VV     =TH1D("h_VV"     ,"h_VV"     +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_VV.Sumw2();
            if tag in ["CRW","CRt"]: 
                c_data=h_data.Clone("c_data"); c_Signal1=h_Signal1.Clone("c_Signal1"); c_Signal2=h_Signal2.Clone("c_Signal2"); c_Signal3=h_Signal3.Clone("c_Signal3"); c_Signal4=h_Signal4.Clone("c_Signal4"); c_WJets=h_WJets.Clone("c_WJets"); c_TTbar=h_TTbar.Clone("c_TTbar"); c_STop=h_STop.Clone("c_STop"); c_VV=h_VV.Clone("c_VV");
        hstack_TotalMC= THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle));                #cstack_TotalMC= THStack("cstack_TotalMC","cstack_TotalMC"+";%s;%s"%(xtitle,ytitle));

        #=================== SET WEIGHTS, SCALE TREES, DEFINE TOTAL AND STACK  =================================================
        if "SR" or "CR" in REGION : options.scale1=10; options.scale2=30;#100000,300000
        ext1="";ext2="";ext3=""; weight="weight";        #        lumi=1;
        if channel=="had":         # if options.loop==1:  ext1=" && %s<jet_mass_puppi_max&&jet_mass_puppi_max<%s "%(options.mass_low,options.mass_high); ext2=" && %s<jet_mass_puppi_mid&&jet_mass_puppi_mid<%s "%(options.mass_low,options.mass_high); ext3=" && %s<jet_mass_puppi_min&&jet_mass_puppi_min<%s "%(options.mass_low,options.mass_high);
            for region in ["h_"]:
                for sample in ["data","QCD","WJets","ZJets","TTbar","STop","VV","Signal1","Signal2","Signal3","Signal4","Signal5","Signal6"] :
                    if sample in ["data"] and "SR" not in REGION            :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample),                                               cut+ext1  );   # No weights on data
                    if sample in ["Signal1","Signal2","Signal3"]                      :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale1,cut+ext1) );   # Extra scaling on signal
                    if sample in ["Signal4","Signal5","Signal6"]                      :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale2,cut+ext1) );   # Extra scaling on signal
                    if sample in ["QCD","WJets","ZJets","TTbar","STop","VV"]:eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,1             ,cut+ext1) );   # MC only weights on MC BKG
                    #if options.loop==1:
                    #    variable2=variable[0:len(variable)-3]+"mid"; print "-->in *loop* with :",variable2,"  with ",ext2;
                    #    tmp =TH1D("tmp","tmp"+";%s;%s"%(xtitle,ytitle),nbin,min,max); tmp.Reset(); tmp.Sumw2();
                    #    if sample in ["data"] and "SR" not in REGION            :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2),                                            cut+ext2  );   # No weights on data
                    #    if sample in ["Signal1","Signal3"]                      :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale1,cut+ext2) );   # Extra scaling on signal
                    #    if sample in ["Signal2","Signal4"]                      :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale2,cut+ext2) );   # Extra scaling on signal
                    #    if sample in ["QCD","WJets","ZJets","TTbar","STop","VV"]:eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,1             ,cut+ext2) );   # MC only weights on MC BKGs
                    #    eval("h_"+sample).Add(tmp); tmp.Delete();
                    #    variable3=variable[0:len(variable)-3]+"min"; print "-->in *loop* with :",variable3,"  with ",ext3;
                    #    tmp =TH1D("tmp","tmp"+";%s;%s"%(xtitle,ytitle),nbin,min,max); tmp.Reset(); tmp.Sumw2();
                    #    if sample in ["data"] and "SR" not in REGION            :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable3),                                            cut+ext3  );   # No weights on data
                    #    if sample in ["Signal1","Signal3"]                      :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable3), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale1,cut+ext3) );   # Extra scaling on signal
                    #    if sample in ["Signal2","Signal4"]                      :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable3), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale2,cut+ext3) );   # Extra scaling on signal
                    #    if sample in ["QCD","WJets","ZJets","TTbar","STop","VV"]:eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable3), "(%s*%s*%s)*(%s)"%(weight,lumi,1             ,cut+ext3) );   # MC only weights on MC BKGs
                    #    eval("h_"+sample).Add(tmp); tmp.Delete();                       
            if REGION=="CR4q":h_QCDWjets_CR4 = h_data.Clone("h_QCDWjets_CR4"); h_QCDWjets_CR4.Add(h_ZJets,-1); h_QCDWjets_CR4.Add(h_TTbar,-1); h_QCDWjets_CR4.Add(h_STop,-1); h_QCDWjets_CR4.Add(h_VV,-1);

            hstack_TotalMC.Add(h_VV);hstack_TotalMC.Add(h_STop);hstack_TotalMC.Add(h_TTbar);hstack_TotalMC.Add(h_ZJets);hstack_TotalMC.Add(h_WJets);hstack_TotalMC.Add(h_QCD);
            h_TotalMC = TH1D("h_TotalMC","h_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);h_TotalMC.Sumw2(); h_TotalMC.Add(h_QCD);h_TotalMC.Add(h_WJets);h_TotalMC.Add(h_ZJets);h_TotalMC.Add(h_TTbar);h_TotalMC.Add(h_STop);h_TotalMC.Add(h_VV);
            h_data.SetBinErrorOption(TH1D.kPoisson); 
            if REGION=="CR1q" or REGION=="CR2q" or REGION=="CR3q" or REGION=="CR4q" or REGION=="CR5q" or REGION=="CR2t" or REGION=="SR1" or REGION=="SR2" or REGION=="SR3" or REGION=="SR4" or REGION=="SR5":                #cstack_TotalMC.Add(c_VV);hstack_TotalMC.Add(c_STop);cstack_TotalMC.Add(c_TTbar);cstack_TotalMC.Add(c_ZJets);cstack_TotalMC.Add(c_WJets);cstack_TotalMC.Add(c_QCD);
                c_TotalMC = TH1D("c_TotalMC","c_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);c_TotalMC.Sumw2(); c_TotalMC.Add(c_QCD);c_TotalMC.Add(c_WJets);c_TotalMC.Add(c_ZJets);c_TotalMC.Add(c_TTbar);c_TotalMC.Add(c_STop);c_TotalMC.Add(c_VV);
                c_Rest    = TH1D("c_Rest","c_Rest"+";%s;%s"%(xtitle,ytitle),nbin,min,max);c_Rest.Sumw2(); c_Rest.Add(c_WJets);c_Rest.Add(c_ZJets);c_Rest.Add(c_TTbar);c_Rest.Add(c_STop);c_Rest.Add(c_VV);
                c_data.SetBinErrorOption(TH1D.kPoisson);
            if REGION=="PS0" or REGION=="PS2" or REGION=="PS3":
                for h in [h_data,h_Signal1,h_Signal2,h_Signal3,h_Signal4,h_Signal5,h_Signal6,h_QCD,h_WJets,h_ZJets,h_TTbar,h_STop,h_VV,h_TotalMC ]:h=UnderOverFlow1D(h);   # UNDEROVERFLOWS for All histos used after reweighting
            else:  
                for h in [h_data,h_Signal1,h_Signal2,h_Signal3,h_Signal4,h_Signal5,h_Signal6,h_QCD,h_WJets,h_ZJets,h_TTbar,h_STop,h_VV,h_TotalMC,    c_data,c_Signal1,c_Signal2,c_Signal3,c_Signal4,c_QCD,c_WJets,c_ZJets,c_TTbar,c_STop,c_VV,c_TotalMC, c_Rest ]:h=UnderOverFlow1D(h);   # UNDEROVERFLOWS for All histos used after reweighting

        if channel=="mu":                   # if options.loop==1:   #ext1=" && %s<Mj_maxc&&Mj_maxc<%s "%(options.mass_low,options.mass_high); ext2=" && %s<Mj_minc&&Mj_minc<%s "%(options.mass_low,options.mass_high);            
            for region in ["h_","c_"]:      # h_ for Signal Region c_ for Control Region
                if region=="c_"             : cut=cut1;
                if region=="c_" and tag=="" : break;
                print "loop for ",region,"histos, tag =",tag
                for sample in ["data","WJets","TTbar","STop","VV","Signal1","Signal2","Signal3","Signal4"] :   #tmp=h_data.Clone("tmp");tmp.Reset();
                    if sample=="data" and "SR" not in REGION:eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample),                                               cut+ext1  );   # No weights on data                
                    if sample in ["Signal1","Signal3"]          :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale1,cut+ext1) );   # Extra scaling on signal
                    if sample in ["Signal2","Signal4"]          :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale2,cut+ext1) );   # Extra scaling on signal
                    if sample in ["WJets","TTbar","STop","VV"]  :eval("t_"+sample).Draw("abs(%s) >> %s%s"%(variable,region,sample), "(%s*%s*%s)*(%s)"%(weight,lumi,1             ,cut+ext1) );   # MC only weights on MC BKG
                    #if options.loop==1:
                    #    variable2=variable[0:len(variable)-3]+"min"; print "-->in *loop* with :",variable2,"  with ",ext2;
                    #    tmp =TH1D("tmp","tmp"+";%s;%s"%(xtitle,ytitle),nbin,min,max); tmp.Reset(); tmp.Sumw2();
                    #    if sample in ["data"] and "SR" not in REGION:eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2),                                            cut+ext2  );   # No weights on data
                    #    if sample in ["Signal1","Signal3"]          :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale1,cut+ext2) );   # Extra scaling on signal
                    #    if sample in ["Signal2","Signal4"]          :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,options.scale2,cut+ext2) );   # Extra scaling on signal
                    #    if sample in ["WJets","TTbar","STop","VV"]  :eval("t_"+sample).Draw("abs(%s) >> tmp"%(variable2), "(%s*%s*%s)*(%s)"%(weight,lumi,1             ,cut+ext2) );   # MC only weights on MC BKGs
                    #    eval("h_"+sample).Add(tmp); tmp.Delete();
            hstack_TotalMC.Add(h_VV);hstack_TotalMC.Add(h_STop);hstack_TotalMC.Add(h_TTbar);hstack_TotalMC.Add(h_WJets);
            h_TotalMC = TH1D("h_TotalMC","h_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);h_TotalMC.Sumw2();h_TotalMC.Add(h_WJets);h_TotalMC.Add(h_TTbar);h_TotalMC.Add(h_STop);h_TotalMC.Add(h_VV);
            h_data.SetBinErrorOption(TH1D.kPoisson); 
            for h in [h_data,h_Signal1,h_Signal2,h_Signal3,h_Signal4,h_Signal5,h_Signal6,h_WJets,h_TTbar,h_STop,h_VV,h_TotalMC,]:h=UnderOverFlow1D(h);   # UNDEROVERFLOWS for All histos used after reweighting
            if tag in ["CRW","CRt"]:                #cstack_TotalMC.Add(c_VV); cstack_TotalMC.Add(c_STop); cstack_TotalMC.Add(c_TTbar); cstack_TotalMC.Add(c_WJets);
                c_TotalMC = TH1D("c_TotalMC","c_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max); c_TotalMC.Sumw2(); c_TotalMC.Add(c_WJets); c_TotalMC.Add(c_TTbar); c_TotalMC.Add(c_STop); c_TotalMC.Add(c_VV);
                if tag=="CRW": 
                    c_Rest= TH1D("c_Rest","c_Rest"+";%s;%s"%(xtitle,ytitle),nbin,min,max); c_Rest.Sumw2(); c_Rest.Add(c_TTbar); c_Rest.Add(c_STop); c_Rest.Add(c_VV); c_Rest=UnderOverFlow1D(c_Rest);
                if tag=="CRt": 
                    c_Rest= TH1D("c_Rest","c_Rest"+";%s;%s"%(xtitle,ytitle),nbin,min,max); c_Rest.Sumw2(); c_Rest.Add(c_WJets); c_Rest.Add(c_VV);   c_Rest=UnderOverFlow1D(c_Rest);
                    c_top = TH1D("c_top" ,"c_top" +";%s;%s"%(xtitle,ytitle),nbin,min,max); c_top.Sumw2();  c_top.Add(c_TTbar);  c_top.Add(c_STop);  c_top =UnderOverFlow1D(c_top);
                    h_top = TH1D("h_top" ,"h_top" +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_top.Sumw2();  h_top.Add(h_TTbar);  h_top.Add(h_STop);  h_top =UnderOverFlow1D(h_top);
                c_data.SetBinErrorOption(TH1D.kPoisson);
                for c in [c_data,c_Signal1,c_Signal2,c_Signal3,c_Signal4,c_WJets,c_TTbar,c_STop,c_VV,c_TotalMC,]: c=UnderOverFlow1D(c); #print c.Integral();   # UNDEROVERFLOWS for All histos used after reweighting

        #============ COSMETICS ====================================================================
        h_data.SetLineColor(self.color_palet["data"]); h_data.SetFillColor(self.color_palet["data"]);
        h_Signal1.SetLineColor(1); h_Signal1.SetFillStyle(0); h_Signal1.SetLineWidth(2); h_Signal1.SetLineStyle(2);
        h_Signal2.SetLineColor(4); h_Signal2.SetFillStyle(0); h_Signal2.SetLineWidth(2); h_Signal2.SetLineStyle(2);
        h_Signal3.SetLineColor(kGreen+2);h_Signal3.SetFillStyle(0); h_Signal3.SetLineWidth(2); h_Signal3.SetLineStyle(2);
        h_Signal4.SetLineColor(93);h_Signal4.SetFillStyle(0); h_Signal4.SetLineWidth(2); h_Signal4.SetLineStyle(2);
        h_Signal5.SetLineColor(8); h_Signal5.SetFillStyle(0); h_Signal5.SetLineWidth(2); h_Signal5.SetLineStyle(2);
        h_Signal6.SetLineColor(17); h_Signal6.SetFillStyle(0); h_Signal6.SetLineWidth(2); h_Signal6.SetLineStyle(2);
        h_WJets.SetLineColor(1);   h_WJets.SetFillColor(self.color_palet["WJets"]); h_WJets.SetLineWidth(0);
        h_TTbar.SetLineColor(1);   h_TTbar.SetFillColor(self.color_palet["TTbar"]); h_TTbar.SetLineWidth(0);
        h_STop.SetLineColor(1);    h_STop.SetFillColor( self.color_palet["STop"] ); h_STop.SetLineWidth(0);
        h_VV.SetLineColor(1);      h_VV.SetFillColor(   self.color_palet["VV"]   ); h_VV.SetLineWidth(0);
        if channel=="had":
            h_QCD.SetLineColor(1);  h_QCD.SetFillColor(  self.color_palet["QCD"]  ); h_QCD.SetLineWidth(0);
            h_ZJets.SetLineColor(1);h_ZJets.SetFillColor(self.color_palet["ZJets"]); h_ZJets.SetLineWidth(0); #h_Signal1.Draw("HIST");
        h_TotalMC.SetLineStyle(3); h_TotalMC.SetMarkerStyle(0); h_TotalMC.SetLineWidth(5); h_TotalMC.SetLineColor(15); #h_TotalMC.SetFillColorAlpha(1,0.3); h_TotalMC.SetFillStyle(3344); #h_TotalMC.SetFillColor(1);

        #============ DRAW COMPARISON PLOTS ==================
        if tag in ["CRW","CRt"]:
            Leg = TLegend(0.6, 0.7, 0.76, 0.92, "", "NDC"); Leg.SetBorderSize(0); Leg.SetLineColor(0); Leg.SetFillColor(0); Leg.SetFillStyle(0); Leg.SetLineWidth(0); Leg.SetLineStyle(0); Leg.SetTextFont(42); Leg.SetTextSize(.04); Leg.SetBorderSize(0); #Leg.SetNColumns(2);
            leg = TLegend(0.13, 0.8, 0.8, 0.92, "", "NDC"); leg.SetBorderSize(0); leg.SetLineColor(0); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetLineWidth(0); leg.SetLineStyle(0); leg.SetTextFont(42); leg.SetTextSize(.08); leg.SetBorderSize(0); leg.SetNColumns(2);
            if tag=="CRW":
                h_WJets.SetFillStyle(0); h_WJets.SetLineWidth(2); h_WJets.SetLineColor(92); h_WJets.SetMarkerSize(0);
                c_WJets.SetFillStyle(0); c_WJets.SetLineWidth(2); c_WJets.SetLineColor(4); c_WJets.SetMarkerSize(0);
                c_data.SetFillStyle(0);  c_data.SetLineWidth(2); c_data.SetLineColor(1); c_data.SetMarkerSize(0); c_data.GetXaxis().SetNdivisions(509);
                c_WJets_sc= c_WJets.Clone("c_WJets_sc"); c_WJets_sc.Scale(h_WJets.Integral()/c_WJets.Integral());
                c_data_m_rest = c_data.Clone("c_data_m_rest"); c_data_m_rest.Add(c_Rest,-1); c_data_m_rest.Scale(h_WJets.Integral()/c_data_m_rest.Integral());  
                NormCRW = round( (c_data.Integral()-c_Rest.Integral())/c_WJets.Integral(), 2);  print "   --> NormCRW for",REGION,", (data-rest)/Wjets =", c_data.Integral(),"-",c_Rest.Integral(),"/",c_WJets.Integral(),"  =  " ,NormCRW," ,   prefit Norm. Sys. Unc. = ", 1-NormCRW;
                c_data_m_rest.Draw("pe"); c_WJets_sc.Draw("pe,same");  h_WJets.Draw("pe,same");  #  <========= DRAW Session ============
                maxY = TMath.Max( h_WJets.GetMaximum(), c_WJets_sc.GetMaximum() ); maxY = TMath.Max( maxY, c_data_m_rest.GetMaximum() ); h_WJets.GetYaxis().SetRangeUser(0, maxY*1.3 );  c_WJets.GetYaxis().SetRangeUser(0, maxY*1.3 ); c_data_m_rest.GetYaxis().SetRangeUser( 0, maxY*1.3 );
                Leg.AddEntry(h_WJets,"W+jets  S"+REGION[1:3],"F"); Leg.AddEntry(c_WJets,"W+jets  C"+REGION[1:4],"F");  Leg.AddEntry(c_data_m_rest, "Data-Rest(MC)  NS: "+str(NormCRW),"F");
                Ratio_SR_o_CR = h_WJets.Clone("Ratio_SR_o_CR");  Ratio_SR_o_CR.Divide(c_WJets_sc); Ratio_SR_o_CR.GetYaxis().SetRangeUser(0,2);  Ratio_SR_o_CR.GetYaxis().SetNdivisions( 504,0); Ratio_SR_o_CR.SetMarkerSize(0); Ratio_SR_o_CR.SetMarkerColor(92); Ratio_SR_o_CR.SetLineColor(92); #Ratio_SR_o_CR.SetErrorX(0.);
                Ratio_DmR_o_CR= c_data_m_rest.Clone( "Ratio_DmR_o_CR"); Ratio_DmR_o_CR.Divide(c_WJets_sc); leg.AddEntry(Ratio_DmR_o_CR, "CR(Data-Rest)/CR(W+jets)","F"); leg.AddEntry(Ratio_SR_o_CR, "SR(W+jets)/CR(W+jets)","F");                                                                #Ratio_DmR_o_CR.GetXaxis().SetError(0.);  
            if tag=="CRt":
                h_top.SetFillStyle(0); h_top.SetLineWidth(2); h_top.SetLineColor(8); h_top.SetMarkerSize(0);
                c_top.SetFillStyle(0); c_top.SetLineWidth(2); c_top.SetLineColor(4); c_top.SetMarkerSize(0);
                c_data.SetFillStyle(0); c_data.SetLineWidth(2); c_data.SetLineColor(1); c_data.SetMarkerSize(0); c_data.GetXaxis().SetNdivisions(509);
                c_top_sc = c_top.Clone("c_top_sc"); c_top_sc.Scale(h_top.Integral()/c_top.Integral());
                c_data_m_rest = c_data.Clone("c_data_m_rest"); c_data_m_rest.Add(c_Rest,-1); c_data_m_rest.Scale( h_top.Integral() / c_data_m_rest.Integral() );
                NormCRt = round( (c_data.Integral()-c_Rest.Integral())/(c_top.Integral()), 2);  print "   --> NormCRt for",REGION,", (data-rest)/(tt+sing.t) =", c_data.Integral(),"-",c_Rest.Integral(),"/",(c_top.Integral()),"  =  " ,NormCRt," ,   prefit Norm. Sys. Unc. = ", 1-NormCRt;
                c_data_m_rest.Draw("p,e"); c_top_sc.Draw("p,e,same"); h_top.Draw("p,e,same"); #  <========= DRAW Session ============
                maxY = TMath.Max( h_top.GetMaximum(), c_top_sc.GetMaximum() ); maxY = TMath.Max( maxY, c_data_m_rest.GetMaximum() ); h_top.GetYaxis().SetRangeUser(0, maxY*1.3 );  c_top.GetYaxis().SetRangeUser(0, maxY*1.3 ); c_data_m_rest.GetYaxis().SetRangeUser( 0, maxY*1.3 );
                Leg.AddEntry(h_top,"tt+sing.t  S"+REGION[1:3],"F"); Leg.AddEntry(c_top,"tt+sing.t  C"+REGION[1:4],"F");  Leg.AddEntry(c_data_m_rest, "Data-Rest(MC)  NS: "+str(NormCRt),"F");
                Ratio_SR_o_CR = h_top.Clone("Ratio_SR_o_CR");  Ratio_SR_o_CR.Divide(c_top_sc); Ratio_SR_o_CR.GetYaxis().SetRangeUser(0,2);  Ratio_SR_o_CR.GetYaxis().SetNdivisions( 504,0); Ratio_SR_o_CR.SetMarkerSize(0); Ratio_SR_o_CR.SetMarkerColor(8); Ratio_SR_o_CR.SetLineColor(8); #Ratio_SR_o_CR.GetXaxis().SetErrorX(0.);
                Ratio_DmR_o_CR= c_data_m_rest.Clone( "Ratio_DmR_o_CR"); Ratio_DmR_o_CR.Divide(c_top_sc); leg.AddEntry(Ratio_DmR_o_CR, "CR(Data-Rest)/CR(top)","F");  leg.AddEntry(Ratio_SR_o_CR, "SR(top)/CR(top)","F");                                                                    #Ratio_DmR_o_CR.GetXaxis().SetErrorX(0.);  
            fPads2.cd();  #self.tdrStyle.SetErrorX(0.);
            Ratio_DmR_o_CR.GetYaxis().SetRangeUser(0,2);  Ratio_DmR_o_CR.GetYaxis().SetNdivisions(504,0); Ratio_DmR_o_CR.SetFillStyle(0); Ratio_DmR_o_CR.SetLineWidth(2); Ratio_DmR_o_CR.SetLineColor(1); Ratio_DmR_o_CR.SetMarkerSize(0);
            Ratio_SR_o_CR.GetYaxis().SetTitle("Ratios   ");  Ratio_SR_o_CR.GetYaxis().SetTitleOffset(0.35); Ratio_SR_o_CR.GetYaxis().SetTitleSize(0.13);  Ratio_SR_o_CR.GetYaxis().SetTitleSize(0.13);
            Ratio_SR_o_CR.GetYaxis().SetLabelSize(0.11);     Ratio_SR_o_CR.GetXaxis().SetLabelSize(0.1);    Ratio_SR_o_CR.GetXaxis().SetTitleOffset(0.7); Ratio_SR_o_CR.GetXaxis().SetTitleSize(0.14);
            Ratio_SR_o_CR.Draw("e0"); Ratio_DmR_o_CR.Draw("e0,same"); Ratio_SR_o_CR.Draw("e0,same"); #  <========= DRAW Session ============
            axis1=TGaxis( min,1,max,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1); axis1.Draw("same");
            Leg.SetY1NDC(0.9-0.05*6-0.005);Leg.SetY1(Leg.GetY1NDC()); fPads1.cd(); Leg.Draw();fPads2.cd(); leg.Draw();fPads1.cd(); #leg.AddEntry(gr_MCStat, "Sys.","F");
            CMS = TLatex(0.11,0.96,"CMS    SR-vs-CR shapes comparison with normalized yields"); CMS.SetNDC(); CMS.SetTextSize(0.042); CMS.SetTextFont(42);CMS.Draw(); #CMS.SetTextAlign(31); CMS.SetLineWidth(2); 
                    

        if tag=="": #============ DRAW for SR PLOTS =====================
            h_TotalMC.Draw("e");h_TotalMC.GetXaxis().SetNdivisions(509);    #h_data.Draw("e same"); h_data.GetXaxis().SetNdivisions(509);
            hstack_TotalMC.Draw("same HIST");
            if "SR" not in REGION : h_data.Draw("same e");                  #h_data.Fit("pol3"); #BLINDING
            h_TotalMC.Draw("same e"   );
            h_Signal1.Draw("same HIST");
            h_Signal2.Draw("same HIST");
            h_Signal3.Draw("same HIST");
            h_Signal4.Draw("same HIST");
            h_Signal5.Draw("same HIST");
            h_Signal6.Draw("same HIST");  
   #h_data.GetXaxis().SetTitleOffset(1.2);h_data.GetYaxis().SetTitleOffset(1.3);h_data.GetXaxis().SetTitleSize(0.06);h_data.GetYaxis().SetTitleSize(0.06);h_data.GetXaxis().SetLabelSize(0.06);h_data.GetYaxis().SetLabelSize(0.06);
        #=========== ADD TEXT SESION ============================
            banner    = TLatex(0.96,0.96,str(lumi)+" fb^{-1} (13 TeV)");   banner.SetNDC();   banner.SetTextSize(0.034);   banner.SetTextFont(42);   banner.SetTextAlign(31);   banner.SetLineWidth(2);   banner.Draw();
            CMS       = TLatex(0.22,0.96,"CMS"                        );      CMS.SetNDC();      CMS.SetTextSize(0.042);      CMS.SetTextFont(42);      CMS.SetTextAlign(31);      CMS.SetLineWidth(2);      CMS.Draw();
            Extratext = TLatex(0.24,0.96,"Preliminary"                );Extratext.SetNDC();Extratext.SetTextSize(0.034);Extratext.SetTextFont(52);Extratext.SetTextAlign(11);Extratext.SetLineWidth(2);Extratext.Draw();
            RegionTxt = TLatex(0.20,0.88,"%s"%(REGION)                );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.044);RegionTxt.SetTextFont(42);RegionTxt.SetTextAlign(31);RegionTxt.SetLineWidth(2);RegionTxt.Draw();
            if channel=="mu"  and REGION[1:3] in ["R1","R2","R3",     "S1"] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 1,  %s-channel"%(channel));     Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
            if channel=="mu"  and REGION[1:3] in ["R4","R5","R6","R7","S2"] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 2,  %s-channel"%(channel));     Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
            if channel=="had" and REGION[1:3] in ["R1","R2","R3",     "S2"] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 2 or 3,  %s-channel"%(channel));     Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
            if channel=="had" and REGION[1:3] in ["R4","R5","R6","R7","S3"] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 3 or 4,  %s-channel"%(channel));     Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
            if channel=="mu"  and REGION      in ["PS0","SR"              ] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 1 or 2,  %s-channel"%(channel));Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
            if channel=="had" and REGION      in ["PS0","SR"              ] : Extratext1=TLatex(0.4, 0.96, "N_{j} = 2 or 3 or 4,  %s-channel"%(channel));Extratext1.SetNDC(); Extratext1.SetTextSize(0.034); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();

        #========== 2ND PAD SESSION ============================================
            if channel=="had": STcut=1250;
            if channel=="mu" : STcut=1000;
            if (REGION[0:2] in ["CR","PS","SR"]) or (options.ST_high<=STcut): #"" -->was "PS"
                fPads2.cd();            #fPads2.SetGridx(); #fPads2.SetGridy();
                h_Ratio = h_data.Clone("h_Ratio"); h_Ratio.Divide( h_TotalMC );  maxY=2;#TMath.Max( 2,  TMath.Min(3,h_Ratio.GetMaximum()*1.1) );
                h_Ratio.SetLineColor(1); h_Ratio.SetLineWidth(2); h_Ratio.SetMarkerStyle(8); h_Ratio.SetMarkerSize(0.7); h_Ratio.GetYaxis().SetRangeUser( 0 , maxY );  h_Ratio.GetYaxis().SetNdivisions(504,0);
                h_Ratio.GetYaxis().SetTitle("Data / MC  ");  h_Ratio.GetYaxis().SetTitleOffset(0.35);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetLabelSize(0.11); h_Ratio.GetXaxis().SetLabelSize(0.1); h_Ratio.GetXaxis().SetTitleOffset(0.7); h_Ratio.GetXaxis().SetTitleSize(0.14); 
                axis1=TGaxis( min,1,max,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1);  #axis1->SetLabelColor(16);
                for i in range(1,h_Ratio.GetNbinsX()+1,1):
                    D  = h_data.GetBinContent(i);    eD = h_data.GetBinError(i);
                    if D==0: eD=0.92;
                    B  = h_TotalMC.GetBinContent(i); eB = h_TotalMC.GetBinError(i);
                    if B<0.1 and eB>=B : eB=0.92;Err= 0.;
                    if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio.SetBinContent(i, D/B   );  h_Ratio.SetBinError(i, Err); #print i,")",h_Ratio.GetNbinsX()+1,")   data:",D," pm ",eD,"     Bkg:",B," pm ",eB,"   R:",D/B," pm ", Err
                    if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio.SetBinContent(i, D/0.92);  h_Ratio.SetBinError(i, Err);
                    if D==0 and B==0:                                                             h_Ratio.SetBinContent(i, -1);      h_Ratio.SetBinError(i, 0  );
                    if h_Ratio.GetBinContent(i)>maxY:h_Ratio.SetBinContent(i, maxY); ### To visualise the points above axis... #h_Ratio.Fit("pol1");
                    h_Ratio.Draw("e0"); axis1.Draw();
            #======= SIGNIFICANCES ON 2nd PAD =====================================
                sig1="";sig2="";sig3="";sig4="";sig5="";sig6="";
                if "SR" in REGION:
                    fPads2.SetLogy();  axis2=TGaxis( min,2,max,2, 0,0,0, "L"); axis2.SetLineColor(2); axis2.SetLineWidth(1);  axis3=TGaxis( min,5,max,5, 0,0,0, "L"); axis3.SetLineColor(3); axis3.SetLineWidth(1);
                    h_SqrtTotalMC=h_TotalMC.Clone("h_SqrtTotalMC"); h_SqrtTotalMC.Sumw2();
                    for i in range(1,h_SqrtTotalMC.GetNbinsX()+1,1): h_SqrtTotalMC.SetBinContent(i, TMath.Sqrt( 1.2*h_SqrtTotalMC.GetBinContent(i)+1 ) ); h_SqrtTotalMC.SetBinError( i, TMath.Sqrt( 1.2*h_SqrtTotalMC.GetBinError(i) ) );
                    h_Signif2=h_Signal2.Clone("h_Signif2");h_Signif2.Sumw2();h_Signif2.Divide( h_SqrtTotalMC ); h_Signif2.SetMarkerSize(0); 
                    if "SR" in REGION: h_Signif2.Draw("hist"); #axis2.Draw(); axis3.Draw(); 
                    h_Signif2.GetYaxis().SetTitle("S/\sqrt{1.2B+1}"); h_Signif2.GetYaxis().SetTitleOffset(0.35);h_Signif2.GetYaxis().SetTitleSize(0.13);h_Signif2.GetYaxis().SetTitleSize(0.13);h_Signif2.GetYaxis().SetLabelSize(0.11);h_Signif2.GetXaxis().SetLabelSize(0.1);h_Signif2.GetXaxis().SetTitleOffset(0.7);h_Signif2.GetXaxis().SetTitleSize(0.14);
                    h_Signif1=h_Signal1.Clone("h_Signif1"); h_Signif1.Sumw2(); h_Signif1.Divide( h_SqrtTotalMC ); h_Signif1.SetMarkerSize(0); 
                    h_Signif3=h_Signal3.Clone("h_Signif3"); h_Signif3.Sumw2(); h_Signif3.Divide( h_SqrtTotalMC ); h_Signif3.SetMarkerSize(0); 
                    h_Signif4=h_Signal4.Clone("h_Signif4"); h_Signif4.Sumw2(); h_Signif4.Divide( h_SqrtTotalMC ); h_Signif4.SetMarkerSize(0);
                    h_Signif5=h_Signal5.Clone("h_Signif5"); h_Signif5.Sumw2(); h_Signif5.Divide( h_SqrtTotalMC ); h_Signif5.SetMarkerSize(0);
                    h_Signif6=h_Signal6.Clone("h_Signif6"); h_Signif6.Sumw2(); h_Signif6.Divide( h_SqrtTotalMC ); h_Signif6.SetMarkerSize(0);
                    MaxY=1.6;MaxY=TMath.Max(MaxY,h_Signif6.GetMaximum());  MaxY=TMath.Max(MaxY,h_Signif5.GetMaximum());  MaxY=TMath.Max(MaxY,h_Signif1.GetMaximum()); MaxY=TMath.Max(MaxY,h_Signif2.GetMaximum());  MaxY=TMath.Max(MaxY,h_Signif3.GetMaximum());  MaxY=TMath.Max(MaxY,h_Signif4.GetMaximum());  #print MaxY;
                    h_Signif2.GetYaxis().SetRangeUser(0.01,MaxY*1.2); 
                    h_Signif2.Draw("hist");h_Signif1.Draw("hist,same");h_Signif2.Draw("hist,same");h_Signif3.Draw("hist,same");h_Signif4.Draw("hist,same");h_Signif5.Draw("hist,same");h_Signif6.Draw("hist,same"); axis2.Draw(); axis3.Draw(); fPads2.RedrawAxis(); fPads2.Update();
                    sig1=round(h_Signif1.Integral()*TMath.Sqrt(1.2)/options.scale1 ,1);
                    sig2=round(h_Signif2.Integral()*TMath.Sqrt(1.2)/options.scale1 ,1);
                    sig3=round(h_Signif3.Integral()*TMath.Sqrt(1.2)/options.scale1 ,1);
                    sig4=round(h_Signif4.Integral()*TMath.Sqrt(1.2)/options.scale2 ,1);
                    sig5=round(h_Signif5.Integral()*TMath.Sqrt(1.2)/options.scale2 ,1);
                    sig6=round(h_Signif6.Integral()*TMath.Sqrt(1.2)/options.scale2 ,1); print "\n --> Sum of Significances: s1, s2, s3, s4, s5, s6:    =   ",sig1,sig2,sig3,sig4,sig5,sig6

        #============= THE LEGEND SESSION =======================
            theLeg = TLegend(0.58, 0.45, 0.76, 0.92, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.04);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            if "SR" not in REGION:theLeg.AddEntry(h_data , "Data"  ,"ep");
            if channel=="had":theLeg.AddEntry(h_QCD    , "QCD"      ,"F");
            theLeg.AddEntry(h_WJets  , "W+jets "   ,"F");
            if channel=="had":theLeg.AddEntry(h_ZJets  , "Z+jets "   ,"F");
            theLeg.AddEntry(h_TTbar  , "t#bar{t}"  ,"F");
            theLeg.AddEntry(h_STop   , "single top","F");
            theLeg.AddEntry(h_VV     , "VV"        ,"F"); end=""; end_=""; ends1=""; ends2=""; ends3=""; ends4="";  ends1=""; ends2=""; ends3=""; ends4=""; 
            if options.scale1+options.scale2 is not 2: end=" #times %s"%(options.scale1);end_=" #times %s"%(options.scale2);
            if REGION in ["SR4","SR5","CR4q","CR5q","CR4t","CR5t"]:
               theLeg.AddEntry(h_Signal1,"(3000, 1200)"+end+"  "+str(sig1),"L");
               theLeg.AddEntry(h_Signal2,"(3000, 1500)"+end+"  "+str(sig2),"L");
               theLeg.AddEntry(h_Signal3,"(3000, 1800)"+end +"  "+str(sig3),"L");
               theLeg.AddEntry(h_Signal4,"(3500, 1400)"+end_ +"  "+str(sig4),"L");
               theLeg.AddEntry(h_Signal5,"(3500, 1750)"+end_ +"  "+str(sig5),"L");
               theLeg.AddEntry(h_Signal6,"(3500, 2100)"+end_ +"  "+str(sig6),"L");
            else:
               theLeg.AddEntry(h_Signal1,"(3000, 180)"+end+"  "+str(sig1),"L");
               theLeg.AddEntry(h_Signal2,"(3000, 240)"+end+"  "+str(sig2),"L");
               theLeg.AddEntry(h_Signal3,"(3000, 300)"+end +"  "+str(sig3),"L");
               theLeg.AddEntry(h_Signal4,"(3500, 210)"+end_ +"  "+str(sig4),"L");
               theLeg.AddEntry(h_Signal5,"(3500, 280)"+end_ +"  "+str(sig5),"L");
               theLeg.AddEntry(h_Signal6,"(3500, 350)"+end_ +"  "+str(sig6),"L");
            theLeg.SetY1NDC(0.9-0.05*6-0.005);theLeg.SetY1(theLeg.GetY1NDC()); fPads1.cd(); theLeg.Draw(); #theLeg.AddEntry(gr_MCStat, "Sys.","F");



        #============ SET MAX Y-AXIS FOR PLOTS ==================
        histsigmax = TMath.Max( h_Signal1.GetMaximum(), h_Signal2.GetMaximum() );
        histsigmax = TMath.Max( histsigmax, h_Signal3.GetMaximum() );
        histsigmax = TMath.Max( histsigmax, h_Signal4.GetMaximum() );
        histsigmax = TMath.Max( histsigmax, h_data.GetMaximum() );
        histsigmax = TMath.Max( histsigmax, h_TotalMC.GetMaximum() );
        h_Signal1.GetYaxis().SetRangeUser(0, histsigmax*1.2 );
        h_TotalMC.GetYaxis().SetRangeUser(0, histsigmax*1.2 );
        h_data.GetYaxis().SetRangeUser(   0, histsigmax*1.2 );


        #============ SAVE PLOTS IN A DIRECTORY ============================
        Directory=TString("%s_%s_%s"%(channel,REGION,REGION[3:4]));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/"); 
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p "+Directory.Data());        #only draw png
        variable1=options.name;
        file0 = TString(Directory.Data()+variable1+".pdf");
        file1 = TString(Directory.Data()+variable1+".png");
        if len(variable)>40:variable="long"; file1 = TString(Directory.Data()+variable+".png"); 
        fPads1.RedrawAxis(); fPads1.Update();
        canvas_controlplot.SaveAs( file0.Data() );
        canvas_controlplot.SaveAs( file1.Data() );        
        if logy: fPads1.Update(); canvas_controlplot.Update(); file0.ReplaceAll(".pdf","_log.pdf"); file1.ReplaceAll(".png","_log.png"); canvas_controlplot.SaveAs( file0.Data() ); canvas_controlplot.SaveAs( file1.Data() );#fPads2.Update(); #canvas_controlplot.SetLogy();  #fPads1.Setlogy();
        print "                   --> display %s_%s_%s/%s.png &"%(channel,REGION,REGION[3:4],variable); os.system(                   "display %s_%s_%s/%s.png &"%(channel,REGION,REGION[3:4],variable));

        #============= SAVE HISTOS TO ROOT FILE FOR LIMITS ====================        
        if channel=="had": STcut=1250;
        if channel=="mu" : STcut=1000;
        outf = TFile("LimitsInput_"+channel+".root","update")#,"recreate");
        if (REGION in ["SR1","SR2","SR3","SR6","SR7","CR1t","CR2t","CR3t","CR6t","CR7t","CR1W","CR2W","CR3W","CR6W","CR7W"] and (variable=="m_jlv" or variable=="MassVV[0]")) or (REGION in ["SR4","SR5","CR4t","CR5t","CR4W","CR5W"] and variable=="m_lvj"): 
            if tag=="":   # Store SR and CR raw spectra
                for h in [ h_WJets,h_TTbar,h_STop,h_VV,    h_TotalMC,     h_Signal1,h_Signal2,h_Signal3,h_Signal4,  ]:  New_h=h.Clone(  channel+"_"+REGION+"_"+variable+"__"+h.GetName()[2:len(h.GetName())]  ); New_h.Write();
                if "SR" in REGION and options.ST_high>STcut: h_data=h_TotalMC.Clone(channel+"_"+REGION+"_"+variable+"__DATA"); Integerization(h_data); h_data.Write();#  BLINDING  !!!
                h_top = h_TTbar.Clone(    channel+"_"+REGION+"_"+variable+"__top" ); h_top.Add(h_STop);  h_top.Write();
                if channel=="had": 
                    for h in [h_ZJets,h_QCD,]:  New_h=h.Clone(  channel+"_"+REGION+"_"+variable+"__"+h.GetName()[2:len(h.GetName())]  ); New_h.Write();
            if tag=="CRW":
                c_data_m_rest.SetBinContent( 0, abs(1-NormCRW) );
                h = c_data_m_rest.Clone( channel+"_"+REGION+"_"+variable+"__WJets__Shape_DmR" ); h.Write();
                h = c_WJets_sc.Clone(    channel+"_"+REGION+"_"+variable+"__WJets__ShapeMCCR" ); h.Write();
            if tag=="CRt":
                c_data_m_rest.SetBinContent( 0, abs(1-NormCRt) );
                h = c_data_m_rest.Clone( channel+"_"+REGION+"_"+variable+"__top__Shape_DmR" ); h.Write();
                h = c_top_sc.Clone(      channel+"_"+REGION+"_"+variable+"__top__ShapeMCCR" ); h.Write();
        outf.Close(); print "---> "+REGION+" observables in LimitsInput_"+channel+".root"



        #============== PIE CHARTS =================================
        if options.piechart:
            num_events, colors = array('d'), array('i');   gStyle.SetOptStat(000000000);            
            piecanvas=TCanvas("PIES","PIES",400,400);  piecanvas.SetTickx(1); piecanvas.SetTicky(1); piecanvas.SetRightMargin(-0.5); piecanvas.SetTopMargin(-0.5); piecanvas.SetLeftMargin(-0.5); piecanvas.SetBottomMargin(-0.5);gStyle.SetOptStat(000000000);
            if channel=="mu" :
                num_events.append(h_WJets.Integral());   num_events.append(h_TTbar.Integral()); num_events.append(h_STop.Integral());  num_events.append(h_VV.Integral());
                colors.append(self.color_palet["WJets"]);colors.append(self.color_palet["TTbar"]);colors.append(self.color_palet["STop"]);colors.append(self.color_palet["VV"]);
                pieplot=TPie("PIE",REGION,4,num_events,colors);
                pieplot.SetEntryLabel(0,"W+jets");pieplot.SetEntryLabel(1,"t#bar{t}");pieplot.SetEntryLabel(2,"single t");pieplot.SetEntryLabel(3,"VV");
            if channel=="had":
                num_events.append(h_QCD.Integral());   num_events.append(h_WJets.Integral());   num_events.append(h_ZJets.Integral()); num_events.append(h_TTbar.Integral()); num_events.append(h_STop.Integral()); num_events.append(h_VV.Integral());
                colors.append(self.color_palet["QCD"]);   colors.append(self.color_palet["WJets"]);   colors.append(self.color_palet["ZJets"]); colors.append(self.color_palet["TTbar"]); colors.append(self.color_palet["STop"]); colors.append(self.color_palet["VV"]);
                pieplot=TPie("PIE","",6,num_events,colors);
                pieplot.SetEntryLabel(0,"QCD");pieplot.SetEntryLabel(1,"W+jets");pieplot.SetEntryLabel(2,"Z+jets");pieplot.SetEntryLabel(3,"t#bar{t}");pieplot.SetEntryLabel(4,"single t");pieplot.SetEntryLabel(5,"VV");
            pieplot.SetTextSize(.045); pieplot.SetAngularOffset(30); pieplot.SetLabelFormat("%val (%perc) %txt"); pieplot.SetRadius(.4); pieplot.SetLabelsOffset(-.33); piecanvas.cd(1);pieplot.Draw("nol <"); file2=TString(Directory.Data()+"PIE.png"); piecanvas.SaveAs(file2.Data());
        

################# Main Code ################################
def Draw_Control_Plot( channel ) :
    if channel in ["had"]           : Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_0lep();
    if channel in ["mu","el","lep"] : Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_1lep();

if __name__ == '__main__' : 
    Beginning = strftime("%H:%M:%S",gmtime())
    print '____START__________channel:[',options.channel,']____ region:[',options.REGION,']______________(',Beginning,')________'
    Draw_Control_Plot( options.channel );
    Finishing = strftime("%H:%M:%S",gmtime());
    #========== CALCULATE DURATION OF THE RUN ===========
    MIN=int(Finishing[3:5])-int(Beginning[3:5]); SEC=int(Finishing[6:8])-int(Beginning[6:8]); 
    if SEC<0 and MIN>0 : SEC=60+SEC; MIN=MIN-1;
    if SEC>0 and MIN<0 : MIN=60+MIN;
    if SEC<0 and MIN<0 : SEC=60+SEC; MIN=60+MIN-1;
    print '_____END_____________________________________(',Finishing,', time:',MIN,'min',SEC,'sec)________\n'
