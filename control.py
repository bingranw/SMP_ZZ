import argparse
import yaml
import os
from array import *
from ROOT import *
import numpy
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-y", "--year", type=str, default='2017', help="which year")
options = parser.parse_args()

if options.year=='2016':
    intLum = 35.9
    sample_file = 'MC_samples_2016.yml'
    yml_file = '../../python/postprocessing/monoZ/tools/ROOTfiles_2016.yml'
elif options.year=='2017':
    intLum = 41.5
    sample_file = 'MC_samples_2017.yml'
    yml_file = '../../python/postprocessing/monoZ/tools/ROOTfiles_2017.yml'
else:
    print '2016 or 2017?'


gROOT.SetBatch(1)
TH1.SetDefaultSumw2()

with open(yml_file, 'r') as f_yml: 
    dict_yml = yaml.load(f_yml)

with open(sample_file, 'r') as f_yml: 
    MC_samples = yaml.load(f_yml)


def VarName(var_to_draw):
    name_dic = {
    'emulatedMET':'Emulated MET [GeV]',
    'emulatedMET':'Emulated MET [GeV]',
    'Z_pt':'Dilepton pT [GeV]',
    }
    return name_dic[var_to_draw] if var_to_draw in name_dic else var_to_draw


def DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name):
    # ==================================================
    # Setup tproof
    # ==================================================

    prf = TProof.Open("lite://")

    # ==================================================
    # MC
    # ==================================================

    hists = []
    hist_integrals = []
    stk_title = ';;Events / bin'
    stk = THStack("stk", ";%s;Events / bin" % VarName(var_to_draw))

    sumMCErrors = TH1F('sumMCErrors_hist', stk_title, nbins_h, xmin_h, xmax_h)
    sumMCErrors.SetFillColorAlpha(kGray, 0.5)
    sumMCErrors.SetMarkerSize(0)

    # pT binned DY
    hist_DY = TH1F('DYToLL', stk_title, nbins_h, xmin_h, xmax_h) 
    hist_DY.SetFillColor( TColor.GetColor('#bebada') )
    hist_DY.SetLineColor( kBlack )
    hist_DY_integral = 0

    # Single top
    hist_ST = TH1F('ST_tW', stk_title, nbins_h, xmin_h, xmax_h) 
    hist_ST.SetFillColor( TColor.GetColor('#b3de69') )
    hist_ST.SetLineColor( kBlack )
    hist_ST_integral = 0

    # all others
    hist_Other = TH1F('Others', stk_title, nbins_h, xmin_h, xmax_h) 
    hist_Other.SetFillColor( TColor.GetColor('#d9d9d9') )
    hist_Other.SetLineColor( kBlack )
    hist_Other_integral = 0

    hist_names = []
    yield_dic = {}
    for dataset in MC_samples.keys():
        # if dataset[:2] != 'DY':
        #     continue
        hist_name = dataset.split('_')[0]
        if dataset[:3] == 'Glu':
            hist_name = 'GG' + hist_name.split('To')[-1]
        while hist_name in hist_names:
            hist_name += '_'
        hist_names.append(hist_name)

        selections = "(%s) * puWeight * weight" % (' && '.join(cuts))
        selections_data = "(%s)" % ' && '.join(cuts)

        hist = TH1F(hist_name, stk_title, nbins_h, xmin_h, xmax_h) 
        chain = TChain("Events")
        for file in dict_yml[dataset]['files']: 
            chain.Add(file)
        chain.SetProof()
        print 'Doing dataset', dataset
        chain.Project(hist_name, var_to_draw, selections)

        # Scaling MC by a number, not sure why 
        if dataset[:2] == 'WZ' and options.year=='2017':
            hist.Scale( intLum * dict_yml[dataset]['GWCorr'] * 0.85 )
        elif dataset[:2] == 'DY' and options.year=='2017':
            hist.Scale( intLum * dict_yml[dataset]['GWCorr'] / 1.35 )
        elif dataset[:3] == 'Glu':
            hist.Scale( intLum * dict_yml[dataset]['GWCorr'] / 1000. )
        else:
            hist.Scale( intLum * dict_yml[dataset]['GWCorr'] )

        # Hist style
        hist_color = TColor.GetColor( MC_samples[dataset]['color'] )
        hist.SetFillColor( hist_color )
        hist.SetLineColor( kBlack )

        # pT binned or not
        if dataset[:2] == 'DY':
            hist_DY.Add( hist )
            hist_DY_integral += hist.Integral() 
        elif dataset[:3] == 'ST_':
            hist_ST.Add( hist )
            hist_ST_integral += hist.Integral() 
        elif dataset[:3] in ['Glu', 'ZZZ', 'WZZ', 'WWZ', 'TTZ', 'TTW', 'ZGT'] or hist_name=='ZZTo2L2Q':
            hist_Other.Add( hist )
            hist_Other_integral += hist.Integral() 
        else:
            hists.append( hist )
            hist_integrals.append( hist.Integral() )

        yields = [hist.GetBinContent(_bin) for _bin in range(hist.GetNbinsX()+2) ]
        yield_dic[dataset] = yields
        chain.Reset()
        chain.Delete()

    writer = csv.writer(open("bingran.csv", "w"))
    for key, val in yield_dic.items():          
        writer.writerow([key]+val)

    hists.append( hist_DY )
    hist_integrals.append( hist_DY_integral )
    hists.append( hist_ST )
    hist_integrals.append( hist_ST_integral )
    hists.append( hist_Other )
    hist_integrals.append( hist_Other_integral )

    for i in numpy.argsort(hist_integrals):
        stk.Add(hists[i])
        sumMCErrors.Add(hists[i])

    # ==================================================
    # Data
    # ==================================================

    datasets = ['SingleElectron', 'SingleMuon', 'DoubleEG', 'DoubleMuon', 'MuonEG']
    hist_data = TH1F('hist_data', stk_title, nbins_h, xmin_h, xmax_h) 
    chain = TChain("Events")
    for dataset in datasets: 
        for file in dict_yml[dataset]['files']: 
            chain.Add(file)
    chain.SetProof()
    chain.Project('hist_data', var_to_draw, selections_data)
    # hist_data.Scale(1.)
    chain.Reset()
    chain.Delete()

    # ==================================================
    # Stack plot
    # ==================================================

    cOutput = TCanvas("cOutput", "cOutput", 500, 500)
    cOutput.SetFillStyle(4000)
    cOutput.SetFrameFillStyle(4000)
    gStyle.SetOptStat(0)

    pad_stack = TPad("pad_stack","", 0., 0.25, 1., 1.)
    pad_ratio = TPad("pad_ratio","", 0., 0., 1., 0.25)
    pad_stack.SetFillStyle(4000)
    pad_stack.SetFrameFillStyle(4000)
    pad_ratio.SetFillStyle(4000)
    pad_ratio.SetFrameFillStyle(4000)
    pad_stack.Draw()
    pad_ratio.Draw()

    pad_stack.cd()
    pad_stack.SetBottomMargin(0.)
    pad_stack.SetTopMargin(0.07)
    pad_stack.SetRightMargin(0.03)
    stk.Draw("hist")
    stk.SetMinimum( 1. )
    stk.GetXaxis().SetTitleOffset(1.1)
    stk.GetYaxis().SetTitleOffset(1.3)

    lg = pad_stack.BuildLegend(0.15, 0.80, 0.93, 0.92, "", "fNDC")
    lg.SetNColumns(4)
    lg.SetFillStyle(0)
    lg.SetBorderSize(0)
    # pad_stack.SetLogx()
    # stk.GetXaxis().SetMoreLogLabels()
    # stk.GetXaxis().SetNoExponent()
    pad_stack.SetLogy()

    # ==================================================
    # Data plot
    # ==================================================

    hist_data.SetMarkerStyle(kFullDotLarge)
    hist_data.SetMarkerSize(0.5)
    hist_data.SetLineColor(kBlack)
    hist_data.Draw("EP Same")
    std_max = 10 ** ( numpy.log10( max(hist_data.GetMaximum(),stk.GetMaximum()) ) * 1.2 )
    stk.SetMaximum( std_max )

    pad_stack.Modified()
    pad_stack.Update()
    xMin = pad_stack.GetFrame().GetX1()
    xMax = pad_stack.GetFrame().GetX2()
    yMin = pad_stack.GetFrame().GetY1()
    yMax = pad_stack.GetFrame().GetY2()
    # axisTop = TGaxis(xMin, yMax, xMax, yMax, 0, 1, 1, "-")
    axisTop = TGaxis(xMin, pow(10,yMax), xMax, pow(10,yMax), 0, 1, 1, "-")
    # axisTop = TGaxis(pow(10,xMin), pow(10,yMax), pow(10,xMax), pow(10,yMax), 0, 1, 1, "-")
    axisTop.ChangeLabel(1,-1,0.035,11,-1,-1, "Preliminary")
    axisTop.ChangeLabel(-1,-1,0.035,31,-1,-1, "%s, %s fb^{-1} (13 TeV), %s" % (region, intLum, options.year))
    axisTop.Draw()

    # oldBottomMargin = pad_stack.GetBottomMargin()
    # pad_stack.SetTopMargin(pad_stack.GetTopMargin()/0.7)


    # ==================================================
    # Ratio plot
    # ==================================================

    pad_ratio.cd()
    pad_ratio.SetTopMargin(0.)
    # pad_ratio.SetTopMargin(0.1)
    pad_ratio.SetBottomMargin(0.25)
    pad_ratio.SetRightMargin(0.03)

    dataOverSumMC = hist_data.Clone('data_ov_MC')
    dataOverSumMC.Divide(sumMCErrors)

    for i in range(hist_data.GetNbinsX()+2): 
        dataOverSumMC.SetBinError(i, hist_data.GetBinError(i)/max(hist_data.GetBinContent(i), 1))
        # sumMCErrors.SetBinError( i, numpy.sqrt(abs(sumMCErrors.GetBinContent(i))) / max(sumMCErrors.GetBinContent(i),1) )
        # sumMCErrors.SetBinError( i, abs(sumMCErrors.GetBinError(i)) / max(sumMCErrors.GetBinContent(i),1) )
        # sumMCErrors.SetBinContent(i, 1.)

    sumMCErrors.Divide(sumMCErrors)

    stk.GetXaxis().Copy(sumMCErrors.GetXaxis())
    sumMCErrors.GetYaxis().SetTitle('Data / MC')
    sumMCErrors.GetYaxis().CenterTitle()
    sumMCErrors.GetYaxis().SetRangeUser(0.3, 1.7)
    sumMCErrors.GetYaxis().SetNdivisions(503)
    sumMCErrors.GetYaxis().SetTitleOffset(0.5)
    sumMCErrors.GetYaxis().SetTitleSize(sumMCErrors.GetYaxis().GetTitleSize()*2.5)
    sumMCErrors.GetYaxis().SetLabelSize(sumMCErrors.GetYaxis().GetLabelSize()*2.5)
    sumMCErrors.GetXaxis().SetTitleSize(sumMCErrors.GetXaxis().GetTitleSize()*2.5)
    sumMCErrors.GetXaxis().SetLabelSize(sumMCErrors.GetXaxis().GetLabelSize()*2.5)

    sumMCErrors.Draw("E2")
    dataOverSumMC.Draw("same")

    xaxis = dataOverSumMC.GetXaxis()
    line = TLine(xaxis.GetBinLowEdge(xaxis.GetFirst()), 1, xaxis.GetBinUpEdge(xaxis.GetLast()), 1)
    line.SetLineStyle(kDotted)
    line.Draw()

    # ==================================================
    # Save plot and clean up
    # ==================================================
    cOutput.SaveAs(save_name)
    line.Delete()
    hist.Delete()
    hist_DY.Delete()
    dataOverSumMC.Delete()
    sumMCErrors.Delete()
    hist_data.Delete()
    pad_stack.Delete()
    pad_ratio.Delete()
    cOutput.Delete()
    prf.Close()
    prf.Delete()
    stk.Delete()


def CR_DY():
    cuts = [
    '(lep_category==1 || lep_category==3)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    'nhad_taus==0',
    'met_pt<40',
    'Z_pt>60',
    'abs(Z_mass-91.1876)<15',
    ]
    # nbins_h = 66
    # xmin_h = -3.3
    # xmax_h = 3.3
    # var_to_draw = 'delta_phi_ZMet'
    nbins_h = 67
    xmin_h = 30
    xmax_h = 700
    var_to_draw = 'Z_pt'
    # nbins_h = 50
    # xmin_h = 0
    # xmax_h = 5.0
    # var_to_draw = 'delta_R_ll'
    region = "ll channel"
    save_name = "plots/CR_%s_DY.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_EMU():
    cuts = [
    # '(lep_category==2)', 
    '(lep_category==2)', 
    'abs(Z_mass-91.1876)<15',
    'Z_pt>60',
    'ngood_jets<2',
    'ngood_bjets==0',
    'met_pt>40',
    # 'nhad_taus==0',
    # 'abs(Z_mass-91.1876)>15',
    ]
    nbins_h = 26
    xmin_h = 40
    xmax_h = 300
    # nbins_h = 20
    # xmin_h = 0
    # xmax_h = 500
    var_to_draw = 'Z_pt'
    region = "e mu channel"
    save_name = "plots/CR_%s_EMU.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_TT():
    cuts = [
    '(lep_category==2)', 
    'ngood_bjets==1',
    # 'ngood_bjets>1',
    'ngood_jets>2',
    # 'ngood_jets>1',
    'nhad_taus==0',
    'met_pt>50',
    'Z_pt>60',
    'abs(Z_mass-91.1876)<15',
    # 'Z_mass>50'
    ]
    # nbins_h = 66
    # xmin_h = -3.3
    # xmax_h = 3.3
    # var_to_draw = 'delta_phi_ZMet'
    # nbins_h = 15
    # xmin_h = 0
    # xmax_h = 300
    # var_to_draw = 'Z_pt'
    nbins_h = 31
    xmin_h = 40
    xmax_h = 350
    var_to_draw = 'Z_pt'

    region = "e mu channel"
    save_name = "plots/CR_%s_TT.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_3L():
    cuts = [
    '(lep_category==4 || lep_category==5)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    'nhad_taus==0',
    'Z_pt>60',
    'met_pt>50',
    'abs(Z_mass-91.1876)<15',
    ]
    nbins_h = 20
    xmin_h = 0
    xmax_h = 400
    var_to_draw = 'emulatedMET'
    # nbins_h = 50
    # xmin_h = 0
    # xmax_h = 5.0
    # var_to_draw = 'delta_R_ll'
    region = "3l channel"
    save_name = "plots/CR_%s_3L.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_4L():
    cuts = [
    '(lep_category==6 || lep_category==7)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    'nhad_taus==0',
    'Z_pt>60',
    # 'met_pt>10',
    'emulatedMET/Z_pt>0.4',
    'emulatedMET/Z_pt<2.5',
    'abs(Z_mass-91.1876)<15',
    ]
    nbins_h = 15
    xmin_h = 0
    xmax_h = 300
    var_to_draw = 'emulatedMET'
    region = "4l channel"
    save_name = "plots/CR_%s_4L.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_Top():
    cuts = [
    '(lep_category==2)', 
    '(ngood_bjets==1 || ngood_bjets==2)',
    # "ngood_bjets==2",
    'ngood_jets>=2',
    'nhad_taus==0',
    'Z_pt>60',
    'met_pt>50',
    # 'abs(Z_mass-91.1876)>15',
    # '(Z_mass>50 && Z_mass<200)',
    # 'abs(Z_mass-91.1876)<15',
    'Z_mass>50'
    ]
    # nbins_h = 50
    # xmin_h = 0
    # xmax_h = 5.0
    # var_to_draw = 'delta_R_ll'
    nbins_h = 46
    xmin_h = 40
    xmax_h = 500
    var_to_draw = 'Z_pt'
    # nbins_h = 72
    # xmin_h = 40
    # xmax_h = 400
    # var_to_draw = 'Z_mass'
    region = "e mu channel"
    save_name = "plots/CR_%s_Top.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_WW():
    cuts = [
    '(lep_category==2)', 
    # 'ngood_bjets==0',
    # 'ngood_jets<2',
    'ngood_jets==0',
    'nhad_taus==0',
    'Z_pt>60',
    'met_pt>50',
    # 'abs(delta_phi_j_met)>0.5',
    # 'sca_balance>0.4',
    # 'sca_balance<2.5',
    'Z_mass>50'
    ]
    # nbins_h = 72
    # xmin_h = 40
    # xmax_h = 400
    # var_to_draw = 'Z_mass'
    nbins_h = 16
    xmin_h = 30
    xmax_h = 350
    var_to_draw = 'Z_pt'
    region = "e mu channel"
    save_name = "plots/CR_%s_WW.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_tmp():
    cuts = [
    '(lep_category==5)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    'nhad_taus==0',
    'Z_pt>60',
    'met_pt>50',
    'abs(Z_mass-91.1876)<15',
    ]
    nbins_h = 20
    xmin_h = 0
    xmax_h = 400
    var_to_draw = 'emulatedMET'
    # nbins_h = 50
    # xmin_h = 0
    # xmax_h = 5.0
    # var_to_draw = 'delta_R_ll'
    region = "mml channel"
    save_name = "plots/CR_%s_3L_mml.pdf" % options.year
    DrawVar(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def main():
    # CR_3L()
    # CR_4L()
    # CR_DY()
    # CR_EMU()
    # CR_TT()
    # CR_WW()
    CR_tmp()
    gROOT.Reset()
    # CR_WZ()
    # CR_Top()


main()



