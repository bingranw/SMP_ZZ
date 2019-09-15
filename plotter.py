
import argparse
import yaml
import os
from array import *
from ROOT import *
import numpy
import random
import csv



# ==================================================
# Setting up
# ==================================================
parser = argparse.ArgumentParser()
parser.add_argument("-y", "--year", type=str, default='2016', help="which year")
options = parser.parse_args()

if options.year=='2016':
    intLum = 35.9
    sample_file = 'Config_MC_2016.yml'
    yml_file = 'ROOTfiles_2016.yml'
    xsec_file = 'xsections_2016.yaml'
elif options.year=='2017':
    intLum = 41.5
    sample_file = 'Config_MC_2017.yml'
    yml_file = 'ROOTfiles_2017.yml'
    xsec_file = 'xsections_2017.yaml'
else:
    print '2016 or 2017?'

with open(yml_file, 'r') as f_yml: 
    dict_files = yaml.load(f_yml)
with open(sample_file, 'r') as f_yml: 
    config_MC = yaml.load(f_yml)
with open(xsec_file, 'r') as f_yml: 
    dict_xsec = yaml.load(f_yml)

gROOT.SetBatch(1)
TH1.SetDefaultSumw2()



# ==================================================
# Get the x-axis label from the variable to draw
# ==================================================
def VarName(var_to_draw):
    name_dic = {
    'emulatedMET':'Emulated MET [GeV]',
    'Z_pt':'Dilepton pT [GeV]',
    }
    return name_dic[var_to_draw] if var_to_draw in name_dic else var_to_draw


# ==================================================
# Return histogram from a list of files
# ==================================================
def MakeHist(var_to_draw, selections, file_list, nbins_h, xmin_h, xmax_h):
    hist_name = str(random.randint(100, 999))
    hist = TH1F(hist_name, '', nbins_h, xmin_h, xmax_h)

    chain = TChain("Events")
    for fn in file_list:
        chain.Add(fn)
    chain.SetProof()
    chain.Project(hist_name, var_to_draw, selections)
    chain.Reset()
    chain.Delete()

    return hist


# ==================================================
# Return a stack of MC histograms
# ==================================================
def MakeStack(var_to_draw, cuts, nbins_h, xmin_h, xmax_h):
    selections = "(%s) * puWeight * weight" % (' && '.join(cuts))
    hists = []
    hist_integrals = []
    stk_title = ';;Events / bin'
    stk = THStack("stk", ";%s;Events / bin" % VarName(var_to_draw))

    for sample_name in config_MC.keys():
        if sample_name=='ZZ_EWK':
            continue
        print '\n=================================================='
        print 'Doing sample:', sample_name
        print '=================================================='

        hist = TH1F(sample_name, stk_title, nbins_h, xmin_h, xmax_h)
        hist_color = TColor.GetColor( config_MC[sample_name]['color'] )
        hist.SetFillColor( hist_color )
        hist.SetLineColor( kBlack )
        for dataset in config_MC[sample_name]['dataset']:
            file_list = dict_files[dataset]['files']
            _hist = MakeHist(var_to_draw, selections, file_list, nbins_h, xmin_h, xmax_h)
            _xsec = dict_xsec[dataset]['xsec']
            _hist.Scale(intLum * dict_files[dataset]['GWCorr'] * 1000. * _xsec)
            hist.Add( _hist )
        hists.append( hist )
        hist_integrals.append( hist.Integral() )

    for i in numpy.argsort(hist_integrals):
        stk.Add(hists[i])

    return stk


# ==================================================
# Return the histogram of signal
# ==================================================
def MakeSignal(var_to_draw, cuts, nbins_h, xmin_h, xmax_h):
    print '\n=================================================='
    print 'Doing Signal: ZZJJ'
    print '=================================================='
    selections = "(%s) * puWeight * weight" % (' && '.join(cuts))
    sample_name='ZZ_EWK'

    hist = TH1F(sample_name, '', nbins_h, xmin_h, xmax_h)
    hist_color = TColor.GetColor( config_MC[sample_name]['color'] )
    hist.SetFillStyle( 0 )
    hist.SetLineColor( hist_color )
    for dataset in config_MC[sample_name]['dataset']:
        file_list = dict_files[dataset]['files']
        _hist = MakeHist(var_to_draw, selections, file_list, nbins_h, xmin_h, xmax_h)
        _xsec = dict_xsec[dataset]['xsec']
        _hist.Scale(intLum * dict_files[dataset]['GWCorr'] * 1000. * _xsec)
        hist.Add( _hist )
    return hist


# ==================================================
# Return the histogram of data
# ==================================================
def MakeData(var_to_draw, cuts, nbins_h, xmin_h, xmax_h, draw_data=False):
    print '\n=================================================='
    print 'Doing Data'
    print '=================================================='
    selections_data = "(%s)" % ' && '.join(cuts)
    hist = TH1F('hist_data', '', nbins_h, xmin_h, xmax_h)
    hist.SetMarkerStyle(kFullDotLarge)
    hist.SetMarkerSize(0.5)
    hist.SetLineColor(kBlack)
    if not draw_data:
        return hist

    datasets = ['SingleElectron', 'SingleMuon', 'DoubleEG', 'DoubleMuon', 'MuonEG']
    chain = TChain("Events")
    for dataset in dict_files.keys(): 
        if not [_dataset for _dataset in datasets if _dataset in dataset]:
            continue
        for fn in dict_files[dataset]['files']: 
            chain.Add(fn)
    chain.SetProof()
    chain.Project('hist_data', var_to_draw, selections_data)
    chain.Reset()
    chain.Delete()

    return hist


# ==================================================
# Make MC and Data plot
# ==================================================
def MakePlot(stk, hist_signal, hist_data, region, save_name, draw_data=False):
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

    # MC plot 
    pad_stack.cd()
    pad_stack.SetBottomMargin(0.)
    pad_stack.SetTopMargin(0.07)
    pad_stack.SetRightMargin(0.03)
    stk.Draw('hist')
    stk.SetMinimum( 0.1 )
    hist_signal.Draw("same hist")
    lg = pad_stack.BuildLegend(0.15, 0.80, 0.93, 0.92, "", "fNDC")
    lg.SetNColumns(4)
    lg.SetFillStyle(0)
    lg.SetBorderSize(0)
    pad_stack.SetLogy()

    # Data plot
    if draw_data:
        hist_data.Draw("EP Same")
        std_max = 10 ** ( numpy.log10( max(hist_data.GetMaximum(),stk.GetMaximum()) ) * 1.2 )
    else:
        std_max = 10 ** ( numpy.log10( stk.GetMaximum() ) * 1.2 )
    stk.SetMaximum( std_max )

    # Put info on top
    pad_stack.Modified()
    pad_stack.Update()
    xMin = pad_stack.GetFrame().GetX1()
    xMax = pad_stack.GetFrame().GetX2()
    yMin = pad_stack.GetFrame().GetY1()
    yMax = pad_stack.GetFrame().GetY2()
    axisTop = TGaxis(xMin, pow(10,yMax), xMax, pow(10,yMax), 0, 1, 1, "-")
    axisTop.ChangeLabel(1,-1,0.035,11,-1,-1, "Preliminary")
    axisTop.ChangeLabel(-1,-1,0.035,31,-1,-1, "%s, %s fb^{-1} (13 TeV), %s" % (region, intLum, options.year))
    axisTop.Draw()

    # Ratio plot
    pad_ratio.cd()
    pad_ratio.SetTopMargin(0.)
    pad_ratio.SetBottomMargin(0.25)
    pad_ratio.SetRightMargin(0.03)

    sumMCErrors = stk.GetHistogram()
    dataOverSumMC = hist_data.Clone('data_ov_MC')
    if draw_data:
        dataOverSumMC.Divide(sumMCErrors)
    for i in range(hist_data.GetNbinsX()+2): 
        dataOverSumMC.SetBinError(i, hist_data.GetBinError(i)/max(hist_data.GetBinContent(i), 1))

    stk.GetXaxis().Copy(dataOverSumMC.GetXaxis())
    dataOverSumMC.GetYaxis().SetTitle('Data / MC')
    dataOverSumMC.GetYaxis().CenterTitle()
    dataOverSumMC.GetYaxis().SetRangeUser(0.3, 1.7)
    dataOverSumMC.GetYaxis().SetNdivisions(503)
    dataOverSumMC.GetYaxis().SetTitleOffset(0.5)
    dataOverSumMC.GetYaxis().SetTitleSize(dataOverSumMC.GetYaxis().GetTitleSize()*2.5)
    dataOverSumMC.GetYaxis().SetLabelSize(dataOverSumMC.GetYaxis().GetLabelSize()*2.5)
    dataOverSumMC.GetXaxis().SetTitleSize(dataOverSumMC.GetXaxis().GetTitleSize()*2.5)
    dataOverSumMC.GetXaxis().SetLabelSize(dataOverSumMC.GetXaxis().GetLabelSize()*2.5)
    dataOverSumMC.Draw()

    xaxis = dataOverSumMC.GetXaxis()
    line = TLine(xaxis.GetBinLowEdge(xaxis.GetFirst()), 1, xaxis.GetBinUpEdge(xaxis.GetLast()), 1)
    line.SetLineStyle(kDotted)
    line.Draw()

    cOutput.SaveAs(save_name)


# ==================================================
# Decide what to draw
# ==================================================
def CR_tmp(draw_data=False):
    cuts = [
    '(lep_category==1 || lep_category==3)', 
    'abs(Z_mass-91.1876)<15',
    'abs(delta_R_ll)<2',
    'sca_balance>0.6',
    'sca_balance<2.5',
    'Z_pt>60',
    'ngood_jets>=2',
    'ngood_bjets==0',
    'met_pt>90',
    ]
    nbins_h = 15
    xmin_h = 100
    xmax_h = 1600
    var_to_draw = 'dijet_Mjj'
    region = "ll channel"
    save_name = "plots/test_%s.pdf" % options.year

    stk = MakeStack(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_signal = MakeSignal(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_data = MakeData(var_to_draw, cuts, nbins_h, xmin_h, xmax_h, draw_data)
    MakePlot(stk, hist_signal, hist_data, region, save_name, draw_data)

def Draw_nm1(draw_data=False):
    dict_nm1 = {
        "Z_mass":               [60,    120,    60],
        "Z_pt":                 [0,     1000,   50],
        "met_pt":               [0,     500,    50],
        "sca_balance":          [0,     10,     50],
        "delta_R_ll":           [0.5,   2.5,    40],
        "ngood_jets":           [0,     10,     10],
        "abs(delta_phi_ZMet)":  [0,     4,      40],
        "abs(delta_phi_j_met)": [0,     4,      40],
        "H_T":                  [0,     2000,   40],
        "CJV_Pt_Sum":           [0,     400,    40],
        "x_jet30":              [0,     1,      20],
        "abs(dijet_Zep)":       [0,     5,      50],
        "dijet_centrality_gg":  [0,     1,      20],
        "Jet_pt_Ratio":         [0,     1,      20],
        "lead_jet_pt":          [0,     500,    50],
        "trail_jet_pt":         [0,     500,    50],
        "Jet_etas_multiplied":  [0,     20,     40]
    }
    cuts = [
    '(lep_category==1 || lep_category==3)', 
    'nhad_taus==0',
    'ngood_jets>=2',
    'ngood_bjets==0',
    'abs(Z_mass-91.1876)<15',
    'Z_pt>50',
    'met_pt>50',
    'sca_balance>0. && sca_balance<6.',
    'delta_R_ll<2',
    'abs(delta_phi_ZMet)>1',
    'abs(delta_phi_j_met)>1',
    # 'H_T<1500',
    # 'CJV_Pt_Sum>0',
    # 'x_jet30>0.1 && x_jet30<0.9',
    # 'abs(dijet_Zep)<4',
    # 'dijet_centrality_gg>0',
    # 'Jet_pt_Ratio>0',
    'lead_jet_pt>70',
    'trail_jet_pt>50',
    'Jet_etas_multiplied>0'
    ]
    var_to_draw = 'delta_R_ll'
    cuts = [cut for cut in cuts if var_to_draw not in cut]
    xmin_h, xmax_h, nbins_h = dict_nm1[var_to_draw]
    region = "ll channel"
    save_name = "plots/VBS_%s.pdf" % var_to_draw

    stk = MakeStack(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_signal = MakeSignal(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_data = MakeData(var_to_draw, cuts, nbins_h, xmin_h, xmax_h, draw_data)
    MakePlot(stk, hist_signal, hist_data, region, save_name, draw_data)

def Draw_VBS(draw_data=True):
    dict_var = {
        'dijet_Mjj' : [0, 1500, 15, ['dijet_Mjj<400']],
        'dijet_abs_delta_eta' : [0, 10, 10, ['dijet_abs_delta_eta<2.4']]
    }
    cuts = [
    '(lep_category==1 || lep_category==3)', 
    'nhad_taus==0',
    'ngood_jets>=2',
    'ngood_bjets==0',
    'abs(Z_mass-91.1876)<15',
    'Z_pt>50',
    'met_pt>100',
    'sca_balance>0. && sca_balance<6.',
    'delta_R_ll<2',
    'abs(delta_phi_ZMet)>1',
    'abs(delta_phi_j_met)>1',
    'H_T<1500',
    'CJV_Pt_Sum>0',
    'x_jet30>0.1 && x_jet30<0.9',
    'abs(dijet_Zep)<4',
    'dijet_centrality_gg>0',
    'Jet_pt_Ratio>0',
    'lead_jet_pt>70',
    'trail_jet_pt>50',
    'Jet_etas_multiplied>0'
    ]
    var_to_draw = 'delta_R_ll'
    xmin_h, xmax_h, nbins_h, cut_add = dict_var[var_to_draw]
    region = "ll channel"
    save_name = "plots/VBS_%s.pdf" % var_to_draw

    stk = MakeStack(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_signal = MakeSignal(var_to_draw, cuts, nbins_h, xmin_h, xmax_h)
    hist_data = MakeData(var_to_draw, cuts+cut_add, nbins_h, xmin_h, xmax_h, draw_data)
    MakePlot(stk, hist_signal, hist_data, region, save_name, draw_data)


# ==================================================
# Main program
# ==================================================
def main():
    prf = TProof.Open("lite://")

    # Draw_VBS(True)
    Draw_nm1(False)

    prf.Close()
    prf.Delete()
    gROOT.Reset()


if __name__ == "__main__":
    main()



