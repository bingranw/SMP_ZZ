
import yaml
import os
from array import *
from ROOT import *


intLum = 41.5
gROOT.SetBatch(1)

with open('ROOTfiles.yml', 'r') as f_yml: 
    dict_yml = yaml.load(f_yml)

with open('ROOTfiles_data.yml', 'r') as f_yml: 
    dict_data_yml = yaml.load(f_yml)

with open('genWeightCorrection.yml', 'r') as f_yml: 
    genWeightCorrection = yaml.load(f_yml)

def Draw_pT(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name):
    # ==================================================
    # Setup chain and cuts
    # ==================================================
    MC_samples = [
    'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',
    'ZZTo2L2Nu_13TeV_powheg_pythia8',
    'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_low',
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
    'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'ZZTo4L_13TeV_powheg_pythia8',
    # 'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',
    # 'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8',
    # 'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',
    # 'GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8',
    # 'GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8'
    ]

    hex_cols = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']
    ndatasets = len(MC_samples)


    prf = TProof.Open("lite://")

    # ==================================================
    # MC
    # ==================================================
    stk_title = ';;Events / bin'
    stk = THStack("stk", ";%s;Events / bin" % var_to_draw)
    for i in range(ndatasets):
        dataset = MC_samples[i]
        hist_name = dataset.split('_')[0]
        selections = "puWeight * lumiWeight * genWeight * %f * (%s)" % (genWeightCorrection[dataset], ' && '.join(cuts))
        # selections = "puWeight * lumiWeight * (%s)" % (' && '.join(cuts))
        selections_data = "(%s)" % ' && '.join(cuts)

        hist = TH1D(hist_name, stk_title, nbins_h, xmin_h, xmax_h) 
        chain = TChain("Events")
        for file in dict_yml[dataset]['files']: 
            chain.Add(file)
        chain.SetProof()
        chain.Project(hist_name, var_to_draw, selections)
        hist.Scale(intLum)
        # if dataset=='TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8':
        #    hist.Scale(687.1/88.29)
        hist.SetFillColor( TColor.GetColor(hex_cols[i]) )
        hist.SetLineColor(kBlack)
        stk.Add(hist)
        chain.Reset()
        chain.Delete()

    # ==================================================
    # Data
    # ==================================================
    hist_data = TH1D('hist_data', stk_title, nbins_h, xmin_h, xmax_h) 
    chain = TChain("Events")
    for dataset in dict_data_yml: 
        # hist_data_chad = TH1D('hist_data_chad', stk_title, nbins_h, xmin_h, xmax_h) 
        # chain = TChain("Events")
        # HLT_Cuts = ''
        # if dataset == 'SingleElectron': 
        #     HLT_Cuts = ''
        # if dataset == 'SingleMuon': 
        #     HLT_Cuts = ' && HLT_Ele35_WPTight_Gsf==0 && HLT_Ele38_WPTight_Gsf==0 && HLT_Ele40_WPTight_Gsf==0'
        # if dataset == 'DoubleEG': 
        #     HLT_Cuts = ' && HLT_Ele35_WPTight_Gsf==0 && HLT_Ele38_WPTight_Gsf==0 && HLT_Ele40_WPTight_Gsf==0 && HLT_IsoMu27==0'
        # if dataset == 'DoubleMuon': 
        #     HLT_Cuts = ' && HLT_Ele35_WPTight_Gsf==0 && HLT_Ele38_WPTight_Gsf==0 && HLT_Ele40_WPTight_Gsf==0 && HLT_IsoMu27==0 && HLT_DoubleEle33_CaloIdL_MW==0 && HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL==0'
        # if dataset == 'MuonEG': 
        #     # HLT_Cuts = ' && HLT_Ele35_WPTight_Gsf==0 && HLT_Ele38_WPTight_Gsf==0 && HLT_Ele40_WPTight_Gsf==0 && HLT_IsoMu27==0 && HLT_DoubleEle33_CaloIdL_MW==0 && HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL==0 && HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8==0 && HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8==0'
        #     HLT_Cuts = ' && HLT_Ele35_WPTight_Gsf==0 && HLT_Ele38_WPTight_Gsf==0 && HLT_Ele40_WPTight_Gsf==0 && HLT_IsoMu27==0 && HLT_DoubleEle33_CaloIdL_MW==0 && HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL==0 && HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8==0'
        for file in dict_data_yml[dataset]['files']: 
            chain.Add(file)
        # chain.SetProof()
        # selections_data = "(%s)" % (' && '.join(cuts) + HLT_Cuts)
        # chain.Project('hist_data_chad', var_to_draw, selections_data)
        # chain.Reset()
        # chain.Delete()
        # hist_data.Add(hist_data_chad)
        # hist_data_chad.SetName("foo")
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
    gStyle.SetOptStat(0)

    pad_stack = TPad("pad_stack","", 0., 0.3, 1., 1.)
    pad_ratio = TPad("pad_ratio","", 0., 0., 1., 0.3)
    pad_stack.Draw()
    pad_ratio.Draw()

    pad_stack.cd()
    stk.Draw("hist")
    # stk.SetMinimum( 1.E-1 )
    stk.GetXaxis().SetTitleOffset(1.1)
    stk.GetYaxis().SetTitleOffset(1.3)

    lg = pad_stack.BuildLegend(0.15, 0.78, 0.9, 0.89, "", "fNDC")
    lg.SetNColumns(3)
    lg.SetFillStyle(0)
    lg.SetBorderSize(0)
    # pad_stack.SetLogx()
    # stk.GetXaxis().SetMoreLogLabels()
    # stk.GetXaxis().SetNoExponent()
    pad_stack.SetLogy()

    hist_data.SetMarkerStyle(kFullDotLarge)
    hist_data.SetLineColor(kBlack)
    hist_data.Draw("EP Same")
    stk.SetMaximum( max(hist_data.GetMaximum()*10.3, stk.GetMaximum()*10.3) )

    pad_stack.Modified()
    pad_stack.Update()
    xMin = pad_stack.GetFrame().GetX1()
    xMax = pad_stack.GetFrame().GetX2()
    yMin = pad_stack.GetFrame().GetY1()
    yMax = pad_stack.GetFrame().GetY2()
    # axisTop = TGaxis(xMin, yMax, xMax, yMax, 0, 1, 1, "-")
    axisTop = TGaxis(xMin, pow(10,yMax), xMax, pow(10,yMax), 0, 1, 1, "-")
    # axisTop = TGaxis(pow(10,xMin), pow(10,yMax), pow(10,xMax), pow(10,yMax), 0, 1, 1, "-")
    axisTop.ChangeLabel(1,-1,0.03,11,-1,-1, "Preliminary")
    axisTop.ChangeLabel(-1,-1,0.03,31,-1,-1, "%s, %s fb^{-1} (13 TeV)" % (region, intLum))
    axisTop.Draw()

    # oldBottomMargin = pad_stack.GetBottomMargin()
    pad_stack.SetBottomMargin(0.)
    # pad_stack.SetTopMargin(pad_stack.GetTopMargin()/0.7)


    # ==================================================
    # Ratio plot
    # ==================================================


    pad_ratio.cd()
    pad_ratio.SetTopMargin(0.)
    pad_ratio.SetBottomMargin(0.25)

    dataOverSumMC = hist_data.Clone('data_ov_MC')
    sumMCErrors = stk.GetHists()[0].Clone('sumMCErrors_hist')
    sumMCErrors.SetFillColorAlpha(kGray, 0.5)
    sumMCErrors.SetMarkerSize(0)

    if len(stk.GetHists()) > 1: 
        map(sumMCErrors.Add, stk.GetHists()[1:])

    dataOverSumMC.Divide(sumMCErrors)
    for i in range(hist_data.GetNbinsX()+2): 
        dataOverSumMC.SetBinError(i, hist_data.GetBinError(i)/max(hist_data.GetBinContent(i), 1))
        sumMCErrors.SetBinError(i, sumMCErrors.GetBinError(i)/max(sumMCErrors.GetBinContent(i), 1))
        sumMCErrors.SetBinContent(i, 1.)

    stk.GetXaxis().Copy(sumMCErrors.GetXaxis())
    sumMCErrors.GetYaxis().SetTitle('Data / MC')
    sumMCErrors.GetYaxis().CenterTitle()
    # sumMCErrors.GetYaxis().SetRangeUser(0.3, 1.7)
    sumMCErrors.GetYaxis().SetRangeUser(0, 4)
    sumMCErrors.GetYaxis().SetNdivisions(003)
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
    # Save plot
    # ==================================================

    cOutput.SaveAs(save_name)
    line.Delete()
    hist.Delete()
    dataOverSumMC.Delete()
    sumMCErrors.Delete()
    hist_data.Delete()
    pad_stack.Delete()
    pad_ratio.Delete()
    cOutput.Delete()
    prf.Close()
    prf.Delete()
    stk.Delete()


def CR_NR():
    cuts = [
    '(lep_category==2)', 
    # 'ngood_bjets==0',
    'ngood_jets<2',
    'nhad_taus==0',
    'met_pt>100',
    # 'abs(delta_R_ll)<3',
    # 'abs(Z_mass-91.1876)>8',
    # 'abs(delta_phi_ZMet)>2.75',
    # 'abs(delta_phi_j_met)>1.2',
    # 'sca_balance>0.4',
    # 'sca_balance<2.5'
    ]
    nbins_h = 20
    xmin_h = 0
    xmax_h = 200
    var_to_draw = 'Z_mass'
    region = "e mu region"
    save_name = "plots/CR_NR_mZ.pdf"
    Draw_pT(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_WZ():
    cuts = [
    '(lep_category==4 || lep_category==5)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    # 'nhad_taus==0',
    # 'met_pt>50',
    'Z_pt>60',
    # 'abs(delta_R_ll)<1.75',
    'abs(Z_mass-91.1876)<15',
    # 'abs(delta_phi_ZMet)>2.75',
    # 'abs(delta_phi_j_met)>1.2',
    # 'sca_balance>0.4',
    # 'sca_balance<2.5'
    ]
    nbins_h = 30
    xmin_h = 0
    xmax_h = 300
    var_to_draw = 'emulatedMET'
    region = "3l region"
    save_name = "plots/CR_WZ_mZ.pdf"
    Draw_pT(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def CR_ZZ():
    cuts = [
    '(lep_category>=6)', 
    'ngood_bjets==0',
    'ngood_jets<2',
    # 'nhad_taus==0',
    'met_pt>40',
    'Z_pt>60',
    # 'abs(delta_R_ll)<1.75',
    'abs(Z_mass-91.1876)<15',
    # 'abs(delta_phi_ZMet)>2.75',
    # 'abs(delta_phi_j_met)>1.2',
    # 'sca_balance>0.4',
    # 'sca_balance<2.5'
    ]
    nbins_h = 15
    xmin_h = 0
    xmax_h = 300
    var_to_draw = 'emulatedMET'
    region = "4l region"
    save_name = "plots/CR_ZZ_mZ.pdf"
    Draw_pT(cuts, nbins_h, xmin_h, xmax_h, var_to_draw, region, save_name)

def main():
    # CR_WZ()
    # CR_NR()
    CR_ZZ()
    gROOT.Reset()


main()



