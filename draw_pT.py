
import yaml
import os
from array import *
from ROOT import *


intLum = 150
gROOT.SetBatch(1)

with open('../../python/postprocessing/monoZ/ROOTfiles.yml', 'r') as f_yml:
    dict_yml = yaml.load(f_yml)


def Draw_pT():
   datasets = [
   'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
   'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',
   'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
   'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',
   'ZZTo2L2Nu_13TeV_powheg_pythia8'
   ]
   hex_cols = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3']
   ndatasets = len(datasets)

   stk = THStack("stk", ";Dilepton pT [GeV];Events / GeV")
   bins_Det = array('d', [40, 60, 80, 100, 200, 400, 1000])
   # bins_Det = array('d', [    50.00, 55.62, 61.25, 66.87, 72.49, 78.12, 84.88, 91.65, 
   #    98.42, 105.19, 111.96, 118.73, 125.50, 134.55, 143.60, 152.66, 161.71, 170.77, 179.82, 191.64, 
   #    203.46, 215.28, 227.10, 241.96, 256.82, 271.69, 286.55, 305.06, 323.58, 342.10, 367.25, 392.40, 
   #    419.66, 446.92, 476.25, 509.53, 542.81, 580.49, 618.98, 657.47, 698.26, 741.03, 790.66, 
   #    845.02, 899.31, 954.91, 1011.28, 1067.66])
   nbins_Det = len(bins_Det)-1

   cuts = [
   '(lep_category==1 || lep_category==3)', 
   'ngood_bjets==0',
   'ngood_jets<2',
   'nhad_taus==0',
   'met_pt>90',
   'abs(delta_R_ll)<1.75',
   'abs(Z_mass-91.1876)<8',
   'abs(delta_phi_ZMet)>2.75',
   'abs(delta_phi_j_met)>1.2',
   'sca_balance>0.4',
   'sca_balance<2.5'
   ]
   selections = "puWeight * lumiWeight * (%s)" % ' && '.join(cuts)

   prf = TProof.Open("lite://")
   for i in range(ndatasets):
      dataset = datasets[i]
      hist_name = dataset.split('_')[0]
      hist = TH1D(hist_name, ';pT [GeV];events / GeV', nbins_Det, bins_Det) 
      chain = TChain("Events")
      for file in dict_yml[dataset]['files']:
         chain.Add(file)
      chain.SetProof()
      chain.Project(hist_name, "Z_pt", selections)
      hist.Scale(intLum, "width")
      hist.SetFillColor( TColor.GetColor(hex_cols[i]) )
      hist.SetLineColor(kBlack)
      stk.Add(hist)
      chain.Reset()
      chain.Delete()

   # ==================================================
   # Plotting
   # ==================================================

   cOutput = TCanvas("cOutput", "cOutput", 500, 500)
   cOutput.SetFillStyle(4000)
   gStyle.SetOptStat(0)

   pad_out = TPad("pad_out","",0,0,1,1)

   pad_out.Draw()
   pad_out.cd()
   stk.Draw("hist")
   stk.SetMaximum( 10 * stk.GetMaximum() )
   stk.SetMinimum( 1.E-3 )
   stk.GetXaxis().SetTitleOffset(1.1)
   stk.GetYaxis().SetTitleOffset(1.3)

   lg = pad_out.BuildLegend(0.15, 0.78, 0.9, 0.89, "", "fNDC")
   lg.SetNColumns(3)
   lg.SetFillStyle(0)
   lg.SetBorderSize(0)
   pad_out.SetLogx()
   stk.GetXaxis().SetMoreLogLabels()
   stk.GetXaxis().SetNoExponent()
   pad_out.SetLogy()

   pad_out.Modified()
   pad_out.Update()
   xMin = pad_out.GetFrame().GetX1()
   xMax = pad_out.GetFrame().GetX2()
   yMin = pad_out.GetFrame().GetY1()
   yMax = pad_out.GetFrame().GetY2()
   # axisTop = TGaxis(xMin, pow(10,yMax), xMax, pow(10,yMax), 0, 1, 1, "-")
   axisTop = TGaxis(pow(10,xMin), pow(10,yMax), pow(10,xMax), pow(10,yMax), 0, 1, 1, "-")
   axisTop.ChangeLabel(1,-1,0.03,11,-1,-1, "Preliminary")
   axisTop.ChangeLabel(-1,-1,0.03,31,-1,-1, "%s fb^{-1} (13 TeV)" % intLum)
   axisTop.Draw()

   cOutput.SaveAs("plots/Rec_bined.pdf")
   stk.Delete()
   pad_out.Delete()
   prf.Close()
   prf.Delete()


def main():
   Draw_pT()
   gROOT.Reset()


main()

