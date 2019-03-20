
import yaml
import numpy as np
from ROOT import *


intLum = 35.9
gROOT.SetBatch(1)

with open('ROOTfiles.yml', 'r') as f_yml: 
   dict_yml = yaml.load(f_yml)

datasets = [
   'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
   'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',
   'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',
   'ZZTo2L2Nu_13TeV_powheg_pythia8',
   'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8'
   ]
ndatasets = len(datasets)

selections_basic = [
   'ngood_jets<2', 
   '(lep_category==1 || lep_category==3)',
   'ngood_bjets==0', 
   'nhad_taus==0']

nbins = 60
xmins = [0,    0,    0,    0,    0.2,  0,    0]
xmaxs = [500,  3.5,  30,   3.5,  2.5,  5,    3.5]
stk_mins = [0.1, 0.01, 10, 1, 0.1, 0.001, 1]
cuts_to_plot = [
   'met_pt>90', 
   'abs(delta_R_ll)<1.75', 
   'abs(Z_mass-91.1876)<13', 
   'abs(delta_phi_ZMet)>2.5', 
   'sca_balance>0.6', 
   'sca_balance<2.5', 
   'abs(delta_phi_j_met)>1.2']


def Draw_Stack(cut_to_plot):
   # ============================================================
   # ============================================================
   prf = TProof.Open("lite://")
   
   cuts_tmp = cuts_to_plot[:]
   ind_element = cuts_tmp.index(cut_to_plot)
   xmin = xmins[ind_element]
   xmax = xmaxs[ind_element]
   cut_var = cuts_tmp.pop(ind_element)
   if '<' in cut_var:
      var_to_plot, cut_value = cut_var.split('<')
      less_cut = True
      snc_title = ";;#diamond Significance Integrated from underflowed bin"
   if '>' in cut_var:
      var_to_plot, cut_value = cut_var.split('>')
      less_cut = False
      snc_title = ";;#diamond Significance Integrated from overflowed bin"
   cut_value = float(cut_value)

   selections = 'puWeight * lumiWeight * (' + ' && '.join(selections_basic + cuts_tmp) + ")"

   hex_colors = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462']

   # ============================================================
   # ============================================================
   hists = [None] * ndatasets
   for i in range(ndatasets):
      dataset = datasets[i]
      chain = TChain("Events") 
      for file in dict_yml[dataset]['files']:
         chain.Add(file)
      chain.SetProof()

      hist_name = dataset.split('_')[0]
      hists[i] = TH1D('tmp_name', ';;Events / Bin', nbins, xmin, xmax) 
      chain.Project('tmp_name', var_to_plot, selections)
      hists[i].SetName(hist_name)
      chain.Reset()
      chain.Delete()

      hists[i].Scale(intLum)
      # hists[i].Scale(intLum, "width")
      hists[i].SetFillColor(TColor.GetColor(hex_colors[i]))
      hists[i].SetLineColor(kBlack)

   # ============================================================
   # ============================================================
   snc = TH1D('snc', snc_title, nbins, xmin, xmax)
   snc.SetMarkerColor(TColor.GetColor(hex_colors[-1]))
   snc.SetMarkerStyle(33)
   snc.SetMarkerSize(0.7)

   for i in range(1, nbins+1):
      integral = []
      for j in range(ndatasets):
         integral.append( hists[j].Integral(0,i) if less_cut else hists[j].Integral(i,nbins+1) )
      sig = integral[ datasets.index('ZZTo2L2Nu_13TeV_powheg_pythia8') ]
      bgr = sum(integral) - sig
      significance = 0 if bgr<=0 else sig / bgr**0.5
      snc.SetBinContent( i, abs(significance) )

   stk = THStack("stk", ";%s;Events / Bin" % var_to_plot)
   for i in range(ndatasets):
      hists[i].Rebin(5)
      stk.Add(hists[i])

   # ============================================================
   # Plotting
   # ============================================================
   cOutput = TCanvas("cOutput", "cOutput", 500, 500)
   cOutput.SetFillStyle(4000)
   gStyle.SetOptStat(0)
   pad1 = TPad("pad1","",0,0,1,1)
   pad2 = TPad("pad2","",0,0,1,1)
   pad2.SetFillStyle(4000)
   pad2.SetFrameFillStyle(0)

   pad1.Draw()
   pad1.cd()
   stk.Draw("hist")
   stk_min = stk_mins[ind_element]
   stk.SetMinimum(stk_min)
   stk_max = stk.GetMaximum()
   stk_min_log = np.log10(stk_min)
   stk_max_log = np.log10(stk_max)
   stk_max_new = pow(10, stk_max_log + (stk_max_log - stk_min_log) * 0.3 )
   stk.SetMaximum(stk_max_new)
   lg = gPad.BuildLegend(0.15, 0.78, 0.9, 0.89, "", "fNDC")
   lg.SetNColumns(3)
   lg.SetFillStyle(0)
   lg.SetBorderSize(0)
   gPad.SetLogy()
   gPad.Update()

   arrOpt = "<|-|" if less_cut else "|-|>" 
   xMin =  gPad.GetFrame().GetX1() if less_cut else cut_value 
   xMax =  cut_value if less_cut else gPad.GetFrame().GetX2() 
   x_frm_Min = pad1.GetFrame().GetX1()
   x_frm_Max = pad1.GetFrame().GetX2()
   y_frm_Min = gPad.GetFrame().GetY1()
   y_frm_Max = gPad.GetFrame().GetY2()
   delta_y = y_frm_Max-y_frm_Min
   
   arCut = TArrow(xMin, pow(10, y_frm_Max-0.18*delta_y), xMax, pow(10, y_frm_Max-0.18*delta_y), 0.05, arrOpt)
   arCut.SetAngle(30)
   arCut.SetLineWidth(1)
   arCut.SetLineColor(kGreen+1)
   arCut.SetFillColor(kGreen+1)
   arCut.Draw()

   lCut = TLine(cut_value, pow(10, y_frm_Min), cut_value, pow(10, y_frm_Max))
   lCut.SetLineColor(kBlack)
   lCut.SetLineWidth(2)
   lCut.Draw()

   axisTop = TGaxis(x_frm_Min, pow(10,y_frm_Max), x_frm_Max, pow(10,y_frm_Max), 0, 1, 1, "-")
   # axisTop = TGaxis(pow(10,x_frm_Min), pow(10,y_frm_Max), pow(10,x_frm_Max), pow(10,y_frm_Max), 0, 1, 1, "-")
   axisTop.ChangeLabel(1,-1,0.03,11,-1,-1, "Preliminary")
   axisTop.ChangeLabel(-1,-1,0.03,31,-1,-1, "%s fb^{-1} (13 TeV)" % intLum)
   axisTop.Draw()

   pad2.Draw()
   pad2.cd()
   snc.GetYaxis().SetTitleOffset(1.4)
   snc_min = snc.GetMinimum()
   snc_max = snc.GetMaximum()
   snc.SetMaximum( snc_max + 0.3 * (snc_max - snc_min) )
   snc.GetYaxis().SetAxisColor(TColor.GetColor(hex_colors[-1]))
   snc.GetYaxis().SetLabelColor(TColor.GetColor(hex_colors[-1]))
   snc.GetYaxis().SetTitleColor(TColor.GetColor(hex_colors[-1]))
   snc.Draw("Y+ hist p")
   # cOutput.SaveAs( "plots/Nm1WithSig_%d.pdf" % (ind_element+1) )

   # ============================================================
   # cleaning things up
   # ============================================================
   snc.Delete()
   for hist in hists:
      hist.Delete()
   stk.Delete()
   prf.Close()
   prf.Delete()

def Draw_Nm1():
   for cut_to_plot in cuts_to_plot:
      Draw_Stack(cut_to_plot)


def main():
   Draw_Nm1()
   gROOT.Reset()


if __name__ == "__main__":
   main()
