from ROOT import *

if __name__ == "__main__":

    file1 = TFile.Open('output.root')
    etabefore = file1.Get('eta_mu_beforeL1')
    etaafter = file1.Get('eta_mu_afterL1')
    etabefore.Sumw2()
    etaafter.Sumw2()
    effeta = etaafter.Clone()
    effeta.Divide(etabefore)
    c1 = TCanvas()
    effeta.Draw()
    c1.SaveAs('AhmedEffeta.png')
    
    twodbefore = file1.Get('genparticles_ETA_PHI_beforeL1')
    twodafter = file1.Get('genparticles_ETA_PHI_afterL1')

    twoeff = twodafter.Clone()
    twoeff.Divide(twodbefore)

    twoeff.Draw('COLZ')
    c1.SaveAs('AhmedEff2D.png')
