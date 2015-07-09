from ROOT import *

if __name__ == "__main__":

    file1 = TFile.Open('output.root')
    phibefore = file1.Get('phi_mu_beforeL1')
    phiafter = file1.Get('phi_mu_afterL1')
    phibefore.Sumw2()
    phiafter.Sumw2()
    effphi = phiafter.Clone()
    effphi.Divide(phibefore)
    c1 = TCanvas()
    effphi.Draw()
    c1.SaveAs('AhmedEffPhi.png')
    
    twodbefore = file1.Get('genparticles_ETA_PHI_beforeL1')
    twodafter = file1.Get('genparticles_ETA_PHI_afterL1')

    twoeff = twodafter.Clone()
    twoeff.Divide(twodbefore)

    twoeff.Draw('COLZ')
    c1.SaveAs('AhmedEff2D.png')
