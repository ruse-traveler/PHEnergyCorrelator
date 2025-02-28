

void MinimalCalc() {

  // announce staart
  std::cout << "\n  Starting minimal calc..." << std::endl;

  // create output file
  TFile* fOutput = new TFile("minimalCalc_step2_randomSpinXY.nIter10k.d28m2y2025.root", "recreate");
  std::cout << "    Created output file." << std::endl;

  // initialize rng
  TDatime*  time  = new TDatime();
  TRandom3* rando = new TRandom3();
  rando -> SetSeed(time -> Get());
  std::cout << "    Initialized RNG." << std::endl;

  // histogram binning
  const int   nAngBins  = 100;
  const int   nMagBins  = 100;
  const float xAngStart = -TMath::TwoPi();
  const float xAngStop  = TMath::TwoPi(); 
  const float xMagStart = -5.;
  const float xMagStop  = 5.;

  // create histograms
  TH1D* hPhiSpin        = new TH1D("hPhiSpin", "", nAngBins, xAngStart, xAngStop);
  TH2D* hPhiSpinVsDot   = new TH2D("hPhiSpinVsDot", "", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHadron      = new TH1D("hPhiHadron", "", nAngBins, xAngStart, xAngStop);
  TH2D* hPhiHadronVsDot = new TH2D("hPhiHadronVsDot", "", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);
  TH1D* hPhicollins     = new TH1D("hPhiCollins", "", nAngBins, xAngStart, xAngStop);
  std::cout << "    Created histograms." << std::endl;

  // set beam direction
  TVector3 beam(0.0, 0.0, 1.0);

  // announce mc loop starting
  const std::size_t nIter = 10000;
  std::cout << "    Staring MC loop: " << nIter << " iterations." << std::endl;

  // mc loop
  for (std::size_t iIter = 0; iIter < nIter; ++iIter) {

    // random spin
    const double xs = rando -> Uniform(-1.0, 1.0);
    const double ys = rando -> Uniform(-1.0, 1.0);
    TVector3 spin_noNorm(xs, ys, 0.0);

    // random jet axis
    const double xj = rando -> Uniform(-1.0, 1.0);
    const double yj = rando -> Uniform(-1.0, 1.0);
    const double zj = rando -> Uniform(-1.0, 1.0);
    TVector3 jet_noNorm(xj, yj, zj);

    // random hadron direction
    const double xh = rando -> Uniform(-1.0, 1.0);
    const double yh = rando -> Uniform(-1.0, 1.0);
    const double zh = rando -> Uniform(-1.0, 1.0);
    TVector3 hadron_noNorm(xh, yh, zh);

    // normalize spin/jet/hadron
    TVector3 spin   = spin_noNorm.Unit();
    TVector3 jet    = jet_noNorm.Unit();
    TVector3 hadron = hadron_noNorm.Unit();

    // get spin-beam, jet-beam, hadron-jet normals
    TVector3 nSpinBeam  = ( beam.Cross(spin) ).Unit();
    TVector3 nJetBeam   = ( beam.Cross(jet) ).Unit();
    TVector3 nHadronJet = ( jet.Cross(hadron) ).Unit();

    // ------------------------------------------------------------------------

    // get normal to spin & jet-beam
    TVector3 nJetBeamSpin = nJetBeam.Cross(nSpinBeam);

    // calculate dot products
    const double dotSpinJet = nJetBeam.Dot( nSpinBeam );
    const double dotSpin    = nJetBeam.Dot( spin );

    // calculate phi spin
    double phiSpin = atan2( nJetBeamSpin.Mag(), dotSpinJet);
    if (dotSpin < 0.0) phiSpin = TMath::TwoPi() - phiSpin;

    // fill histograms
    hPhiSpin      -> Fill(phiSpin);
    hPhiSpinVsDot -> Fill(dotSpin, phiSpin);

    // ------------------------------------------------------------------------

    // get normal to hadron-jet & jet-beam
    TVector3 nJetBeamHadron = nJetBeam.Cross(nHadronJet);

    // calculate dot products
    const double dotHadJet = nJetBeam.Dot( nHadronJet );
    const double dotHadron = nJetBeam.Dot( hadron );

    // calculate phi hadron
    double phiHadron = atan2( nJetBeamHadron.Mag(), dotHadJet );
    if (dotHadron < 0.0) phiHadron = TMath::TwoPi() - phiHadron;

    // fill histograms
    hPhiHadron      -> Fill(phiHadron);
    hPhiHadronVsDot -> Fill(dotHadron, phiHadron);

    // ------------------------------------------------------------------------

    // calculate phi collins
    double phiCollins = phiSpin - phiHadron;
    if (phiCollins < 0.0)            phiCollins += TMath::TwoPi();
    if (phiCollins > TMath::TwoPi()) phiCollins -= TMath::TwoPi();

    // fill histogram
    hPhiCollins -> Fill(phiCollins);

  }
  std::cout << "    MC loop finished." << std::endl;

  // save histograms
  fOutput         -> cd();
  hPhiSpin        -> Write();
  hPhiSpinVsDot   -> Write();
  hPhiHadron      -> Write();
  hPhiHadronVsDot -> Write();
  hPhicollins     -> Write();
  std::cout << "    Saved histograms." << std::endl;

  // announce end
  std::cout << "  Ending minimal calc.\n" << std::endl;
  return;

}

// end ========================================================================
