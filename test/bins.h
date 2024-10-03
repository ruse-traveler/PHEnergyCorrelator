
// Set this to use the reduced 2D binning for Run-8
#define RUN8_BINS

// Run-13 response matrix parameters

const int NPTBINS_RECO_13 = 11;
const double PTBINS_RECO_13[NPTBINS_RECO_13 + 1] =
  {
    5.0,
    7.5,
    10.5,
    15.0,
    22.0,
    31.5,
    46.0,
    66.0,
    95.0,
    138.0,
    200.0,
    255.0
  };

const int NPTBINS_TRUE_13 = 11;
const double PTBINS_TRUE_13[NPTBINS_TRUE_13 + 1] =
  {
    5.0,
    7.5,
    10.5,
    15.0,
    22.0,
    31.5,
    46.0,
    66.0,
    95.0,
    138.0,
    200.0,
    255.0
  };

// Run-12 response matrix parameters
// new binning 1/17/2020

const int NPTBINS_RECO_12 = 19;
const double PTBINS_RECO_12[NPTBINS_RECO_12 + 1] =
  {
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    12.0,
    14.5,
    17.5,
    20.5,
    24.5,
    29.0,
    35.0,
    42.0,
    50.0,
    60.0,
    70.0, 
    80.0,
    90.0,
   100.0
  };

const int NPTBINS_TRUE_12 = 19;
const double PTBINS_TRUE_12[NPTBINS_TRUE_12 + 1] =
  {
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    12.0,
    14.5,
    17.5,
    20.5,
    24.5,
    29.0,
    35.0,
    42.0,
    50.0,
    60.0,
    70.0, 
    80.0,
    90.0,
   100.0
  };


// zg unfolding binning (same for all runs)
  
const int NZGBINS = 8; 
double ZGBINS[NZGBINS+1] = 
  {
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.5
  }; 

// oang unfolding binning (same for all runs)
  
const int NOANGBINS = 14; 
double OANGBINS[NOANGBINS+1] = 
  {
    0.0,
    0.025,
    0.050,
    0.075,
    0.10,
    0.125,
    0.15,
    0.175,
    0.20,
    0.225,
    0.25,
    0.275,
    0.30,
    0.325,
    0.35
  }; 

// FF z unfolding binning (same for all runs)

const int NFFZBINS = 9; 
double FFZBINS[NFFZBINS+1] = 
  {
    0.00,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    1.00
  }; 

const int NTFFZBINS = 9; 
double TFFZBINS[NTFFZBINS+1] = 
  {
    0.00,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    1.00
  }; 


#ifndef RUN8_BINS

// Run-12,13 binning

const int NFFXIBINS = 19; 
double FFXIBINS[NFFXIBINS+1] = 
  {
    0.00,
    0.25,
    0.50,
    0.75,
    1.00,
    1.25,
    1.50,
    1.75,
    2.00,
    2.25,
    2.50,
    2.75,
    3.00,
    3.25,
    3.50,
    3.75,
    4.00,
    4.25,
    4.50,
    4.75
  }; 

const int NTFFXIBINS = 19; 
double TFFXIBINS[NTFFXIBINS+1] = 
  {
    0.00,
    0.25,
    0.50,
    0.75,
    1.00,
    1.25,
    1.50,
    1.75,
    2.00,
    2.25,
    2.50,
    2.75,
    3.00,
    3.25,
    3.50,
    3.75,
    4.00,
    4.25,
    4.50,
    4.75
  }; 

// FF JT unfolding binning

//const int NFFJTBINS = 13; 
//double FFJTBINS[NFFJTBINS+1] = 
//  {
//    0.00,
//    0.10,
//    0.20,
//    0.30,
//    0.40,
//    0.50,
//    0.60,
//    0.70,
//    0.80,
//    1.0,
//    1.1,
//    1.2,
//    1.3,
//    1.4
//  }; 

const int NFFJTBINS = 17; 
double FFJTBINS[NFFJTBINS+1] = 
  {
    0.00,
    0.025,
    0.050,
    0.075,
    0.10,
    0.125,
    0.150,
    0.175,
    0.20,
    0.225,
    0.250,
    0.275,
    0.30,
    0.325,
    0.350,
    0.375,
    0.40,
    0.50,
  }; 

// FF dR unfolding binning 
  
const int NFFDRBINS = 12; 
double FFDRBINS[NFFDRBINS+1] = 
  {
    0.00,
    0.025,
    0.05,
    0.075,
    0.10,
    0.125,
    0.150,
    0.175,
    0.20,
    0.225,
    0.250,
    0.275,
    0.30
  }; 

#else

// Run-8 reduced binning

const int NFFXIBINS = 12; 
double FFXIBINS[NFFXIBINS+1] = 
  {
    0.00,
    0.50,
    1.00,
    1.50,
    2.00,
    2.50,
    3.00,
    3.50,
    3.75,
    4.00,
    4.25,
    4.50,
    4.75
  }; 

const int NTFFXIBINS = 12; 
double TFFXIBINS[NTFFXIBINS+1] = 
  {
    0.00,
    0.50,
    1.00,
    1.50,
    2.00,
    2.50,
    3.00,
    3.50,
    3.75,
    4.00,
    4.25,
    4.50,
    4.75
  }; 

// FF JT unfolding binning

//const int NFFJTBINS = 7; 
//double FFJTBINS[NFFJTBINS+1] = 
//  {
//    0.00,
//    0.20,
//    0.40,
//    0.60,
//    0.80,
//    1.00,
//    1.20,
//    1.40
//  }; 

const int NFFJTBINS = 8; 
double FFJTBINS[NFFJTBINS+1] = 
  {
    0.00,
    0.0125,
    0.025,
    0.0375,
    0.050,
    0.0625,
    0.075,
    0.0875,
    0.10
  }; 

// FF dR unfolding binning 
  
const int NFFDRBINS = 6; 
double FFDRBINS[NFFDRBINS+1] = 
  {
    0.00,
    0.05,
    0.10,
    0.150,
    0.20,
    0.250,
    0.30
  }; 


#endif

// Phi bins for Collins distribution

const int NFFPHIBINS = 8; 
double FFPHIBINS[NFFPHIBINS+1] = 
  {
    0.0,
    0.39269875,
    0.7853975,
    1.17809625,
    1.570795,
    1.96349375,
    2.3561925,
    2.74889125,
    3.14159
  }; 

// Phi bins for Collins distribution

const int NJCBINS = 40; 
double JCBINS[NJCBINS+1] = 
  {
    -2.0,
    -1.9,
    -1.8,
    -1.7,
    -1.6,
    -1.5,
    -1.4,
    -1.3,
    -1.2,
    -1.1,
    -1.0,
    -0.9,
    -0.8,
    -0.7,
    -0.6,
    -0.5,
    -0.4,
    -0.3,
    -0.2,
    -0.1,
    0.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.2,
    1.3,
    1.4,
    1.5,
    1.6,
    1.7,
    1.8,
    1.9,
    2.0
  }; 

