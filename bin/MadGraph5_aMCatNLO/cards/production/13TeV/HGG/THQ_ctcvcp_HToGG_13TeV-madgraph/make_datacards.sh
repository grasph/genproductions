#!/bin/bash

for MH in 120 123 124 125 126 127 130; do

d=THQ_ctcvcp_HToGG_M${MH}_TuneCP5_13TeV-madgraph-pythia8

mkdir "$d"

cat > $d/thq_4f_ckm_LO_ctcvcp_MH${MH}_customizecards.dat <<EOF
#put card customizations here (change top and higgs mass for example)
set param_card mass 6 172.5
set param_card mass 25 ${MH}
set param_card yukawa 6 -172.5
EOF

cat > $d/thq_4f_ckm_LO_ctcvcp_MH${MH}_extramodels.dat <<EOF
#higgs characterization NLO model version 1.3 from http://feynrules.irmp.ucl.ac.be/wiki/HiggsCharacterisation
HC_NLO_X0_UFO-v1.3.zip

EOF

cat > $d/thq_4f_ckm_LO_ctcvcp_MH${MH}_proc_card.dat <<EOF
set group_subprocesses Auto
set ignore_six_quark_processes False
set loop_optimized_output True
set complex_mass_scheme False
define l+ = e+ mu+ ta+
define l- = e- mu- ta-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define wdec = u c d s u~ c~ d~ s~ l+ l- vl vl~
import model HC_NLO_X0_UFO
generate p p > x0 t~ b j \$\$ w+ w-, (t~ > b~ w-, w- > wdec wdec)
add process p p > x0 t b~ j \$\$ w+ w-, (t > b w+, w+ > wdec wdec)
output thq_4f_ckm_LO_ctcvcp_MH${MH} -nojpeg
EOF

cat > $d/thq_4f_ckm_LO_ctcvcp_MH${MH}_reweight_card.dat <<EOF
change rwgt_dir rwgt
launch
set kHtt 3.0
launch
set kHtt 2.0
launch
set kHtt 1.5
launch
set kHtt 1.25
launch
set kHtt 0.75 
launch
set kHtt 0.5
launch
set kHtt 0.25
launch
set kHtt 0.0001
launch
set kHtt -0.25
launch
set kHtt -0.5
launch
set kHtt -0.75
launch
set kHtt -1.0
launch
set kHtt -1.25
launch
set kHtt -1.5
launch
set kHtt -2.0
launch
set kHtt -3.0
launch
set kHtt 3.0
set kSM 1.5
launch
set kHtt 2.0
set kSM 1.5
launch
set kHtt 1.5
set kSM 1.5
launch
set kHtt 1.25
set kSM 1.5
launch
set kHtt 1.0
set kSM 1.5
launch
set kHtt 0.75
set kSM 1.5
launch
set kHtt 0.5 
set kSM 1.5
launch
set kHtt 0.25
set kSM 1.5
launch
set kHtt 0.0001
set kSM 1.5
launch
set kHtt -0.25
set kSM 1.5
launch
set kHtt -0.5
set kSM 1.5
launch
set kHtt -0.75
set kSM 1.5
launch
set kHtt -1.0
set kSM 1.5
launch
set kHtt -1.25
set kSM 1.5
launch
set kHtt -1.5
set kSM 1.5
launch
set kHtt -2.0
set kSM 1.5
launch
set kHtt -3.0
set kSM 1.5
launch
set kHtt 3.0
set kSM 0.5
launch
set kHtt 2.0
set kSM 0.5
launch
set kHtt 1.5
set kSM 0.5
launch
set kHtt 1.25
set kSM 0.5
launch
set kHtt 1.0
set kSM 0.5
launch
set kHtt 0.75
set kSM 0.5
launch
set kHtt 0.5 
set kSM 0.5
launch
set kHtt 0.25
set kSM 0.5
launch
set kHtt 0.0001
set kSM 0.5
launch
set kHtt -0.25
set kSM 0.5
launch
set kHtt -0.5
set kSM 0.5
launch
set kHtt -0.75
set kSM 0.5
launch
set kHtt -1.0
set kSM 0.5
launch
set kHtt -1.25
set kSM 0.5
launch
set kHtt -1.5
set kSM 0.5
launch
set kHtt -2.0
set kSM 0.5
launch
set kHtt -3.0
set kSM 0.5
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.9
set kSM -1.1111
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.8
set kSM -1.25
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.7
set kSM -1.42857
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.6
set kSM -1.6667
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.5
set kSM -2.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.4
set kSM -2.5
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.3
set kSM -3.333
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.2
set kSM -5.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.1
set kSM -10.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa -0.0001
set kSM -10000
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.1
set kSM 10.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.2
set kSM 5.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.3
set kSM 3.333
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.4
set kSM 2.5
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.5
set kSM 2.0
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.6
set kSM 1.6667
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.7
set kSM 1.42857
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.8
set kSM 1.25
launch
set ymt 172.5
set kAtt 0.666667
set kHll 0.000000
set kAll 0.000000
set kHaa 0.000000
set kAaa 0.000000
set kHza 0.000000
set kAza 0.000000
set cosa 0.9
set kSM 1.1111
EOF

cat > $d/thq_4f_ckm_LO_ctcvcp_MH${MH}_run_card.dat <<EOF
*******************                                                 
# Running parameters
#*******************                                                 
#                                                                    
#*********************************************************************
# Tag name for the run (one word)                                    *
#*********************************************************************
  tag_1	=  run_tag  ! name of the run 
#*********************************************************************
# Run to generate the grid pack                                      *
#*********************************************************************
  False	=  gridpack   !True = setting up the grid pack
#*********************************************************************
# Number of events and rnd seed                                      *
# Warning: Do not generate more than 1M events in a single run       *
# If you want to run Pythia, avoid more than 50k events in a run.    *
#*********************************************************************
  1000	=  nevents  ! Number of unweighted events requested 
  0	=  iseed    ! rnd seed (0=assigned automatically=default))
#*********************************************************************
# Collider type and energy                                           *
# lpp: 0=No PDF, 1=proton, -1=antiproton, 2=photon from proton,      *
#                                         3=photon from electron     *
#*********************************************************************
  1	=  lpp1     ! beam 1 type 
  1	=  lpp2     ! beam 2 type
  6500.0	=  ebeam1   ! beam 1 total energy in GeV
  6500.0	=  ebeam2   ! beam 2 total energy in GeV
#*********************************************************************
# Beam polarization from -100 (left-handed) to 100 (right-handed)    *
#*********************************************************************
  0.0	=  polbeam1  ! beam polarization for beam 1
  0.0	=  polbeam2  ! beam polarization for beam 2
#*********************************************************************
# PDF CHOICE: this automatically fixes also alpha_s and its evol.    *
#*********************************************************************
  'lhapdf'    = pdlabel     ! PDF set
  \$DEFAULT_PDF_SETS = lhaid
  \$DEFAULT_PDF_MEMBERS = reweight_PDF     ! if pdlabel=lhapdf, this is the lhapdf number
#*********************************************************************
# Renormalization and factorization scales                           *
#*********************************************************************
  False	=  fixed_ren_scale   ! if .true. use fixed ren scale
  False	=  fixed_fac_scale   ! if .true. use fixed fac scale
  91.188	=  scale             ! fixed ren scale
  91.188	=  dsqrt_q2fact1     ! fixed fact scale for pdf1
  91.188	=  dsqrt_q2fact2     ! fixed fact scale for pdf2
  -1	=  dynamical_scale_choice  ! Select one of the preselect dynamical choice
  0.33 =  scalefact         ! scale factor for event-by-event scales
#*********************************************************************
# Matching - Warning! ickkw > 1 is still beta
#*********************************************************************
  0	=  ickkw             ! 0 no matching, 1 MLM, 2 CKKW matching
  1	=  highestmult       ! for ickkw=2, highest mult group
  1	=  ktscheme          ! for ickkw=1, 1 Durham kT, 2 Pythia pTE
  1	=  alpsfact          ! scale factor for QCD emission vx
  False	=  chcluster         ! cluster only according to channel diag
  False=  pdfwgt            ! for ickkw=1, perform pdf reweighting
  5	=  asrwgtflavor      ! highest quark flavor for a_s reweight
  True	=  clusinfo          ! include clustering tag in output
  3.0	=  lhe_version        ! Change the way clustering information pass to shower.        
#*********************************************************************
#**********************************************************
#
#**********************************************************
# Automatic ptj and mjj cuts if xqcut > 0
# (turn off for VBF and single top processes)
#**********************************************************
  False	=  auto_ptj_mjj   ! Automatic setting of ptj and mjj
#**********************************************************
#                                                                    
#**********************************
# BW cutoff (M+/-bwcutoff*Gamma)
#**********************************
  15.0	=  bwcutoff       ! (M+/-bwcutoff*Gamma)
#**********************************************************
# Apply pt/E/eta/dr/mij cuts on decay products or not
# (note that etmiss/ptll/ptheavy/ht/sorted cuts always apply)
#**********************************************************
  False	=  cut_decays     ! Cut decay products 
#*************************************************************
# Number of helicities to sum per event (0 = all helicities)
# 0 gives more stable result, but longer run time (needed for
# long decay chains e.g.).
# Use >=2 if most helicities contribute, e.g. pure QCD.
#*************************************************************
  0	=  nhel           ! Number of helicities used per event
#*******************                                                 
# Standard Cuts
#*******************                                                 
#                                                                    
#*********************************************************************
# Minimum and maximum pt's (for max, -1 means no cut)                *
#*********************************************************************
  0.1	=  ptj        ! minimum pt for the jets 
  0.1	=  ptb        ! minimum pt for the b 
  0.1	=  pta        ! minimum pt for the photons 
  0.1	=  ptl        ! minimum pt for the charged leptons 
  0.0	=  misset     ! minimum missing Et (sum of neutrino's momenta)
  0.0	=  ptheavy    ! minimum pt for one heavy final state
  -1.0	=  ptjmax     ! maximum pt for the jets
  -1.0	=  ptbmax     ! maximum pt for the b
  -1.0	=  ptamax     ! maximum pt for the photons
  -1.0	=  ptlmax     ! maximum pt for the charged leptons
  -1.0	=  missetmax  ! maximum missing Et (sum of neutrino's momenta)
#*********************************************************************
# Minimum and maximum E's (in the center of mass frame)              *
#*********************************************************************
  0.0	=  ej      ! minimum E for the jets 
  0.0	=  eb      ! minimum E for the b 
  0.0	=  ea      ! minimum E for the photons 
  0.0	=  el      ! minimum E for the charged leptons 
  -1.0	=  ejmax  ! maximum E for the jets
  -1.0	=  ebmax  ! maximum E for the b
  -1.0	=  eamax  ! maximum E for the photons
  -1.0	=  elmax  ! maximum E for the charged leptons
#*********************************************************************
# Maximum and minimum absolute rapidity (for max, -1 means no cut)   *
#*********************************************************************
  -1.0	=  etaj     ! max rap for the jets 
  -1.0	=  etab     ! max rap for the b
  -1.0	=  etaa     ! max rap for the photons 
  -1.0	=  etal     ! max rap for the charged leptons 
  0.0	=  etajmin  ! min rap for the jets
  0.0	=  etabmin  ! min rap for the b
  0.0	=  etaamin  ! min rap for the photons
  0.0	=  etalmin  ! main rap for the charged leptons
#*********************************************************************
# Minimum and maximum DeltaR distance                                *
#*********************************************************************
  0.01	=  drjj     ! min distance between jets 
  0.0	=  drbb     ! min distance between b's 
  0.01	=  drll     ! min distance between leptons 
  0.01	=  draa     ! min distance between gammas 
  0.0	=  drbj     ! min distance between b and jet 
  0.01	=  draj     ! min distance between gamma and jet 
  0.01	=  drjl     ! min distance between jet and lepton 
  0.0	=  drab     ! min distance between gamma and b 
  0.0	=  drbl     ! min distance between b and lepton 
  0.01	=  dral     ! min distance between gamma and lepton 
  -1.0	=  drjjmax  ! max distance between jets
  -1.0	=  drbbmax  ! max distance between b's
  -1.0	=  drllmax  ! max distance between leptons
  -1.0	=  draamax  ! max distance between gammas
  -1.0	=  drbjmax  ! max distance between b and jet
  -1.0	=  drajmax  ! max distance between gamma and jet
  -1.0	=  drjlmax  ! max distance between jet and lepton
  -1.0	=  drabmax  ! max distance between gamma and b
  -1.0	=  drblmax  ! max distance between b and lepton
  -1.0	=  dralmax  ! maxdistance between gamma and lepton
#*********************************************************************
# Minimum and maximum invariant mass for pairs                       *
# WARNING: for four lepton final state mmll cut require to have      *
#          different lepton masses for each flavor!                  *           
#*********************************************************************
  0.0	=  mmjj     ! min invariant mass of a jet pair 
  0.0	=  mmbb     ! min invariant mass of a b pair 
  0.0	=  mmaa     ! min invariant mass of gamma gamma pair
  0.0	=  mmll     ! min invariant mass of l+l- (same flavour) lepton pair
  -1.0	=  mmjjmax  ! max invariant mass of a jet pair
  -1.0	=  mmbbmax  ! max invariant mass of a b pair
  -1.0	=  mmaamax  ! max invariant mass of gamma gamma pair
  -1.0	=  mmllmax  ! max invariant mass of l+l- (same flavour) lepton pair
#*********************************************************************
# Minimum and maximum invariant mass for all letpons                 *
#*********************************************************************
  0.0	=  mmnl     ! min invariant mass for all letpons (l+- and vl) 
  -1.0	=  mmnlmax  ! max invariant mass for all letpons (l+- and vl) 
#*********************************************************************
# Minimum and maximum pt for 4-momenta sum of leptons                *
#*********************************************************************
  0.0	=  ptllmin   ! Minimum pt for 4-momenta sum of leptons(l and vl)
  -1.0	=  ptllmax   ! Maximum pt for 4-momenta sum of leptons(l and vl)
#*********************************************************************
# Inclusive cuts                                                     *
#*********************************************************************
  0.0	=  xptj  ! minimum pt for at least one jet  
  0.0	=  xptb  ! minimum pt for at least one b 
  0.0	=  xpta  ! minimum pt for at least one photon 
  0.0	=  xptl  ! minimum pt for at least one charged lepton 
#*********************************************************************
# Control the pt's of the jets sorted by pt                          *
#*********************************************************************
  0.0	=  ptj1min  ! minimum pt for the leading jet in pt
  0.0	=  ptj2min  ! minimum pt for the second jet in pt
  0.0	=  ptj3min  ! minimum pt for the third jet in pt
  0.0	=  ptj4min  ! minimum pt for the fourth jet in pt
  -1.0	=  ptj1max  ! maximum pt for the leading jet in pt 
  -1.0	=  ptj2max  ! maximum pt for the second jet in pt
  -1.0	=  ptj3max  ! maximum pt for the third jet in pt
  -1.0	=  ptj4max  ! maximum pt for the fourth jet in pt
  0	=  cutuse   ! reject event if fails any (0) / all (1) jet pt cuts
#*********************************************************************
# Control the pt's of leptons sorted by pt                           *
#*********************************************************************
  0.0	=  ptl1min  ! minimum pt for the leading lepton in pt
  0.0	=  ptl2min  ! minimum pt for the second lepton in pt
  0.0	=  ptl3min  ! minimum pt for the third lepton in pt
  0.0	=  ptl4min  ! minimum pt for the fourth lepton in pt
  -1.0	=  ptl1max  ! maximum pt for the leading lepton in pt 
  -1.0	=  ptl2max  ! maximum pt for the second lepton in pt
  -1.0	=  ptl3max  ! maximum pt for the third lepton in pt
  -1.0	=  ptl4max  ! maximum pt for the fourth lepton in pt
#*********************************************************************
# Control the Ht(k)=Sum of k leading jets                            *
#*********************************************************************
  0.0	=  htjmin  ! minimum jet HT=Sum(jet pt)
  -1.0	=  htjmax  ! maximum jet HT=Sum(jet pt)
  0.0	=  ihtmin   !inclusive Ht for all partons (including b)
  -1.0	=  ihtmax   !inclusive Ht for all partons (including b)
  0.0	=  ht2min  ! minimum Ht for the two leading jets
  0.0	=  ht3min  ! minimum Ht for the three leading jets
  0.0	=  ht4min  ! minimum Ht for the four leading jets
  -1.0	=  ht2max  ! maximum Ht for the two leading jets
  -1.0	=  ht3max  ! maximum Ht for the three leading jets
  -1.0	=  ht4max  ! maximum Ht for the four leading jets
#***********************************************************************
# Photon-isolation cuts, according to hep-ph/9801442                   *
# When ptgmin=0, all the other parameters are ignored                  *
# When ptgmin>0, pta and draj are not going to be used                 *
#***********************************************************************
  0.0	=  ptgmin  ! Min photon transverse momentum
  0.4	=  R0gamma  ! Radius of isolation code
  1.0	=  xn  ! n parameter of eq.(3.4) in hep-ph/9801442
  1.0	=  epsgamma  ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442
  True	=  isoEM  ! isolate photons from EM energy (photons and leptons)
#*********************************************************************
# WBF cuts                                                           *
#*********************************************************************
  0.0	=  xetamin  ! minimum rapidity for two jets in the WBF case  
  0.0	=  deltaeta  ! minimum rapidity for two jets in the WBF case 
#*********************************************************************
# KT DURHAM CUT                                                      *
#*********************************************************************
  -1.0	=   ktdurham        
   0.4	=   dparameter 
 #*********************************************************************
# maximal pdg code for quark to be considered as a light jet         *
# (otherwise b cuts are applied)                                     *
#*********************************************************************
  4	=  maxjetflavor     ! Maximum jet pdg code
#*********************************************************************
# Jet measure cuts                                                   *
#*********************************************************************
  0.0	=  xqcut    ! minimum kt jet measure between partons
#*********************************************************************
EOF

done