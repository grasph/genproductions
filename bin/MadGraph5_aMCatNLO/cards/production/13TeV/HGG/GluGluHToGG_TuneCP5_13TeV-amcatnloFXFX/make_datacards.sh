#!/bin/bash

for MH in 60 65 70 75 80 85 90 95 100 105 110 115 120 123 124 125 126 127 130; do
d=GluGluHToGG_M${MH}_TuneCP5_13TeV-amcatnloFXFX-pythia8
mkdir "$d"

cat > $d/ggh012j_5f_NLO_FXFX_${MH}_customizecards.dat <<EOF
set param_card mass 25 ${MH}.
EOF

cat > $d/ggh012j_5f_NLO_FXFX_${MH}_extramodels.dat <<EOF
HC_NLO_X0_UFO-v1.3.zip
EOF

cat > $d/ggh012j_5f_NLO_FXFX_${MH}_proc_card.dat <<EOF
import model HC_NLO_X0_UFO-heft

generate p p > x0 / t [QCD] @0
add process p p > x0 j / t [QCD] @1
add process p p > x0 j j / t [QCD] @2

output ggh012j_5f_NLO_FXFX_${MH} -nojpeg
EOF

cat > $d/ggh012j_5f_NLO_FXFX_${MH}_run_card.dat <<EOF
#*******************                                                 
# Running parameters
#*******************                                                 
#
#***********************************************************************
# Tag name for the run (one word)                                      *
#***********************************************************************
  tag_1     = run_tag ! name of the run 
#***********************************************************************
# Number of events (and their normalization) and the required          *
# (relative) accuracy on the Xsec.                                     *
# These values are ignored for fixed order runs                        *
#***********************************************************************
  1500 = nevents ! Number of unweighted events requested 
 0.001 = req_acc ! Required accuracy (-1=auto determined from nevents)
    20 = nevt_job! Max number of events per job in event generation. 
                 !  (-1= no split).
average = event_norm ! Normalize events to sum or average to the X sect.
#***********************************************************************
# Number of points per itegration channel (ignored for aMC@NLO runs)   *
#***********************************************************************
 0.01   = req_acc_FO       ! Required accuracy (-1=ignored, and use the 
 	                   ! number of points and iter. below)
# These numbers are ignored except if req_acc_FO is equal to -1
 5000   = npoints_FO_grid  ! number of points to setup grids
 4      = niters_FO_grid   ! number of iter. to setup grids
 10000  = npoints_FO       ! number of points to compute Xsec
 6      = niters_FO        ! number of iter. to compute Xsec
#***********************************************************************
# Random number seed                                                   *
#***********************************************************************
     0    = iseed       ! rnd seed (0=assigned automatically=default))
#***********************************************************************
# Collider type and energy                                             *
#***********************************************************************
    1   = lpp1    ! beam 1 type (0 = no PDF)
    1   = lpp2    ! beam 2 type (0 = no PDF)
 6500   = ebeam1  ! beam 1 energy in GeV
 6500   = ebeam2  ! beam 2 energy in GeV
#***********************************************************************
# PDF choice: this automatically fixes also alpha_s(MZ) and its evol.  *
#***********************************************************************
 lhapdf    = pdlabel   ! PDF set                                     
 \$DEFAULT_PDF_SETS = lhaid     ! if pdlabel=lhapdf, this is the lhapdf number
#***********************************************************************
# Include the NLO Monte Carlo subtr. terms for the following parton    *
# shower (HERWIG6 | HERWIGPP | PYTHIA6Q | PYTHIA6PT | PYTHIA8)         *
# WARNING: PYTHIA6PT works only for processes without FSR!!!!          *
#***********************************************************************
  PYTHIA8   = parton_shower
#***********************************************************************
# Renormalization and factorization scales                             *
# (Default functional form for the non-fixed scales is the sum of      *
# the transverse masses of all final state particles and partons. This *
# can be changed in SubProcesses/set_scales.f)                         *
#***********************************************************************
 F        = fixed_ren_scale  ! if .true. use fixed ren scale
 F        = fixed_fac_scale  ! if .true. use fixed fac scale
 91.188   = muR_ref_fixed    ! fixed ren reference scale 
 91.188   = muF1_ref_fixed   ! fixed fact reference scale for pdf1
 91.188   = muF2_ref_fixed   ! fixed fact reference scale for pdf2
#***********************************************************************
# Renormalization and factorization scales (advanced and NLO options)  *
#***********************************************************************
 F        = fixed_QES_scale  ! if .true. use fixed Ellis-Sexton scale
 91.188   = QES_ref_fixed    ! fixed Ellis-Sexton reference scale
 1        = muR_over_ref     ! ratio of current muR over reference muR
 1        = muF1_over_ref    ! ratio of current muF1 over reference muF1
 1        = muF2_over_ref    ! ratio of current muF2 over reference muF2
 1        = QES_over_ref     ! ratio of current QES over reference QES
#*********************************************************************** 
# Reweight flags to get scale dependence and PDF uncertainty           *
# For scale dependence: factor rw_scale_up/down around central scale   *
# For PDF uncertainty: use LHAPDF with supported set                   *
#***********************************************************************
 .true.   = reweight_scale   ! reweight to get scale dependence
  0.5     = rw_Rscale_down   ! lower bound for ren scale variations
  2.0     = rw_Rscale_up     ! upper bound for ren scale variations
  0.5     = rw_Fscale_down   ! lower bound for fact scale variations
  2.0     = rw_Fscale_up     ! upper bound for fact scale variations
 \$DEFAULT_PDF_MEMBERS = reweight_PDF     ! reweight to get PDF uncertainty
      ! First of the error PDF sets 
      ! Last of the error PDF sets
#***********************************************************************
# Merging - WARNING! Applies merging only at the hard-event level.     *
# After showering an MLM-type merging should be applied as well.       *
# See http://amcatnlo.cern.ch/FxFx_merging.htm for more details.       *
#***********************************************************************
 3        = ickkw            ! 0 no merging, 3 FxFx merging
#***********************************************************************
#
#***********************************************************************
# BW cutoff (M+/-bwcutoff*Gamma)                                       *
#***********************************************************************
 15  = bwcutoff
#***********************************************************************
# Cuts on the jets                                                     *
# Jet clustering is performed by FastJet.
# When matching to a parton shower, these generation cuts should be    *
# considerably softer than the analysis cuts.                          *
# (more specific cuts can be specified in SubProcesses/cuts.f)         *
#***********************************************************************
   1  = jetalgo   ! FastJet jet algorithm (1=kT, 0=C/A, -1=anti-kT)
 1.0  = jetradius ! The radius parameter for the jet algorithm
  10  = ptj       ! Min jet transverse momentum
  -1  = etaj      ! Max jet abs(pseudo-rap) (a value .lt.0 means no cut)
#***********************************************************************
# Cuts on the charged leptons (e+, e-, mu+, mu-, tau+ and tau-)        *
# (more specific gen cuts can be specified in SubProcesses/cuts.f)     *
#***********************************************************************
   0  = ptl     ! Min lepton transverse momentum
  -1  = etal    ! Max lepton abs(pseudo-rap) (a value .lt.0 means no cut)
   0  = drll    ! Min distance between opposite sign lepton pairs
   0  = drll_sf ! Min distance between opp. sign same-flavor lepton pairs
   0  = mll     ! Min inv. mass of all opposite sign lepton pairs
  30  = mll_sf  ! Min inv. mass of all opp. sign same-flavor lepton pairs
#***********************************************************************
# Photon-isolation cuts, according to hep-ph/9801442                   *
# When ptgmin=0, all the other parameters are ignored                  *
#***********************************************************************
  20  = ptgmin    ! Min photon transverse momentum
  -1  = etagamma  ! Max photon abs(pseudo-rap)
 0.4  = R0gamma   ! Radius of isolation code
 1.0  = xn        ! n parameter of eq.(3.4) in hep-ph/9801442
 1.0  = epsgamma  ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442
 .true.  = isoEM  ! isolate photons from EM energy (photons and leptons)
#***********************************************************************
# Maximal PDG code for quark to be considered a jet when applying cuts.*
# At least all massless quarks of the model should be included here.   *
#***********************************************************************
 5 = maxjetflavor
#***********************************************************************
EOF

done*.dat <<EOF
EOF
