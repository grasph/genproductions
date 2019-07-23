#!/bin/bash

for MH in 120 125 130; do

d=VBFHToEEG_M${MH}_MLL-0To60_Dalitz_13TeV-amcatnlo

mkdir "$d"

cat > $d/vbfH${MH}_NLO_HtoEEGamma_customizecards.dat <<EOF
set param_card mass 25 ${MH}.
EOF

cat > $d/vbfH${MH}_NLO_HtoEEGamma_extramodels.dat <<EOF
HC_NLO_X0_UFO-v1.3.zip

HC_NLO_X0_UFO-lepton_masses_no_lepton_yukawas-v1.3.tar.gz

EOF

cat > $d/vbfH${MH}_NLO_HtoEEGamma_FKS_params.dat <<EOF
! ==========================================================================
! This file sets the different technical parameters intrinsic to the
! FKS program and which controls the behaviour of the code at run
! time.  The common user should not edit this file and only experts
! should venture editing these parameters.
! ==========================================================================
!
! ==========================================================================
! Arbitrary numerical parameters used in the FKS formalism
! ==========================================================================
!
! To be implemented by the FKS authors
!
! ==========================================================================
! Parameters controlling the tests based on the IR poles comparison
! ==========================================================================
!
! This threshold sets the limiting value for the comparison of the
! relative difference of the IR pole from the OLP and the one computed
! internaly by MadFKS.  The value below is used for the first PS
! points to assess the sanity of the computation. A value ten times
! smaller will be used for the systematic check of IR poles at
! runtime.
! Notice that the systematic comparison of IR poles at run time is
! only performed when the Monte-Carlo over helicity configurations
! method is not used.  Set this value to '-1.0d0' if you want the
! check to always pass.
#IRPoleCheckThreshold
-1d0
! Default :: 1.0d-5
!
! ==========================================================================
! OLP (virtuals) behavior at run time
! ==========================================================================
!
!
! Set the precision required from the OLP code. The IR poles check
! will be performed at run time with a threshold ten times loser than
! the value below. When equal to '-1d0' the default value of the OLP
! is used, and the poles check is disabled at run time
!
#PrecisionVirtualAtRunTime
-1d0
! Default :: 1.0d-3
!
! ==========================================================================
! Parameters defining the techniques used for the MC integration
! ==========================================================================
!
! This integer sets what is the minimum number of contributing
! helicities (in a given subrpocess) which is necessary for MadFKS to
! switch to the Monte-Carlo over helicity configurations method. Set
! this to '-1' if you want to forbid the use of this method
! altogether.
#NHelForMCoverHels
4
! Default :: 4
!
! This parameter sets for which fraction of the events the virtual
! matrix elements should be included. When using MINT, during the
! grid-setup phase, this number will be updated automatically after
! each iteration depending on the relative MC uncertainties.
#VirtualFraction
1.0d0
! Default :: 1.0d0
!
! This parameter sets the minimal fraction of the events for which the
! virtual matrix elements should be included.
#MinVirtualFraction
0.005d0
! Default :: 0.005d0
!
! ==========================================================================
! End of FKS_params.dat file
! ==========================================================================
EOF

cat > $d/vbfH${MH}_NLO_HtoEEGamma_madspin_card.dat <<EOF
import model HC_NLO_X0_UFO-lepton_masses_no_lepton_yukawas --bypass_check

set ms_dir ./madspingrid

set Nevents_for_max_weigth 250 # number of events for the estimate of the max. weight
set max_weight_ps_point 400  # number of PS to estimate the maximum for each event

set max_running_process 1

set spinmode=none
set run_card mmll=0
set run_card mmllmax=60

define l+ = e+
define l- = e-

decay x0 > l+ l- a

launch
EOF

cat > $d/vbfH${MH}_NLO_HtoEEGamma_proc_card.dat <<EOF
import model loop_sm-no_b_mass

generate p p > h j j \$\$ w+ w- z [QCD]

output vbfH${MH}_NLO_HtoEEGamma -nojpeg
EOF

cat > $d/vbfH${MH}_NLO_HtoEEGamma_run_card.dat <<EOF
*******************                                                 
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
   250 = nevents ! Number of unweighted events requested 
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
 0        = ickkw            ! 0 no merging, 3 FxFx merging
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
 0.7  = jetradius ! The radius parameter for the jet algorithm
 0.01  = ptj       ! Min jet transverse momentum
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

done
