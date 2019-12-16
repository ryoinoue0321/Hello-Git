      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                               Katja Fennel   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This routine computes the  biological sources and sinks for the     !
!  Fennel et at. (2006) ecosystem model. Then, it adds those terms     !
!  to the global biological fields.                                    !
!                                                                      !
!  This model is loosely based on the model by Fasham et al. (1990)    !
!  but it differs in many respects.  The detailed equations of the     !
!  nitrogen cycling component  are given in  Fennel et al. (2006).     !
!  Nitrogen is the  fundamental elemental  currency in this model.     !
!  This model was adapted from a code written originally  by  John     !
!  Moisan and Emanule DiLorenzo.                                       !
!                                                                      !
!  It is recommended to activate always the  "BIO_SEDIMENT" option     !
!  to ensure conservation of mass by converting the organic matter     !
!  that is sinking out of the bottom most grid cell into inorganic     !
!  nutrients (i.e.,  instantanaous remineralization  at the water-     !
!  sediment interface). Additionally, the "DENITRIFICATION" option     !
!  can be activated.  Hence, a fraction of the instantenous bottom     !
!  remineralization is  assumed to  occur  through  the  anearobic     !
!  (denitrification)  pathway  and  thus  lost  from the  pool  of     !
!  biologically availalbe fixed nitrogen. See Fennel et al. (2006)     !
!  for details.                                                        !
!                                                                      !
!  Additional  options can be  activated to  enable  simulation of     !
!  inorganic carbon and dissolved oxygen.  Accounting of inorganic     !
!  carbon is activated by the "CARBON" option,  and results in two     !
!  additional  biological  tracer  variables:  DIC and alkalinity.     !
!  See Fennel et al. (2008) for details.                               !
!                                                                      !
!  If the "pCO2_RZ" options is activated, in addition to "CARBON",     !
!  the carbonate system  routines by Zeebe and Wolf-Gladrow (2001)     !
!  are used,  while the  OCMIP  standard routines are the default.     !
!  There are two different ways of treating alkalinity.  It can be     !
!  treated diagnostically (default),  in this case alkalinity acts     !
!  like a passive tracer  that is  not affected  by changes in the     !
!  concentration of  nitrate or ammonium.  However,  if the option     !
!  "TALK_NONCONSERV" is used,  the alkalinity  will be affected by     !
!  sources and sinks in nitrate. See Fennel et al. (2008) for more     !
!  details.                                                            !
!                                                                      !
!  If the "OXYGEN" option is activated,  one additional biological     !
!  tracer variable for dissolved oxygen. "OXYGEN" can be activated     !
!  independently of the  "CARBON"  option. If "OCMIP_OXYGEN_SC" is     !
!  used, in addition to "OXYGEN",  the Schmidt number of oxygen in     !
!  seawater will be  computed  using the  formulation  proposed by     !
!  Keeling et al. (1998, Global Biogeochem. Cycles,  12, 141-163).     !
!  Otherwise, the Wanninkhof's (1992) formula will be used.            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fennel, K., Wilkin, J., Levin, J., Moisan, J., O'Reilly, J.,      !
!      Haidvogel, D., 2006: Nitrogen cycling in the Mid Atlantic       !
!      Bight and implications for the North Atlantic nitrogen          !
!      budget: Results from a three-dimensional model.  Global         !
!      Biogeochemical Cycles 20, GB3007, doi:10.1029/2005GB002456.     !
!                                                                      !
!    Fennel, K., Wilkin, J., Previdi, M., Najjar, R. 2008:             !
!      Denitrification effects on air-sea CO2 flux in the coastal      !
!      ocean: Simulations for the Northwest North Atlantic.            !
!      Geophys. Res. Letters 35, L24608, doi:10.1029/2008GL036147.     !
!                                                                      !
!***********************************************************************
!
!
! "fennel_2019.h" was edited by Inoue in 2019.
!ver3は植物プランクトンの栄養塩依存項を陽的に示していたが，
!ver4ではそれを陰的になるよう修正している
!vev3だとそのあたりを陽的に解きすぎており，
!数値誤差が大きくなる可能性がある
!一方でver4だと，PHYT1の増殖傾向が一番大きくなりそうな
!表現ではあるが，数値誤差は小さくなる可能性がある
!ver5では，データ同化におけるTLMコード，ADMコード作成がやりやすくなるような改良を加えている
!ver6では，主にtl_fennel_2_param.hが組みやすいような細かい修正が加わっている
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
#ifdef ADJUST_PARAM
     &                   Nparam(ng),                                    &
#endif
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#if defined CARBON || defined OXYGEN
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif
#endif
#ifdef ADJUST_PARAM
     &                   OCEAN(ng) % p,                                 &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif

      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
#ifdef ADJUST_PARAM
     &                         UBp,                                     &
#endif
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#if defined CARBON || defined OXYGEN
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, svstr,                            &
# endif
#endif
#ifdef ADJUST_PARAM
     &                         p,                                       &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel  ! (okada)
      USE mod_iounits   ! (okada)
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
#ifdef ADJUST_PARAM
      integer, intent(in) :: UBp
#endif
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# if defined CARBON || defined OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
# endif
# ifdef ADJUST_PARAM
      real(r8), intent(in) :: p(:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)

#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# if defined CARBON || defined OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef ADJUST_PARAM
      real(r8), intent(in) :: p(2,UBp)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
#ifdef PHOSPHORUS
# if defined PHYT2 || defined PHYT3
      integer, parameter :: Nsink = 8
# elif defined PHYT3
      integer, parameter :: Nsink = 10
# else
      integer, parameter :: Nsink = 6
# endif
#else
      integer, parameter :: Nsink = 4
#endif

      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#if defined CARBON || defined OXYGEN
      real(r8) :: u10squ
#endif
#ifdef OXYGEN
      real(r8), parameter :: OA0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: OA1 = 3.22014_r8       ! saturation
      real(r8), parameter :: OA2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: OA3 = 4.94457_r8
      real(r8), parameter :: OA4 =-0.256847_r8
      real(r8), parameter :: OA5 = 3.88767_r8
      real(r8), parameter :: OB0 =-0.00624523_r8
      real(r8), parameter :: OB1 =-0.00737614_r8
      real(r8), parameter :: OB2 =-0.0103410_r8
      real(r8), parameter :: OB3 =-0.00817083_r8
      real(r8), parameter :: OC0 =-0.000000488682_r8
      real(r8), parameter :: rOxNO3= 8.625_r8       ! 138/16
      real(r8), parameter :: rOxNH4= 6.625_r8       ! 106/16
      real(r8) :: l2mol = 1000.0_r8/22.3916_r8      ! liter to mol (mol/l*l/m3)
      real(r8) :: molv = 22.3916_r8                 ! l/mol
      real(r8) :: rho_O2 = 1.42903_r8               ! g/l
      real(r8) :: mol2g_O2 = 1.42903_r8*22.3916_r8  ! g/mol
      real(r8), dimension(115,122)  :: FSOD_distribution  !inoue, 2018
      integer :: FI,FJ
#endif
      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: fac1, fac2, fac3, fac4, fac5, fac6
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: fac

      real(r8) :: total_N

      real(r8) :: Epp1, L_NH4_1, L_NO3_1, LTOT1, Vp1
      real(r8) :: Chl2C1, t_PPmax1, inhNH4_1
      real(r8) :: Att_T1, Att_L1, Att_N1

      real(r8) :: cff01, cff11, cff21, cff31, cff41, cff51, cff61
      real(r8) :: fac11, fac21, fac31, fac41, fac51, fac61
#if defined PHYT2 || defined PHYT3
      real(r8) :: Epp2, L_NH4_2, L_NO3_2, LTOT2, Vp2
      real(r8) :: Chl2C2, t_PPmax2, inhNH4_2
      real(r8) :: Att_T2, Att_L2, Att_N2

      real(r8) :: cff02, cff12, cff22, cff32, cff42, cff52, cff62
      real(r8) :: fac12, fac22, fac32, fac42, fac52, fac62
#endif
#ifdef PHYT3
      real(r8) :: Epp3, L_NH4_3, L_NO3_3, LTOT3, Vp3
      real(r8) :: Chl2C3, t_PPmax3, inhNH4_3
      real(r8) :: Att_T3, Att_L3, Att_N3

      real(r8) :: cff03, cff13, cff23, cff33, cff43, cff53, cff63
      real(r8) :: fac13, fac23, fac33, fac43, fac53, fac63
#endif

#ifdef OXYGEN
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA
#endif

      real(r8) :: N_Flux_Assim
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_Egest
      real(r8) :: N_Flux_NewProd, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_Remine
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo
      real(r8) :: N_Flux_Presp

      real(r8) :: N_Flux_Assim1
      real(r8) :: N_Flux_CoagD1, N_Flux_CoagP1
      real(r8) :: N_Flux_Egest1
      real(r8) :: N_Flux_NewProd1, N_Flux_RegProd1
      real(r8) :: N_Flux_Nitrifi1
      real(r8) :: N_Flux_Pmortal1, N_Flux_Zmortal1
      real(r8) :: N_Flux_Remine1
      real(r8) :: N_Flux_Zexcret1, N_Flux_Zmetabo1
      real(r8) :: N_Flux_Presp1
#if defined PHYT2 || defined PHYT3
      real(r8) :: N_Flux_Assim2
      real(r8) :: N_Flux_CoagD2, N_Flux_CoagP2
      real(r8) :: N_Flux_Egest2
      real(r8) :: N_Flux_NewProd2, N_Flux_RegProd2
      real(r8) :: N_Flux_Nitrifi2
      real(r8) :: N_Flux_Pmortal2, N_Flux_Zmortal2
      real(r8) :: N_Flux_Remine2
      real(r8) :: N_Flux_Zexcret2, N_Flux_Zmetabo2
      real(r8) :: N_Flux_Presp2
#endif
#ifdef PHYT3
      real(r8) :: N_Flux_Assim3
      real(r8) :: N_Flux_CoagD3, N_Flux_CoagP3
      real(r8) :: N_Flux_Egest3
      real(r8) :: N_Flux_NewProd3, N_Flux_RegProd3
      real(r8) :: N_Flux_Nitrifi3
      real(r8) :: N_Flux_Pmortal3, N_Flux_Zmortal3
      real(r8) :: N_Flux_Remine3
      real(r8) :: N_Flux_Zexcret3, N_Flux_Zmetabo3
      real(r8) :: N_Flux_Presp3
#endif

      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

#ifdef PHOSPHORUS
      real(r8) :: L_PO4, LMIN
      real(r8) :: P_Flux

      real(r8) :: L_PO4_1, LMIN1
      real(r8) :: P_Flux1
# if defined PHYT2 || defined PHYT3
      real(r8) :: L_PO4_2, LMIN2
      real(r8) :: P_Flux2
# endif
# ifdef PHYT3
      real(r8) :: L_PO4_3, LMIN3
      real(r8) :: P_Flux3
# endif

      real(r8), parameter :: rOxPO4 = 106.0_r8   ! 106/1
#endif
#ifdef H2S
      real(r8) :: S_Flux

      real(r8), parameter :: rOxH2S = 2.0_r8     !?
#endif
#define BLOWINGUP_CHECKER
#ifdef BLOWINGUP_CHECKER
      integer :: ii
      real(r8), dimension(NT(ng)) :: val
      character (len=8) :: valchar
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
#ifdef DIAGNOSTICS_BIO
!
!  A factor to account for the number of iterations in accumulating
!  diagnostic rate variables.
!
      fiter=1.0_r8/REAL(BioIter(ng),r8)
#endif
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt1
      idsink(2)=iChlo1
      idsink(3)=iSDeN
      idsink(4)=iLDeN
#ifdef PHOSPHORUS
      idsink(5)=iSDeP
      idsink(6)=iLDeP
#endif
#if defined PHYT2 || defined PHYT3
      idsink(7)=iChlo2
      idsink(8)=iPhyt2
#endif
#ifdef PHYT3
      idsink(9)=iChlo3
      idsink(10)=iPhyt3
#endif
#ifdef ADJUST_PARAM
!
!  Convert parameters.
!
# ifdef EXP_PARAM
      wPhy(ng)=wPhy(ng)*EXP(p(nstp,iwPhy))
      wSDet(ng)=wSDet(ng)*EXP(p(nstp,iwSDet))
      wLDet(ng)=wLDet(ng)*EXP(p(nstp,iwLDet))
      Chl2C_m1(ng)=Chl2C_m1(ng)*EXP(p(nstp,iChl2C_m1))
      Vp0(ng)=Vp0(ng)*EXP(p(nstp,iVp0))
      R_SODf(ng)=R_SODf(ng)*EXP(p(nstp,iR_SODf))
      R_NH4f(ng)=R_NH4f(ng)*EXP(p(nstp,iR_NH4f))
      R_PO4f(ng)=R_PO4f(ng)*EXP(p(nstp,iR_PO4f))
# else
      AttSW(ng)=p(nstp,iAttSW)
      AttChl(ng)=p(nstp,iAttChl)
      Vp0(ng)=p(nstp,iVp0)
      I_thNH4(ng)=p(nstp,iI_thNH4)
      D_p5NH4(ng)=p(nstp,iD_p5NH4)
      K_NO3_1(ng)=p(nstp,iK_NO3_1)
      K_NH4_1(ng)=p(nstp,iK_NH4_1)
      K_Phy1(ng)=p(nstp,iK_Phy1)
      Chl2C_m1(ng)=p(nstp,iChl2C_m1)
      g_max1(ng)=p(nstp,ig_max1)

      PhyCN(ng)=p(nstp,iPhyCN)
      PhyIP(ng)=p(nstp,iPhyIP)
      PhyIS1(ng)=p(nstp,iPhyIS1)
      PhyMR(ng)=p(nstp,iPhyMR)
      PhyPR(ng)=p(nstp,iPhyPR)	! Noda, 2018
      PhyBR(ng)=p(nstp,iPhyBR)	! Noda, 2018
      ZooAE_N(ng)=p(nstp,iZooAE_N)
      ZooBM(ng)=p(nstp,iZooBM)
      ZooCN(ng)=p(nstp,iZooCN)
      ZooER(ng)=p(nstp,iZooER)
      ZooGR(ng)=p(nstp,iZooGR)
      ZooMR(ng)=p(nstp,iZooMR)
      LDeRRN(ng)=p(nstp,iLDeRRN)
      CoagR(ng)=p(nstp,iCoagR)
      SDeRRN(ng)=p(nstp,iSDeRRN)
      wPhy1(ng)=p(nstp,iwPhy1)
      wSDet(ng)=p(nstp,iwSDet)
      wLDet(ng)=p(nstp,iwLDet)

      K_Nitri(ng)=p(nstp,iK_Nitri)
      NitriR(ng)=p(nstp,iNitriR)
      K_Denit(ng)=p(nstp,iK_Denit)
      DenitR(ng)=p(nstp,iDenitR)
      K_PO4_1(ng)=p(nstp,iK_PO4_1)
      PhyPN(ng)=p(nstp,iPhyPN)
      ZooPN(ng)=p(nstp,iZooPN)
      K_DO(ng)=p(nstp,iK_DO)
      LDeRRP(ng)=p(nstp,iLDeRRP)
      SDeRRP(ng)=p(nstp,iSDeRRP)
      R_SODf(ng)=p(nstp,iR_SODf)
      R_NH4f(ng)=p(nstp,iR_NH4f)
      R_PO4f(ng)=p(nstp,iR_PO4f)
      R_NH4f_max(ng)=p(nstp,iR_NH4f_m)
      R_PO4f_max(ng)=p(nstp,iR_PO4f_m)
      K_DO_npflux(ng)=p(nstp,iK_DO_npf)
      t_SODf(ng)=p(nstp,it_SODf)

! added parameters by Inoue(2019)
#  if defined OXYGEN_SPATIAL_DIST1 || defined OXYGEN_SPATIAL_DIST2
      R_SODfa(ng)=p(nstp,iR_SODfa)
      R_SODfb(ng)=p(nstp,iR_SODfb)
      R_SODfc(ng)=p(nstp,iR_SODfc)
      R_SODfd(ng)=p(nstp,iR_SODfd)
      Alpha_SODf(ng)=p(nstp,iAlpha_SODf)
#  endif
#  if defined PHYT2 || defined PHYT3
      g_max2(ng)=p(nstp,ig_max2)
      PhyIS2(ng)=p(nstp,iPhyIS2)
      K_NO3_2(ng)=p(nstp,iK_NO3_2)
      K_NH4_2(ng)=p(nstp,iK_NH4_2)
      K_NO3_2(ng)=p(nstp,iK_NO3_2)
      K_Phy2(ng)=p(nstp,iK_Phy2)
      Chl2C_m2(ng)=p(nstp,iChl2C_m2)
      wPhy2(ng)=p(nstp,iwPhy2)
#   ifdef PHYT3
      g_max3(ng)=p(nstp,ig_max3)
      PhyIS3(ng)=p(nstp,iPhyIS3)
      K_NO3_3(ng)=p(nstp,iK_NO3_3)
      K_NH4_3(ng)=p(nstp,iK_NH4_3)
      K_NO3_3(ng)=p(nstp,iK_NO3_3)
      K_Phy3(ng)=p(nstp,iK_Phy3)
      Chl2C_m3(ng)=p(nstp,iChl2C_m3)
      wPhy3(ng)=p(nstp,iwPhy3)
#   endif
#  endif
#  if defined GROWTH_TOPT1 || defined GROWTH_TOPT2 ||defined GROWTH_TOPT3
      t_opt1(ng)=p(nstp,it_opt1)
#   if defined PHYT2 || defined PHYT3
      t_opt2(ng)=p(nstp,it_opt2)
#    ifdef PHYT3
      t_opt3(ng)=p(nstp,it_opt3)
#    endif
#   endif
#  endif
#  if defined GROWTH_TOPT2 ||defined GROWTH_TOPT3
      t_max1(ng)=p(nstp,it_max1)
#   if defined PHYT2 || defined PHYT3
      t_max2(ng)=p(nstp,it_max2)
#    ifdef PHYT3
      t_max3(ng)=p(nstp,it_max3)
#    endif
#   endif
#  endif
#  ifdef GROWTH_TOPT2 
      beta01(ng)=p(nstp,ibeta01)
#   if defined PHYT2 || defined PHYT3
      beta02(ng)=p(nstp,ibeta02)
#    ifdef PHYT3
      beta03(ng)=p(nstp,ibeta03)
#    endif
#   endif
#  endif
#  ifdef GROWTHL2
      I_opt1(ng)=p(nstp,iI_opt1)
#   if defined PHYT2 || defined PHYT3
      I_opt2(ng)=p(nstp,iI_opt2)
#    ifdef PHYT3
      I_opt3(ng)=p(nstp,iI_opt3)
#    endif
#   endif
#  endif

# endif
#endif
#ifdef CHECKER
!
! stdout 
!
      if (master.and.(mod(iic(ng)-1,ninfo(ng)).eq.0)) then
# ifdef ADJUST_PARAM
        write(stdout,101) 'NL 01:08', p(nstp,1:8)
        write(stdout,101) 'NL 09:16', p(nstp,9:16)
        write(stdout,101) 'NL 17:24', p(nstp,17:24)
        write(stdout,101) 'NL 25:32', p(nstp,25:32)
        write(stdout,101) 'NL 33:40', p(nstp,33:40)
        write(stdout,101) 'NL 41:  ', p(nstp,41:)
        write(stdout,*) ('-',i=1,78)
# endif
      end if
 101  FORMAT (a,8(1pe9.1))
#endif
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy1(ng)                ! phytoplankton
      Wbio(2)=wPhy1(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus
#ifdef PHOSPHORUS
      Wbio(5)=wSDet(ng)               ! small Carbon- or Phosphorus-detritus
      Wbio(6)=wLDet(ng)               ! large Carbon- or Phosphorus-detritus
#endif
#if defined PHYT2 || defined PHYT3
      Wbio(7)=wPhy2(ng)                ! chlorophyll2
      Wbio(8)=wPhy2(ng)                ! phytoplankton2
#endif
#ifdef PHYT3
      Wbio(9)=wPhy3(ng)                ! chlorophyll3
      Wbio(10)=wPhy3(ng)               ! phytoplankton3
#endif
#ifdef BIO_SED_DIAGENESIS
      dia_count=dia_count+1
#endif
# ifdef OXYGEN_SPATIAL_DIST1
      OPEN(2117,FILE='FSOD_dist1.csv',STATUS='old')
      REWIND(2117)
      DO FI=1,115
	  READ(2117,*) (FSOD_distribution(FI,FJ),FJ=1,122)
      END DO
      CLOSE(2117)

# elif defined OXYGEN_SPATIAL_DIST2
      OPEN(2117,FILE='FSOD_dist2.csv',STATUS='old')
      REWIND(2117)
      DO FI=1,115
	  READ(2117,*) (FSOD_distribution(FI,FJ),FJ=1,122)
      END DO
      CLOSE(2117)
# endif
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
#ifdef BLOWINGUP_CHECKER
          val(ibio)=0.0_r8
#endif
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            AttFac=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  1. Compute average light attenuation for each grid cell. 
!     To include other attenuation contributions like suspended sediment 
!     or CDOM modify AttFac.
!

                fac=                           &
#if defined PHYT2 || defined PHYT3
     &              Bio(i,k,iChlo2)+           &
# ifdef PHYT3
     &              Bio(i,k,iChlo3)+           &
# endif
#endif
     &              Bio(i,k,iChlo1)
!
                fac1=AttSW(ng)+AttChl(ng)*fac+AttFac
                Att=fac1*(z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                fac2=Itop*(1.0_r8-ExpAtt)    
                PAR=fac2/Att    ! average at cell center
!
!  2. Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8
                fac11=Bio(i,k,iPhyt1)*cff+eps
                fac21=Bio(i,k,iChlo1)/fac11
                Chl2C1=MIN(fac21,Chl2C_m1(ng))
# if defined PHYT2 || defined PHYT3
!
                fac12=Bio(i,k,iPhyt2)*cff+eps
                fac22=Bio(i,k,iChlo2)/fac12
                Chl2C2=MIN(fac22,Chl2C_m2(ng))
# endif
# ifdef PHYT3
!
                fac13=Bio(i,k,iPhyt3)*cff+eps
                fac23=Bio(i,k,iChlo3)/fac13
                Chl2C3=MIN(fac23,Chl2C_m3(ng))
# endif

!
!  3. Compute attenuations for phytoplankton growth.
!
#ifdef GROWTH_TOPT1
!
!  Temperature attenuation with t_opt, beta1 and beta2.
!
                cff01=(Bio(i,k,itemp)-t_opt1(ng))**2.0_r8
                IF (Bio(i,k,itemp).le.t_opt1(ng)) THEN
                  Att_T1=EXP(-beta1(ng)*cff01)
                ELSE
                  Att_T1=EXP(-beta2(ng)*cff01)
                END IF
# if defined PHYT2 || defined PHYT3
!
                cff02=(Bio(i,k,itemp)-t_opt2(ng))**2.0_r8
                IF (Bio(i,k,itemp).le.t_opt2(ng)) THEN
                  Att_T2=EXP(-beta1(ng)*cff02)
                ELSE
                  Att_T2=EXP(-beta2(ng)*cff02)
                END IF
# endif
# ifdef PHYT3
!
                cff03=(Bio(i,k,itemp)-t_opt3(ng))**2.0_r8
                IF (Bio(i,k,itemp).le.t_opt3(ng)) THEN
                  Att_T3=EXP(-beta1(ng)*cff03)
                ELSE
                  Att_T3=EXP(-beta2(ng)*cff03)
                END IF
# endif
#elif defined GROWTH_TNORMAL
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                 Att_T1=1.066_r8**(Bio(i,k,itemp)-20.0_r8)
# if defined PHYT2 || defined PHYT3
!
                 Att_T2=1.066_r8**(Bio(i,k,itemp)-20.0_r8)
# endif
# ifdef PHYT3
!
                 Att_T3=1.066_r8**(Bio(i,k,itemp)-20.0_r8)
# endif
#elif defined GROWTH_TOPT2
!
!  Temperature attenuation with t_opt, t_max and beta. 
!  [applied to MICRO LOOP for Ise Bay]
!
                fac11=beta01(ng)*(Bio(i,k,itemp)-t_opt1(ng))
                fac21=EXP(fac11)
                IF (Bio(i,k,itemp).le.t_max1(ng)) THEN
                  fac31=t_max1(ng)-Bio(i,k,itemp)
                  cff11=t_max1(ng)-t_opt1(ng)
                  cff21=beta01(ng)*cff11
                  fac41=(fac31/cff11)**cff21
                  Att_T1=fac21*fac41
                ELSE
                  Att_T1=0.0_r8
                END IF
# if defined PHYT2 || defined PHYT3
!
                fac12=beta02(ng)*(Bio(i,k,itemp)-t_opt2(ng))
                fac22=EXP(fac12)
                IF (Bio(i,k,itemp).le.t_max2(ng)) THEN
                  fac32=t_max2(ng)-Bio(i,k,itemp)
                  cff12=t_max2(ng)-t_opt2(ng)
                  cff22=beta02(ng)*cff12
                  fac42=(fac32/cff12)**cff22
                  Att_T2=fac22*fac42
                ELSE
                  Att_T2=0.0_r8
                END IF
# endif
# ifdef PHYT3
!
                fac13=beta03(ng)*(Bio(i,k,itemp)-t_opt3(ng))
                fac23=EXP(fac13)
                IF (Bio(i,k,itemp).le.t_max3(ng)) THEN
                  fac33=t_max3(ng)-Bio(i,k,itemp)
                  cff13=t_max3(ng)-t_opt3(ng)
                  cff23=beta03(ng)*cff13
                  fac43=(fac33/cff13)**cff23
                  Att_T3=fac23*fac43
                ELSE
                  Att_T3=0.0_r8
                END IF
# endif
#elif defined GROWTH_TOPT3
!
!  Temperature attenuation with t_opt and t_max. 
!  [applied to BOX model for Tokyo Bay]
!
                IF (Bio(i,k,itemp).lt.t_opt1(ng)) THEN
                  fac11=Bio(i,k,itemp)-t_opt1(ng)*2.0_r8
                  Att_T1=Bio(i,k,itemp)/(t_opt1(ng)**2.0_r8)*fac11
                ELSE
	           fac11=Bio(i,k,itemp)-t_opt1(ng)
                  fac21=t_max1(ng)+t_opt1(ng)
                  Att_T1=1.0_r8-(fac11/fac21)**2.0_r8
                END IF
# if defined PHYT2 || defined PHYT3
!
                IF (Bio(i,k,itemp).lt.t_opt2(ng)) THEN
                  fac12=Bio(i,k,itemp)-t_opt2(ng)*2.0_r8
                  Att_T2=Bio(i,k,itemp)/(t_opt2(ng)**2.0_r8)*fac12
                ELSE
	           fac12=Bio(i,k,itemp)-t_opt2(ng)
                  fac22=t_max2(ng)+t_opt2(ng)
                  Att_T2=1.0_r8-(fac12/fac22)**2.0_r8
                END IF
# endif
# ifdef PHYT3
!
                IF (Bio(i,k,itemp).lt.t_opt3(ng)) THEN
                  fac13=Bio(i,k,itemp)-t_opt3(ng)*2.0_r8
                  Att_T3=Bio(i,k,itemp)/(t_opt3(ng)**2.0_r8)*fac13
                ELSE
	           fac13=Bio(i,k,itemp)-t_opt3(ng)
                  fac23=t_max3(ng)+t_opt3(ng)
                  Att_T3=1.0_r8-(fac13/fac23)**2.0_r8
                END IF
# endif
#endif
#ifdef GROWTH_L1
                fac11=PAR*PhyIS1(ng)
                fac21=SQRT(g_max1(ng)*g_max1(ng)+fac11*fac11)
                Att_L1=fac11/fac21
# if defined PHYT2 || defined PHYT3
!
                fac12=PAR*PhyIS2(ng)
                fac22=SQRT(g_max2(ng)*g_max2(ng)+fac12*fac12)
                Att_L2=fac12/fac22
# endif
# ifdef PHYT3
!
                fac13=PAR*PhyIS3(ng)
                fac23=SQRT(g_max3(ng)*g_max3(ng)+fac13*fac13)
                Att_L3=fac13/fac23
# endif
#elif defined GROWTH_L2
!
!  Light attenuation with I_opt.
!  [applied to RCA model]
!
                Att=Att/2.0_r8
                cff01=Itop/I_opt1(ng)
                Att_L1=EXP(1.0_r8)/Att*(EXP(-cff01*Exp(-Att))-EXP(-cff01))
# if defined PHYT2 || defined PHYT3
!
                cff02=Itop/I_opt2(ng)
                Att_L2=EXP(1.0_r8)/Att*(EXP(-cff02*Exp(-Att))-EXP(-cff02))
# endif
# ifdef PHYT3
!
                cff03=Itop/I_opt3(ng)
                Att_L3=EXP(1.0_r8)/Att*(EXP(-cff03*Exp(-Att))-EXP(-cff03))
# endif
#elif defined GROWTH_L3
!
!  Photo-inhabitation of phytoplankton photosynthesis
!  (c) 2016-7-4 Teruhisa Okada
!
                cff01=PAR/I_opt1(ng)
                Att_L1=cff01*exp(1.0_r8-cff01)
# if defined PHYT2 || defined PHYT3
!
                cff02=PAR/I_opt2(ng)
                Att_L2=cff02*exp(1.0_r8-cff02)
# endif
# ifdef PHYT3
!
                cff03=PAR/I_opt3(ng)
                Att_L3=cff03*exp(1.0_r8-cff03)
# endif
#endif
!
                t_PPmax1=g_max1(ng)*Att_T1*Att_L1  
#if defined PHYT2 || defined PHYT3
                t_PPmax2=g_max2(ng)*Att_T2*Att_L2
#endif
#ifdef PHYT3
                t_PPmax3=g_max3(ng)*Att_T3*Att_L3
#endif
!
#ifdef GROWTH_N1
!
!  Nutrient-limitation terms (Laurent et al. 2012).
!
                cff11=Bio(i,k,iNH4_)*K_NH4_1(ng)
                cff21=Bio(i,k,iNO3_)*K_NO3_1(ng)
                inhNH4_1=1.0_r8/(1.0_r8+cff11)
                L_NH4_1=cff11/(1.0_r8+cff11)
                L_NO3_1=cff21*inhNH4_1/(1.0_r8+cff21)
                LTOT1=L_NO3_1+L_NH4_1
#  ifdef PHOSPHORUS
!
                cff31=Bio(i,k,iPO4_)*K_PO4_1(ng)
                L_PO4_1=cff31/(1.0_r8+cff31)
                Att_N1=MIN(LTOT1,L_PO4_1)
#  else
!
                Att_N1=LTOT1
#  endif
# if defined PHYT2 || defined PHYT3
                cff12=Bio(i,k,iNH4_)*K_NH4_2(ng)
                cff22=Bio(i,k,iNO3_)*K_NO3_2(ng)
                inhNH4_2=1.0_r8/(1.0_r8+cff12)
                L_NH4_2=cff12/(1.0_r8+cff12)
                L_NO3_2=cff22*inhNH4_2/(1.0_r8+cff22)
                LTOT2=L_NO3_2+L_NH4_2
#  ifdef PHOSPHORUS
!
                cff32=Bio(i,k,iPO4_)*K_PO4_2(ng)
                L_PO4_2=cff32/(1.0_r8+cff32)
                Att_N2=MIN(LTOT2,L_PO4_2)
#  else
!
                Att_N2=LTOT2
#  endif
# endif
# ifdef PHYT3
                cff13=Bio(i,k,iNH4_)*K_NH4_3(ng)
                cff23=Bio(i,k,iNO3_)*K_NO3_3(ng)
                inhNH4_3=1.0_r8/(1.0_r8+cff13)
                L_NH4_3=cff13/(1.0_r8+cff13)
                L_NO3_3=cff23*inhNH4_3/(1.0_r8+cff23)
                LTOT3=L_NO3_3+L_NH4_3
#  ifdef PHOSPHORUS
!
                cff33=Bio(i,k,iPO4_)*K_PO4_3(ng)
                L_PO4_3=cff33/(1.0_r8+cff33)
                Att_N3=MIN(LTOT3,L_PO4_3)
#  else
!
                Att_N3=LTOT3
#  endif
# endif
#elif GROWTH_N2
!
! Arranged attenuation (Juan.A et al. 2011)
!
#endif
!
!  4. Compute nutrients uptake by Phytoplankton.
!
                fac11=dtdays*t_PPmax1*Bio(i,k,iPhyt1)
                cff41=fac11*K_NO3_1(ng)*inhNH4_1/(1.0_r8+cff21)
                cff51=fac11*K_NH4_1(ng)/(1.0_r8+cff11)
# if defined PHYT2 || defined PHYT3
!
                fac12=dtdays*t_PPmax2*Bio(i,k,iPhyt2)
                cff42=fac12*K_NO3_2(ng)*inhNH4_2/(1.0_r8+cff22)
                cff52=fac12*K_NH4_2(ng)/(1.0_r8+cff12)
# endif
# ifdef PHYT3
!
                fac13=dtdays*t_PPmax3*Bio(i,k,iPhyt3)
                cff43=fac13*K_NO3_3(ng)*inhNH4_3/(1.0_r8+cff23)
                cff53=fac13*K_NH4_3(ng)/(1.0_r8+cff13)
# endif
!
# ifdef PHOSPHORUS
                IF (Att_N1.eq.L_PO4_1) THEN
                  cff61=fac11*PhyPN(ng)*K_PO4_1(ng)/(1.0_r8+cff31)
                  cff41=0.0_r8
                  cff51=0.0_r8
                ELSE
                  cff61=0.0_r8
                END IF
#  if defined PHYT2 || defined PHYT3
!
                IF (Att_N2.eq.L_PO4_2) THEN
                  cff62=fac12*PhyPN(ng)*K_PO4_2(ng)/(1.0_r8+cff32)
                  cff42=0.0_r8
                  cff52=0.0_r8
                ELSE
                  cff62=0.0_r8
                END IF          
#  endif
#  ifdef PHYT3
!
                IF (Att_N3.eq.L_PO4_3) THEN
                  cff63=fac13*PhyPN(ng)*K_PO4_3(ng)/(1.0_r8+cff33)
                  cff43=0.0_r8
                  cff53=0.0_r8
                ELSE
                  cff63=0.0_r8
                END IF
#  endif
# endif
!
!  4-1. Compute phosphorus-preponderate uptake by Phytoplankton.
!
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff61)        
                P_Flux1=Bio(i,k,iPO4_)*cff61
                N_Flux_NewProd1=P_Flux1/PhyPN(ng)*L_NO3_1/MAX(LTOT1,eps)
                N_Flux_RegProd1=P_Flux1/PhyPN(ng)*L_NH4_1/MAX(LTOT1,eps)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd1
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd1
# if defined PHYT2 || defined PHYT3
!
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff62)
                P_Flux2=Bio(i,k,iPO4_)*cff62
                N_Flux_NewProd2=P_Flux2/PhyPN(ng)*L_NO3_2/MAX(LTOT2,eps)
                N_Flux_RegProd2=P_Flux2/PhyPN(ng)*L_NH4_2/MAX(LTOT2,eps)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd2
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd2
# endif
# ifdef PHYT3
!
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff63)
                P_Flux3=Bio(i,k,iPO4_)*cff63
                N_Flux_NewProd3=P_Flux3/PhyPN(ng)*L_NO3_3/MAX(LTOT3,eps)
                N_Flux_RegProd3=P_Flux3/PhyPN(ng)*L_NH4_3/MAX(LTOT3,eps)
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd3
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd3
# endif
!
!  4-2. Compute nitrate and ammonium-preponderate uptake by Phytoplankton.
!
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff41)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff51)
                N_Flux_NewProd1=Bio(i,k,iNO3_)*cff41
                N_Flux_RegProd1=Bio(i,k,iNH4_)*cff51
                P_Flux1=(N_Flux_NewProd1+N_Flux_RegProd1)*PhyPN(ng) 
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux1
# if defined PHYT2 || defined PHYT3
!
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff42)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff52)
                N_Flux_NewProd2=Bio(i,k,iNO3_)*cff42
                N_Flux_RegProd2=Bio(i,k,iNH4_)*cff52
                P_Flux2=(N_Flux_NewProd2+N_Flux_RegProd2)*PhyPN(ng)
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux2
# endif
# ifdef PHYT3
!
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff43)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff53)
                N_Flux_NewProd3=Bio(i,k,iNO3_)*cff43
                N_Flux_RegProd3=Bio(i,k,iNH4_)*cff53
                P_Flux3=(N_Flux_NewProd3+N_Flux_RegProd3)*PhyPN(ng)
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux3
# endif

!
!  5. Compute source/Sink about Phyt, Chlo and Oxyg.
!
                Bio(i,k,iPhyt1)=Bio(i,k,iPhyt1)+                          &
     &                         N_Flux_NewProd1+N_Flux_RegProd1
                Bio(i,k,iChlo1)=Bio(i,k,iChlo1)+                          &
     &                         (dtdays*t_PPmax1*t_PPmax1*Att_N1*Att_N1*   &
     &                          Chl2C_m1(ng)*Bio(i,k,iChlo1))/            &
     &                         (PhyIS1(ng)*MAX(Chl2C1,eps)*PAR+eps)
# if defined PHYT2 || defined PHYT3
!
                Bio(i,k,iPhyt2)=Bio(i,k,iPhyt2)+                          &
     &                         N_Flux_NewProd2+N_Flux_RegProd2
                Bio(i,k,iChlo2)=Bio(i,k,iChlo2)+                          &
     &                         (dtdays*t_PPmax2*t_PPmax2*Att_N2*Att_N2*   &
     &                          Chl2C_m2(ng)*Bio(i,k,iChlo2))/            &
     &                         (PhyIS2(ng)*MAX(Chl2C2,eps)*PAR+eps)
# endif
# ifdef PHYT3
!
                Bio(i,k,iPhyt3)=Bio(i,k,iPhyt3)+                          &
     &                         N_Flux_NewProd3+N_Flux_RegProd3
                Bio(i,k,iChlo3)=Bio(i,k,iChlo3)+                          &
     &                         (dtdays*t_PPmax3*t_PPmax3*Att_N3*Att_N3*   &
     &                          Chl2C_m3(ng)*Bio(i,k,iChlo3))/            &
     &                         (PhyIS3(ng)*MAX(Chl2C3,eps)*PAR+eps)
#endif
#ifdef OXYGEN
!
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)                             &
# if defined PHYT2 || defined PHYT3
     &                                       +N_Flux_NewProd2*rOxNO3      &
     &                                       +N_Flux_RegProd2*rOxNH4      &
# endif
# ifdef PHYT3
     &                                       +N_Flux_NewProd3*rOxNO3      &
     &                                       +N_Flux_RegProd3*rOxNH4      &
# endif
     &                                       +N_Flux_NewProd1*rOxNO3      &
     &                                       +N_Flux_RegProd1*rOxNH4      
#endif
!
!  6. Compute the Nitrification of NH4 ==> NO3 is thought to occur only in dark 
!     and only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
                fac=dtdays*NitriR(ng)
#ifdef OXYGEN
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                fac1=fac*fac3*cff
# else
                fac1=fac*fac3
# endif
#else
                fac1=fac
#endif
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
!
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
#endif

!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE
              fac1=dtdays*NitriR(ng)
              DO k=N(ng),1,-1
#if defined OXYGEN
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)

# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                cff3=fac1*fac3*cff

# else
                cff3=fac1*fac3
# endif
#else
                cff3=fac1
#endif
!
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
#endif
              END DO
            END IF
          END DO
#if defined OXYGEN && defined DENITRIFICATION
!
!  7. Compute denitrification in anoxic water.   (Okada 2014/02/13)
!
          DO i=Istr,Iend
            DO k=1,N(ng)
              fac1=dtdays*DenitR(ng)
              cff1=MAX(Bio(i,k,iOxyg),0.0_r8)/K_Denit(ng)
              cff2=1.0_r8/(1.0_r8+cff1)

# ifdef TDEPENDANCE
              fac2=DenitR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              cff3=cff2*fac1*fac2
# else
              cff3=cff2*fac1
# endif
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff3)
            END DO
          END DO
#endif

! x =================================================================== 
!    Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!    assimilated to zooplankton (fraction: ZooAE_N) and phytoplankton
!    respiration(rate: PhyBR and PhyPR), and phytoplankton mortality 
!    (rate: PhyMR) to small detritus. [Landry 1993 L&O 38:468-472]
! x =================================================================== 

          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  8. Compute phytoplankton grazing by zooplankton.
!
#ifdef TDEPENDANCE
              fac2=ZooGR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              fac1=dtdays*ZooGR(ng)*fac2
#else
              fac1=dtdays*ZooGR(ng)
#endif
!
              cff11=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt1)/                  &
     &             (K_Phy1(ng)+Bio(i,k,iPhyt1)*Bio(i,k,iPhyt1))
              cff31=1.0_r8/(1.0_r8+cff11)
              Bio(i,k,iPhyt1)=cff31*Bio(i,k,iPhyt1)
              Bio(i,k,iChlo1)=cff31*Bio(i,k,iChlo1)
# if defined PHYT2 || defined PHYT3
!
              cff12=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt2)/                  &
     &             (K_Phy2(ng)+Bio(i,k,iPhyt2)*Bio(i,k,iPhyt2))
              cff32=1.0_r8/(1.0_r8+cff12)
              Bio(i,k,iPhyt2)=cff32*Bio(i,k,iPhyt2)
              Bio(i,k,iChlo2)=cff32*Bio(i,k,iChlo2)
# endif
# ifdef PHYT3
!
              cff13=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt3)/                  &
     &             (K_Phy3(ng)+Bio(i,k,iPhyt3)*Bio(i,k,iPhyt3))
              cff33=1.0_r8/(1.0_r8+cff13)
              Bio(i,k,iPhyt3)=cff33*Bio(i,k,iPhyt3)
              Bio(i,k,iChlo3)=cff33*Bio(i,k,iChlo3)
# endif	
!
!  9. Compute phytoplankton assimilated to zooplankton and 
!       egested to small detritus.
!
              N_Flux_Assim1=Bio(i,k,iPhyt1)*cff11*ZooAE_N(ng)
              N_Flux_Egest1=Bio(i,k,iPhyt1)*cff11*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim1
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest1
# if defined PHYT2 || defined PHYT3
!
              N_Flux_Assim2=Bio(i,k,iPhyt2)*cff12*ZooAE_N(ng)
              N_Flux_Egest2=Bio(i,k,iPhyt2)*cff12*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest2
# endif
# ifdef PHYT3
!
              N_Flux_Assim3=Bio(i,k,iPhyt3)*cff13*ZooAE_N(ng)
              N_Flux_Egest3=Bio(i,k,iPhyt3)*cff13*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim3
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest3
# endif	

!
!  10. Compute phytoplankton respiration.	
!
#ifdef TDEPENDANCE
              fac4=PhyBR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
#else
              fac4=PhyBR_t(ng)
#endif
!
	        cff41=dtdays*PhyBR(ng)*fac4+dtdays*t_PPmax1*Att_N1*PhyPR(ng)
              N_Flux_Presp1=cff41*MAX(Bio(i,k,iPhyt1)-PhyMin1(ng),0.0_r8)
              Bio(i,k,iPhyt1)=Bio(i,k,iPhyt1)-N_Flux_Presp1
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Presp1
              Bio(i,k,iChlo1)=Bio(i,k,iChlo1)-                            &
     &                        cff41*MAX(Bio(i,k,iChlo1)-ChlMin1(ng),0.0_r8)
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Presp1*rOxNH4
# if defined PHYT2 || defined PHYT3
!
	        cff42=dtdays*PhyBR(ng)*fac4+dtdays*t_PPmax2*Att_N2*PhyPR(ng)
              N_Flux_Presp2=cff42*MAX(Bio(i,k,iPhyt2)-PhyMin2(ng),0.0_r8)
              Bio(i,k,iPhyt2)=Bio(i,k,iPhyt2)-N_Flux_Presp2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Presp2
              Bio(i,k,iChlo2)=Bio(i,k,iChlo2)-                            &
     &                        cff42*MAX(Bio(i,k,iChlo2)-ChlMin2(ng),0.0_r8)
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Presp2*rOxNH4
# endif
# ifdef PHYT3
!
	        cff43=dtdays*PhyBR(ng)*fac4+dtdays*t_PPmax3*Att_N3*PhyPR(ng)
              N_Flux_Presp3=cff43*MAX(Bio(i,k,iPhyt3)-PhyMin3(ng),0.0_r8)
              Bio(i,k,iPhyt3)=Bio(i,k,iPhyt3)-N_Flux_Presp3
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Presp3
              Bio(i,k,iChlo3)=Bio(i,k,iChlo3)-                            &
     &                        cff43*MAX(Bio(i,k,iChlo3)-ChlMin3(ng),0.0_r8)
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Presp3*rOxNH4
# endif	
!
!  11. Compute phytoplankton mortality (limited by a phytoplankton minimum).
!
#ifdef TDEPENDANCE
              fac3=PhyMR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
#else
              fac3=PhyMR_t(ng)
#endif
              cff2=dtdays*PhyMR(ng)*fac3
!
              N_Flux_Pmortal1=cff2*MAX(Bio(i,k,iPhyt1)-PhyMin1(ng),0.0_r8)
              Bio(i,k,iPhyt1)=Bio(i,k,iPhyt1)-N_Flux_Pmortal1
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal1
              Bio(i,k,iChlo1)=Bio(i,k,iChlo1)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo1)-ChlMin1(ng),0.0_r8)
#  ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)                               &
     &                       +PhyPN(ng)*(N_Flux_Egest1+N_Flux_Pmortal1)   &
     &                       +(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim1  
#  endif

# if defined PHYT2 || defined PHYT3
!
              N_Flux_Pmortal2=cff2*MAX(Bio(i,k,iPhyt2)-PhyMin2(ng),0.0_r8)
              Bio(i,k,iPhyt2)=Bio(i,k,iPhyt2)-N_Flux_Pmortal2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal2
              Bio(i,k,iChlo2)=Bio(i,k,iChlo2)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo2)-ChlMin2(ng),0.0_r8)
#  ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)                               &
     &                       +PhyPN(ng)*(N_Flux_Egest2+N_Flux_Pmortal2)   &
     &                       +(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim2  
#  endif
# endif

# ifdef PHYT3
!
              N_Flux_Pmortal3=cff2*MAX(Bio(i,k,iPhyt3)-PhyMin3(ng),0.0_r8)
              Bio(i,k,iPhyt3)=Bio(i,k,iPhyt3)-N_Flux_Pmortal3
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal3
              Bio(i,k,iChlo3)=Bio(i,k,iChlo3)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo3)-ChlMin3(ng),0.0_r8)
#  ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)                               &
     &                       +PhyPN(ng)*(N_Flux_Egest3+N_Flux_Pmortal3)   &
     &                       +(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim3  
#  endif
# endif	

            END DO
          END DO

! x ==================================================================
!    Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!    mortality to small detritus (rate: ZooMR), zooplankton ingestion
!    related excretion (rate: ZooER).
! x ==================================================================

          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  12. Compute zooplankton mortality and excretion.
!
              fac2=dtdays*ZooMR(ng)
              fac3=dtdays*ZooER(ng)
              cff2=fac2*Bio(i,k,iZoop)

              fac11=fac3*Bio(i,k,iPhyt1)*Bio(i,k,iPhyt1)/                  &
     &              (K_Phy1(ng)+Bio(i,k,iPhyt1)*Bio(i,k,iPhyt1))
              fac4=fac11
#if defined PHYT2 || defined PHYT3
              fac12=fac3*Bio(i,k,iPhyt2)*Bio(i,k,iPhyt2)/                  &
     &              (K_Phy2(ng)+Bio(i,k,iPhyt2)*Bio(i,k,iPhyt2))
              fac4=fac4+fac12
#endif
#ifdef PHYT3
              fac13=fac3*Bio(i,k,iPhyt3)*Bio(i,k,iPhyt3)/                  &
     &              (K_Phy3(ng)+Bio(i,k,iPhyt3)*Bio(i,k,iPhyt3))
              fac4=fac4+fac13
#endif
              cff3=fac4*ZooAE_N(ng)
!
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
              N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
#ifdef PHOSPHORUS
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zexcret
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+ZooPN(ng)*N_Flux_Zmortal
#endif
!
!  13. Compute zooplankton basal metabolism (limited by a zooplankton minimum).
!
              cff1=dtdays*ZooBM(ng)
!
              N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef PHOSPHORUS
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zmetabo
#endif
#ifdef OXYGEN
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
#endif
            END DO
          END DO

!
!  14. Compute coagulation of small detritus to large detritus.
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*Bio(i,k,iSDeN)
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)/(1.0_r8+cff1)
              N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+N_Flux_CoagD                             
#ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+PhyPN(ng)*N_Flux_CoagD
#endif
            END DO
          END DO

!
!  15. Compute detritus recycling to NH4, remineralization.
!
#ifdef OXYGEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
# ifdef TDEPENDANCE
              fac2=fac2*(RR_t(ng)**(Bio(i,k,iTemp)-20.0_r8))
# endif
              cff1=dtdays*SDeRRN(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
!
              cff1=dtdays*SDeRRP(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRP(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
# endif
            END DO
          END DO
#else
          cff1=dtdays*SDeRRN(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRN(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine

            END DO
          END DO
# ifdef PHOSPHORUS
!
          cff1=dtdays*SDeRRP(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRP(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux

            END DO
          END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
# ifdef TDEPENDANCE
              fac2=fac2*(1.05_r8**(Bio(i,k,iTemp)-20.0_r8))
# endif
              cff1=dtdays*H2SOR(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
              S_Flux=Bio(i,k,iH2S_)*cff1
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
              Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
            END DO
          END DO
#endif
#ifdef OXYGEN
!
!  16. Compute surface O2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
            SchmidtN_Ox=1638.0_r8-                                      &
     &                  Bio(i,k,itemp)*(81.83_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (1.483_r8-                      &
     &                                   Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
            SchmidtN_Ox=1953.4_r8-                                      &
     &                  Bio(i,k,itemp)*(128.0_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (3.9918_r8-                     &
     &                                   Bio(i,k,itemp)*0.050091_r8))
# endif

            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
            O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
            O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
     &                     O2_Flux*Hz_inv(i,k)

          END DO
#endif
!
!  17. Vertical sinking terms.
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
            cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
            cff3=115.0_r8/16.0_r8
            cff4=106.0_r8/16.0_r8
# endif
            IF ((ibio.eq.iPhyt1).or.                                    &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
#  ifdef OXYGEN
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
		  Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
#  endif
# else
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
#  ifdef OXYGEN
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
		   Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
#  endif
# endif
              END DO
            END IF
# ifdef PHOSPHORUS
            IF ((ibio.eq.iLDeP).or.                                     &
     &          (ibio.eq.iSDeP)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
              END DO
            END IF
            IF (ibio.eq.iPhyt1)THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
              END DO
            END IF
# endif
# if defined H2S && defined OXYGEN
            IF ((ibio.eq.iPhyt1).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
                Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                          &
     &                         0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
              END DO
            END IF
# endif
#endif
          END DO SINK_LOOP

#ifdef BIO_SED_CONSTANT
!
!  18. Elution and oxygen consumption from/by sediment. 
!
          fac1=dtdays
          fac2=1.0_r8
          fac3=1.0_r8
# ifdef NPFLUX_BY_DO
          cff2=R_NH4f_max(ng)/14.0_r8
#  ifdef PHOSPHORUS
          cff3=R_PO4f_max(ng)/31.0_r8
#  endif
# else
          cff2=R_NH4f(ng)/14.0_r8     !NH4 elution flux from sediment
#  ifdef PHOSPHORUS
          cff3=R_PO4f(ng)/31.0_r8     !PO4 elution flux from sediment
#  endif
# endif
   	   DO i=Istr,Iend
# ifdef OXYGEN
            cff1=R_SODf(ng)/mol2g_O2    !SOD flux
#  ifdef OXYGEN_SPATIAL_DIST1
		IF(FSOD_distribution(i,j)==4) THEN
			cff1=R_SODfa(ng)/mol2g_O2   ! at area A
		ELSE IF(FSOD_distribution(i,j)==3) THEN
			cff1=R_SODfb(ng)/mol2g_O2   ! at area B
		ELSE IF(FSOD_distribution(i,j)==2) THEN
			cff1=R_SODfc(ng)/mol2g_O2   ! at area C
		ELSE IF(FSOD_distribution(i,j)==1) THEN
		       cff1=R_SODfd(ng)/mol2g_O2  ! at area D (at the head of Osaka Bay)
		ENDIF

	      cff1=cff1*Alpha_SODf(ng)
#  endif
# endif
            cff=fac1*Hz_inv(i,1)

# if defined PCE_PARAM || defined PCE_2012_inoue || defined TDEPENDANCE
            fac2=t_SODf(ng)**(Bio(i,1,itemp)-20.0_r8)
# endif
# ifdef NPFLUX_BY_DO
            fac4=K_DO_npflux(ng)/mol2g_O2*1000.0_r8
	     Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
            fac5=MAX(Bio(i,1,iOxyg),0.0_r8)/fac4
            fac3=1.0_r8/(1.0_r8+fac5)
# endif
            Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff*cff2*fac3
# ifdef PHOSPHORUS
            Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff*cff3*fac3
# endif

# ifdef OXYGEN 
#  ifdef OXYGEN_SPATIAL_DIST1
!
            fac6=MIN(Bio(i,1,iOxyg),cff*cff1*fac2)
            cff4=MAX(fac6,0.0_r8)
            Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff4
	      Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)

#  elif defined OXYGEN_SPATIAL_DIST2
!
            fac6 = FSOD_distribution(i,j)
            cff4 = Bio(i,1,itemp)*95.52_r8
            cff4 = MIN(cff4-fac6,0.0_r8)
            Bio(i,1,iOxyg) = Bio(i,1,iOxyg)-cff4*cff
#  else
!
            fac6=MIN(Bio(i,1,iOxyg),cff*cff1*fac2)
            cff4=MAX(fac6,0.0_r8)
            Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff4
	      Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
#  endif
# endif

#  ifdef H2S
            cff5=MIN(fac6,0.0_r8)
            Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-cff5*rOxH2S
#  endif

          END DO

#endif

        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
#ifdef BLOWINGUP_CHECKER
              val(ibio)=val(ibio)+t(i,j,k,nnew,ibio)
#endif
            END DO
          END DO
        END DO
      END DO J_LOOP

#ifdef BLOWINGUP_CHECKER
!
!  If blowing-up, set exit_flag to stop computations. (okada)
!
      DO itrc=1,NBT
        ibio=idbio(itrc)
        WRITE (valchar,'(1pe8.1)') val(ibio)
        DO ii=1,8
          IF ((valchar(ii:ii).eq.'N').or.(valchar(ii:ii).eq.'n').or.    &
     &        (valchar(ii:ii).eq.'*')) THEN
            IF (Master) WRITE (stdout,100) ibio
            exit_flag=1
          END IF
        END DO
      END DO
 100  FORMAT ('Blowing-up in fennel.h, varid=',i2)
#endif

      RETURN
      END SUBROUTINE biology_tile
