#include "cppdefs.h"
      MODULE ini_param_mod
#if defined NONLINEAR && defined ADJUST_PARAM
!
!svn $Id$
!=======================================================================
!  (c) 2016-01-20 Teruhisa Okada                                       !
!=======================================================================
!                                                                      !
!  This routine initializes parameter to adjust.                       !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: ini_param
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_param (ng)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
      USE mod_biology
# ifdef CHECKER
      USE mod_parallel
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: it
      integer :: Tindex
      logical, dimension(Nparam(ng)) :: Lrestart
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
# endif
!
!  Initialize time stepping indices and counters.
!
      Tindex=1
      DO it=1,Nparam(ng)
        Lrestart(it)=.FALSE.
      END DO
!
!  Initialize parameters
!
# ifdef CHECKER
      if (master) print*, 'ini_param', Lparam(iwPhy,ng), Ldef_p(iwPhy,ng)
      if (master) print*, 'ini_param', OCEAN(ng)%p(:,iwPhy)
# endif
# ifdef RESTART_PARAM
      DO it=1,Nparam(ng)
        IF (.not.Ldef_p(it,ng).and.Lparam(it,ng)) THEN
          Lrestart(it)=.TRUE.
        END IF
      END DO
# endif
# ifdef EXP_PARAM
      DO it=1,Nparam(ng)
        IF (.not.Lrestart(it)) THEN
          OCEAN(ng)%p(Tindex,it)=0.0_r8
        END IF
      END DO
# else


      IF (.not.Lrestart(iK_PO4)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4)=K_PO4(ng)
      END IF
      IF (.not.Lrestart(ig_max)) THEN
        OCEAN(ng)%p(Tindex,ig_max)=g_max(ng)
      END IF
      IF (.not.Lrestart(it_opt)) THEN
        OCEAN(ng)%p(Tindex,it_opt)=t_opt(ng)
      END IF
      IF (.not.Lrestart(iI_opt)) THEN
        OCEAN(ng)%p(Tindex,iI_opt)=I_opt(ng)
      END IF
      IF (.not.Lrestart(iK_NO3)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3)=K_NO3(ng)
      END IF
      IF (.not.Lrestart(iK_NH4)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4)=K_NH4(ng)
      END IF
      IF (.not.Lrestart(iK_Phy)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy)=K_Phy(ng)
      END IF
      IF (.not.Lrestart(iChl2C_m)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m)=Chl2C_m(ng)
      END IF
      IF (.not.Lrestart(iChlMin)) THEN
        OCEAN(ng)%p(Tindex,iChlMin)=ChlMin(ng)
      END IF
      IF (.not.Lrestart(iPhyMin)) THEN
        OCEAN(ng)%p(Tindex,iPhyMin)=PhyMin(ng)
      END IF
      IF (.not.Lrestart(iPhyIS)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS)=PhyIS(ng)
      END IF
! ====================================================
!   New parameters for the model (Inoue, 2019)
! ====================================================
      IF (.not.Lrestart(iK_PO4_1)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4_1)=K_PO4_1(ng)
      END IF
      IF (.not.Lrestart(iK_PO4_2)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4_2)=K_PO4_2(ng)
      END IF
      IF (.not.Lrestart(iK_PO4_3)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4_3)=K_PO4_3(ng)
      END IF
      IF (.not.Lrestart(ig_max1)) THEN
        OCEAN(ng)%p(Tindex,ig_max1)=g_max1(ng)
      END IF
      IF (.not.Lrestart(ig_max2)) THEN
        OCEAN(ng)%p(Tindex,ig_max2)=g_max2(ng)
      END IF
      IF (.not.Lrestart(ig_max3)) THEN
        OCEAN(ng)%p(Tindex,ig_max3)=g_max3(ng)
      END IF
      IF (.not.Lrestart(it_opt1)) THEN
        OCEAN(ng)%p(Tindex,it_opt1)=t_opt1(ng)
      END IF
      IF (.not.Lrestart(it_opt2)) THEN
        OCEAN(ng)%p(Tindex,it_opt2)=t_opt2(ng)
      END IF
      IF (.not.Lrestart(it_opt3)) THEN
        OCEAN(ng)%p(Tindex,it_opt3)=t_opt3(ng)
      END IF
      IF (.not.Lrestart(it_max1)) THEN
        OCEAN(ng)%p(Tindex,it_max1)=t_max1(ng)
      END IF
      IF (.not.Lrestart(it_max2)) THEN
        OCEAN(ng)%p(Tindex,it_max2)=t_max2(ng)
      END IF
      IF (.not.Lrestart(it_max3)) THEN
        OCEAN(ng)%p(Tindex,it_max3)=t_max3(ng)
      END IF
      IF (.not.Lrestart(iI_opt1)) THEN
        OCEAN(ng)%p(Tindex,iI_opt1)=I_opt1(ng)
      END IF
      IF (.not.Lrestart(iI_opt2)) THEN
        OCEAN(ng)%p(Tindex,iI_opt2)=I_opt3(ng)
      END IF
      IF (.not.Lrestart(iI_opt3)) THEN
        OCEAN(ng)%p(Tindex,iI_opt3)=I_opt3(ng)
      END IF
      IF (.not.Lrestart(iK_NO3_1)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3_1)=K_NO3_1(ng)
      END IF
      IF (.not.Lrestart(iK_NO3_2)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3_2)=K_NO3_2(ng)
      END IF
      IF (.not.Lrestart(iK_NO3_3)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3_3)=K_NO3_3(ng)
      END IF
      IF (.not.Lrestart(iK_NH4_1)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4_1)=K_NH4_1(ng)
      END IF
      IF (.not.Lrestart(iK_NH4_2)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4_2)=K_NH4_2(ng)
      END IF
      IF (.not.Lrestart(iK_NH4_3)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4_3)=K_NH4_3(ng)
      END IF
      IF (.not.Lrestart(iK_Phy1)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy1)=K_Phy1(ng)
      END IF
      IF (.not.Lrestart(iK_Phy2)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy2)=K_Phy2(ng)
      END IF
      IF (.not.Lrestart(iK_Phy3)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy3)=K_Phy3(ng)
      END IF
      IF (.not.Lrestart(iChl2C_m1)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m1)=Chl2C_m1(ng)
      END IF
      IF (.not.Lrestart(iChl2C_m2)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m2)=Chl2C_m2(ng)
      END IF
      IF (.not.Lrestart(iChl2C_m3)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m3)=Chl2C_m3(ng)
      END IF
      IF (.not.Lrestart(iChlMin1)) THEN
        OCEAN(ng)%p(Tindex,iChlMin1)=ChlMin1(ng)
      END IF
      IF (.not.Lrestart(iChlMin2)) THEN
        OCEAN(ng)%p(Tindex,iChlMin2)=ChlMin2(ng)
      END IF
      IF (.not.Lrestart(iChlMin3)) THEN
        OCEAN(ng)%p(Tindex,iChlMin3)=ChlMin3(ng)
      END IF
      IF (.not.Lrestart(iPhyMin1)) THEN
        OCEAN(ng)%p(Tindex,iPhyMin1)=PhyMin1(ng)
      END IF
      IF (.not.Lrestart(iPhyMin2)) THEN
        OCEAN(ng)%p(Tindex,iPhyMin2)=PhyMin2(ng)
      END IF
      IF (.not.Lrestart(iPhyMin3)) THEN
        OCEAN(ng)%p(Tindex,iPhyMin3)=PhyMin3(ng)
      END IF
      IF (.not.Lrestart(iPhyIS1)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS1)=PhyIS1(ng)
      END IF
      IF (.not.Lrestart(iPhyIS2)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS2)=PhyIS2(ng)
      END IF
      IF (.not.Lrestart(iPhyIS3)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS3)=PhyIS3(ng)
      END IF
      IF (.not.Lrestart(it_Phy1)) THEN
        OCEAN(ng)%p(Tindex,it_Phy1)=t_Phy1(ng)
      END IF
      IF (.not.Lrestart(it_Phy2)) THEN
        OCEAN(ng)%p(Tindex,it_Phy2)=t_Phy2(ng)
      END IF
      IF (.not.Lrestart(it_Phy3)) THEN
        OCEAN(ng)%p(Tindex,it_Phy3)=t_Phy3(ng)

      IF (.not.Lrestart(ibeta01)) THEN
        OCEAN(ng)%p(Tindex,beta01)=beta01(ng)
      END IF

      IF (.not.Lrestart(ibeta02)) THEN
        OCEAN(ng)%p(Tindex,beta02)=beta02(ng)
      END IF

      IF (.not.Lrestart(ibeta03)) THEN
        OCEAN(ng)%p(Tindex,beta03)=beta03(ng)
      END IF

! ======================================================

      END IF
      IF (.not.Lrestart(ibeta1)) THEN
        OCEAN(ng)%p(Tindex,beta1)=beta1(ng)
      END IF
      IF (.not.Lrestart(ibeta2)) THEN
        OCEAN(ng)%p(Tindex,beta2)=beta2(ng)
      END IF
      IF (.not.Lrestart(ibeta)) THEN
        OCEAN(ng)%p(Tindex,beta)=beta(ng)
      END IF
      IF (.not.Lrestart(iPhyMR_t)) THEN
        OCEAN(ng)%p(Tindex,iPhyMR_t)=PhyMR_t(ng)
      END IF
      IF (.not.Lrestart(iPhyMR)) THEN
        OCEAN(ng)%p(Tindex,iPhyMR)=PhyMR(ng)
      END IF
      IF (.not.Lrestart(iPhyPR)) THEN
        OCEAN(ng)%p(Tindex,iPhyPR)=PhyPR(ng)
      END IF
      IF (.not.Lrestart(iPhyBR)) THEN
        OCEAN(ng)%p(Tindex,iPhyBR)=PhyBR(ng)
      END IF
      IF (.not.Lrestart(iPhyBR_t)) THEN
        OCEAN(ng)%p(Tindex,iPhyBR_t)=PhyBR_t(ng)
      END IF
      IF (.not.Lrestart(iAttSW)) THEN
        OCEAN(ng)%p(Tindex,iAttSW)=AttSW(ng)
      END IF
      IF (.not.Lrestart(iAttChl)) THEN
        OCEAN(ng)%p(Tindex,iAttChl)=AttChl(ng)
      END IF
      IF (.not.Lrestart(iPARfrac)) THEN
        OCEAN(ng)%p(Tindex,iPARfrac)=PARfrac(ng)
      END IF
      IF (.not.Lrestart(iPhyIP)) THEN
        OCEAN(ng)%p(Tindex,iPhyIP)=PhyIP(ng)
      END IF
      IF (.not.Lrestart(iwPhy)) THEN
        OCEAN(ng)%p(Tindex,iwPhy)=wPhy(ng)
      END IF
      IF (.not.Lrestart(iPhyPN)) THEN
        OCEAN(ng)%p(Tindex,iPhyPN)=PhyPN(ng)
      END IF
      IF (.not.Lrestart(iPhyCN)) THEN
        OCEAN(ng)%p(Tindex,iPhyCN)=PhyCN(ng)
      END IF
      IF (.not.Lrestart(iZooGR_t)) THEN
        OCEAN(ng)%p(Tindex,iZooGR_t)=ZooGR_t(ng)
      END IF
      IF (.not.Lrestart(iZooGR)) THEN
        OCEAN(ng)%p(Tindex,iZooGR)=ZooGR(ng)
      END IF
      IF (.not.Lrestart(iZooAE_N)) THEN
        OCEAN(ng)%p(Tindex,iZooAE_N)=ZooAE_N(ng)
      END IF
      IF (.not.Lrestart(iZooBM)) THEN
        OCEAN(ng)%p(Tindex,iZooBM)=ZooBM(ng)
      END IF
      IF (.not.Lrestart(iZooER)) THEN
        OCEAN(ng)%p(Tindex,iZooER)=ZooER(ng)
      END IF
      IF (.not.Lrestart(iZooMR)) THEN
        OCEAN(ng)%p(Tindex,iZooMR)=ZooMR(ng)
      END IF
      IF (.not.Lrestart(iZooPN)) THEN
        OCEAN(ng)%p(Tindex,iZooPN)=ZooPN(ng)
      END IF
      IF (.not.Lrestart(iZooCN)) THEN
        OCEAN(ng)%p(Tindex,iZooCN)=ZooCN(ng)
      END IF
      IF (.not.Lrestart(iLDeRRN)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRN)=LDeRRN(ng)
      END IF
      IF (.not.Lrestart(iSDeRRN)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRN)=SDeRRN(ng)
      END IF
      IF (.not.Lrestart(iLDeRRP)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRP)=LDeRRP(ng)
      END IF
      IF (.not.Lrestart(iSDeRRP)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRP)=SDeRRP(ng)
      END IF
      IF (.not.Lrestart(iwLDet)) THEN
        OCEAN(ng)%p(Tindex,iwLDet)=wLDet(ng)
      END IF
      IF (.not.Lrestart(iwSDet)) THEN
        OCEAN(ng)%p(Tindex,iwSDet)=wSDet(ng)
      END IF
      IF (.not.Lrestart(iK_Nitri)) THEN
        OCEAN(ng)%p(Tindex,iK_Nitri)=K_nitri(ng)
      END IF
      IF (.not.Lrestart(iDenitR)) THEN
        OCEAN(ng)%p(Tindex,iDenitR)=DenitR(ng)
      END IF
      IF (.not.Lrestart(iK_Denit)) THEN
        OCEAN(ng)%p(Tindex,iK_Denit)=K_Denit(ng)
      END IF
      IF (.not.Lrestart(iDenitR_t)) THEN
        OCEAN(ng)%p(Tindex,iDenitR_t)=DenitR_t(ng)
      END IF
      IF (.not.Lrestart(iI_thNH4)) THEN
        OCEAN(ng)%p(Tindex,iI_thNH4)=I_thNH4(ng)
      END IF
      IF (.not.Lrestart(iD_p5NH4)) THEN
        OCEAN(ng)%p(Tindex,iD_p5NH4)=D_p5NH4(ng)
      END IF
      IF (.not.Lrestart(iNitriR)) THEN
        OCEAN(ng)%p(Tindex,iNitriR)=NitriR(ng)
      END IF
      IF (.not.Lrestart(iNitriR_t)) THEN
        OCEAN(ng)%p(Tindex,iNitriR_t)=NitriR_t(ng)
      END IF
      IF (.not.Lrestart(iR_NH4f)) THEN
        OCEAN(ng)%p(Tindex,iR_NH4f)=R_NH4f(ng)
      END IF
      IF (.not.Lrestart(iR_NH4f_m)) THEN
        OCEAN(ng)%p(Tindex,iR_NH4f_m)=R_NH4f_max(ng)
      END IF
      IF (.not.Lrestart(iK_DO)) THEN
        OCEAN(ng)%p(Tindex,iK_DO)=K_DO(ng)
      END IF
      IF (.not.Lrestart(iRR_t)) THEN
        OCEAN(ng)%p(Tindex,iRR_t)=RR_t(ng)
      END IF
      IF (.not.Lrestart(iK_DO_npf)) THEN
        OCEAN(ng)%p(Tindex,iK_DO_npf)=K_DO_npflux(ng)
      END IF
      IF (.not.Lrestart(iR_SODf)) THEN
        OCEAN(ng)%p(Tindex,iR_SODf)=R_SODf(ng)
      END IF
      IF (.not.Lrestart(it_SODf)) THEN
        OCEAN(ng)%p(Tindex,it_SODf)=t_SODf(ng)
      END IF
      IF (.not.Lrestart(iR_SODfa)) THEN
        OCEAN(ng)%p(Tindex,iR_SODfa)=R_SODfa(ng)
      END IF
      IF (.not.Lrestart(iR_SODfb)) THEN
        OCEAN(ng)%p(Tindex,iR_SODfb)=R_SODfb(ng)
      END IF
      IF (.not.Lrestart(iR_SODfc)) THEN
        OCEAN(ng)%p(Tindex,iR_SODfc)=R_SODfc(ng)
      END IF
      IF (.not.Lrestart(iR_SODfd)) THEN
        OCEAN(ng)%p(Tindex,iR_SODfd)=R_SODfd(ng)
      END IF
      IF (.not.Lrestart(iAlpha_SODf)) THEN
        OCEAN(ng)%p(Tindex,iAlpha_SODf)=Alpha_SODf(ng)
      END IF
      IF (.not.Lrestart(iR_PO4f)) THEN
        OCEAN(ng)%p(Tindex,iR_PO4f)=R_PO4f(ng)
      END IF
      IF (.not.Lrestart(iR_PO4f_m)) THEN
        OCEAN(ng)%p(Tindex,iR_PO4f_m)=R_PO4f_max(ng)
      END IF


# endif
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
          OCEAN(ng)%p(:,it)=OCEAN(ng)%p(Tindex,it)
        END IF
      END DO
# ifdef CHECKER
      if (master) print*, 'ini_param', OCEAN(ng)%p(:,iwPhy)
# endif
!
!  Set parameters standard deviation.
!
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
# ifdef EXP_PARAM
          OCEAN(ng)%e_p(:,it)=pfactor(it,ng)
# else
!         OCEAN(ng)%e_p(:,it)=(pmax(it,ng)-pmin(it,ng))*pfactor(it,ng)
          OCEAN(ng)%e_p(:,it)=OCEAN(ng)%p(Tindex,it)*pfactor(it,ng)
# endif
        END IF
      END DO
!
!  Set parameters normalization factor.
!
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
          OCEAN(ng)%b_p(:,it)=1.0_r8
        END IF
      END DO
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
# endif
      RETURN
      END SUBROUTINE ini_param
#endif
      END MODULE ini_param_mod
