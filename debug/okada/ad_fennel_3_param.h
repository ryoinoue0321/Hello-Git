#include "fennel_bs3.h"

#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
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
! not yet (okada)
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
            O2_Flux=cff3*(O2satu-Bio3(i,k,iOxyg))
!!          Bio3(i,k,iOxyg)=Bio(i,k,iOxyg)
!>          Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
!>   &                     O2_Flux*Hz_inv(i,k)
!
!=======================================================================
! Switch Back
!=======================================================================
!
!  Add in O2 gas exchange.
!
!>          tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)+                        &
!>   &                        tl_O2_Flux*Hz_inv(i,k)+                   &
!>   &                        O2_Flux*tl_Hz_inv(i,k)
            ad_Hz_inv(i,k)=ad_Hz_inv(i,k)+O2_Flux*ad_Bio(i,k,iOxyg)
            ad_O2_Flux=ad_O2_Flux+ad_Bio(i,k,iOxyg)*Hz_inv(i,k)

!>          tl_O2_Flux=tl_cff3*(O2satu-Bio3(i,k,iOxyg))+                &
!>   &                 cff3*(tl_O2satu-tl_Bio(i,k,iOxyg))
            ad_Bio(i,k,iOxyg)=ad_Bio(i,k,iOxyg)-cff3*ad_O2_Flux
            ad_O2satu=ad_O2satu+cff3*ad_O2_Flux
            ad_cff3=ad_cff3+ad_O2_Flux*(O2satu-Bio3(i,k,iOxyg))
            ad_O2_Flux=0.0_r8
!
!  Convert from ml/l to mmol/m3.
!
!>          tl_O2satu=l2mol*EXP(AA)*tl_AA
            ad_AA=ad_AA+l2mol*EXP(AA)*ad_O2satu
            ad_O2satu=0.0_r8
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
!>          tl_AA=tl_TS*(OA1+TS*(2.0_r8*OA2+TS*(3.0_r8*OA3+TS*          &
!>   &            (4.0_r8*OA4+TS*5.0_r8*OA5))))+                        &
!>   &            tl_Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+     &
!>   &            Bio(i,k,isalt)*(tl_TS*(OB1+TS*(2.0_r8*OB2+            &
!>   &                                           TS*3.0_r8*OB3)))+      &
!>   &            OC0*2.0_r8*Bio(i,k,isalt)*tl_Bio(i,k,isalt)
            ad_Bio(i,k,isalt)=ad_Bio(i,k,isalt)+                        &
     &                        OC0*2.0_r8*Bio(i,k,isalt)*ad_AA
            ad_TS=ad_TS+Bio(i,k,isalt)*(ad_AA*(OB1+TS*(2.0_r8*OB2+      &
     &                                                 TS*3.0_r8*OB3)))
            ad_Bio(i,k,isalt)=ad_Bio(i,k,isalt)+                        &
     &                        ad_AA*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))
            ad_TS=ad_TS+ad_AA*(OA1+TS*(2.0_r8*OA2+TS*(3.0_r8*OA3+TS*    &
     &                         (4.0_r8*OA4+TS*5.0_r8*OA5))))
            ad_AA=0.0_r8

!>          tl_TS=-tl_Bio(i,k,itemp)/(298.15_r8-Bio(i,k,itemp))-        &
!>   &             tl_Bio(i,k,itemp)/(273.15_r8+Bio(i,k,itemp))
# ifndef UV_FIXED_TL
            ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)-                        &
     &                        ad_TS/(273.15_r8+Bio(i,k,itemp))
            ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)-                        &
     &                        ad_TS/(298.15_r8-Bio(i,k,itemp))
# endif
            ad_TS=0.0_r8

!>          tl_cff3=-0.5_r8*cff3*tl_SchmidtN_Ox/SchmidtN_Ox
            ad_SchmidtN_Ox=ad_SchmidtN_Ox-                              &
                           0.5_r8*cff3*ad_cff3/SchmidtN_Ox
            ad_cff3=0.0_r8
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
!>          tl_SchmidtN_Ox=-tl_Bio(i,k,itemp)*                          &
!>   &                     (81.83_r8-Bio(i,k,itemp)*                    &
!>   &                     (2.966_r8-Bio(i,k,itemp)*0.024012_r8))
#  ifndef UV_FIXED_TL
            ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)-ad_SchmidtN_Ox*         &
     &                        (81.83_r8-Bio(i,k,itemp)*                 &
     &                        (2.966_r8-Bio(i,k,itemp)*0.024012_r8))
#  endif
            ad_SchmidtN_Ox=0.0_r8
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
!>          tl_SchmidtN_Ox=tl_Bio(i,k,itemp)*                           &
!>                         (128.0_r8-Bio(i,k,itemp)*                    &
!>   &                     (7.9836_r8-Bio(i,k,itemp)*0.150273_r8))
#  ifndef UV_FIXED_TL
            ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+ad_SchmidtN_Ox*         &
                              (128.0_r8-Bio(i,k,itemp)*                 &
     &                        (7.9836_r8-Bio(i,k,itemp)*0.150273_r8))
#  endif
            ad_SchmidtN_Ox=0.0_r8
# endif
          END DO
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
! not yet (okada)
#endif
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#ifdef OXYGEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!             ==========================================================
              fac1=MAX(Bio2(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
              cff1=dtdays*SDeRRN(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!!            Bio2(i,k,iSDeN)=Bio(i,k,iSDeN)
!!            Bio2(i,k,iLDeN)=Bio(i,k,iLDeN)
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
!!            Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!!            Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
              cff1=dtdays*SDeRRP(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRP(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!!            Bio2(i,k,iSDeP)=Bio(i,k,iSDeP)
!!            Bio2(i,k,iLDeP)=Bio(i,k,iLDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
!!            Bio2(i,k,iPO4_)=Bio(i,k,iPO4_)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
# endif
!             ==========================================================
# ifdef PHOSPHORUS
!>            tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+tl_P_Flux
              ad_P_Flux=ad_P_Flux+ad_Bio(i,k,iPO4_)

!>            tl_P_Flux=tl_Bio(i,k,iSDeP)*cff1+Bio(i,k,iSDeP)*tl_cff1+  &
!>   &                  tl_Bio(i,k,iLDeP)*cff3+Bio(i,k,iLDeP)*tl_cff3
              ad_cff3=ad_cff3+Bio(i,k,iLDeP)*ad_P_Flux
              ad_Bio(i,k,iLDeP)=ad_Bio(i,k,iLDeP)+ad_P_Flux*cff3
              ad_cff1=ad_cff1+Bio(i,k,iSDeP)*ad_P_Flux
              ad_Bio(i,k,iSDeP)=ad_Bio(i,k,iSDeP)+ad_P_Flux*cff1
              ad_P_Flux=0.0_r8

!>            tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)*cff4+                 &
!>   &                          Bio2(i,k,iLDeP)*tl_cff4
              ad_cff4=ad_cff4+Bio2(i,k,iLDeP)*ad_Bio(i,k,iLDeP)
              ad_Bio(i,k,iLDeP)=ad_Bio(i,k,iLDeP)*cff4

!>            tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)*cff2+                 &
!>   &                          Bio2(i,k,iSDeP)*tl_cff2
              ad_cff2=ad_cff2+Bio2(i,k,iSDeP)*ad_Bio(i,k,iSDeP)
              ad_Bio(i,k,iSDeP)=ad_Bio(i,k,iSDeP)*cff2

!>            tl_cff4=-cff4*cff4*tl_cff3
              ad_cff3=ad_cff3-cff4*cff4*ad_cff4
              ad_cff4=0.0_r8

!>            tl_cff3=dtdays*(tl_LDeRRP*fac2+LDeRRP(ng)*tl_fac2)
              ad_fac2=ad_fac2+dtdays*LDeRRP(ng)*ad_cff3
              ad_LDeRRP=ad_LDeRRP+dtdays*fac2*ad_cff3
              ad_cff3=0.0_r8

!>            tl_cff2=-cff2*cff2*tl_cff1
              ad_cff1=ad_cff1-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8

!>            tl_cff1=dtdays*(tl_SDeRRP*fac2+SDeRRP(ng)*tl_fac2)
              ad_fac2=ad_fac2+dtdays*SDeRRP(ng)*ad_cff1
              ad_SDeRRP=ad_SDeRRP+dtdays*fac2*ad_cff1
              ad_cff1=0.0_r8
# endif
!             ==========================================================
              cff1=dtdays*SDeRRN(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
!             ==========================================================

!>            tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                      &
!>   &                          tl_N_Flux_Remine*rOxNH4
              ad_N_Flux_Remine=ad_N_Flux_Remine-ad_Bio(i,k,iOxyg)*rOxNH4

!>            tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Remine
              ad_N_Flux_Remine=ad_N_Flux_Remine+ad_Bio(i,k,iNH4_)

!>            tl_N_Flux_Remine=tl_Bio(i,k,iSDeN)*cff1+                  &
!>   &                         Bio(i,k,iSDeN)*tl_cff1+                  &
!>   &                         tl_Bio(i,k,iLDeN)*cff3+                  &
!>   &                         Bio(i,k,iLDeN)*tl_cff3
              ad_cff3=ad_cff3+Bio(i,k,iLDeN)*ad_N_Flux_Remine
              ad_Bio(i,k,iLDeN)=ad_Bio(i,k,iLDeN)+ad_N_Flux_Remine*cff3
              ad_cff1=ad_cff1+Bio(i,k,iSDeN)*ad_N_Flux_Remine
              ad_Bio(i,k,iSDeN)=ad_Bio(i,k,iSDeN)+ad_N_Flux_Remine*cff1
              ad_N_Flux_Remine=0.0_r8

!>            tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)*cff4+                 &
!>   &                          Bio2(i,k,iLDeN)*tl_cff4
              ad_cff4=ad_cff4+Bio2(i,k,iLDeN)*ad_Bio(i,k,iLDeN)
              ad_Bio(i,k,iLDeN)=ad_Bio(i,k,iLDeN)*cff4

!>            tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2+                 &
!>   &                          Bio2(i,k,iSDeN)*tl_cff2
              ad_cff2=ad_cff2+Bio2(i,k,iSDeN)*ad_Bio(i,k,iSDeN)
              ad_Bio(i,k,iSDeN)=ad_Bio(i,k,iSDeN)*cff2

!>            tl_cff4=-cff4*cff4*tl_cff3
              ad_cff3=ad_cff3-cff4*cff4*ad_cff4
              ad_cff4=0.0_r8

!>            tl_cff3=dtdays*(tl_LDeRRN*fac2+LDeRRN(ng)*tl_fac2)
              ad_fac2=ad_fac2+dtdays*LDeRRN(ng)*ad_cff3
              ad_LDeRRN=ad_LDeRRN+dtdays*fac2*ad_cff3
              ad_cff3=0.0_r8

!>            tl_cff2=-cff2*cff2*tl_cff1
              ad_cff1=ad_cff1-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8

!>            tl_cff1=dtdays*(tl_SDeRRN*fac2+SDeRRN(ng)*tl_fac2)
              ad_fac2=ad_fac2+dtdays*SDeRRN(ng)*ad_cff1
              ad_SDeRRN=ad_SDeRRN+dtdays*fac2*ad_cff1
              ad_cff1=0.0_r8

!>            tl_fac2=(tl_fac1-fac2*(tl_K_DO+tl_fac1))/(K_DO(ng)+fac1)
              adfac=ad_fac2/(K_DO(ng)+fac1)
              ad_fac1=ad_fac1-fac2*adfac
              ad_K_DO=ad_K_DO-fac2*adfac
              ad_fac1=ad_fac1-adfac
              ad_fac2=0.0_r8

!>            tl_fac1=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)-6.0_r8))*     &
!>   &                tl_Bio(i,k,iOxyg)
              adfac=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)-6.0_r8))
              ad_Bio(i,k,iOxyg)=ad_Bio(i,k,iOxyg)+adfac*ad_fac1
              ad_fac1=0.0_r8
            END DO
          END DO
#else
! pass okada
#endif
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
!!            Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!!            Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!!            Bio1(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              N_Flux_CoagP=Bio2(i,k,iPhyt)*cff1
              N_Flux_CoagD=Bio2(i,k,iSDeN)*cff1
!!            Bio1(i,k,iLDeN)=Bio(i,k,iLDeN)
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
!>   &                       N_Flux_CoagP+N_Flux_CoagD
#ifdef PHOSPHORUS
!!            Bio1(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
!!            Bio1(i,k,iLDeP)=Bio(i,k,iLDeP)
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
#endif
!             ==========================================================
#ifdef PHOSPHORUS
!>            tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)+PhyPN(ng)*            &
!>   &                          (tl_N_Flux_CoagP+tl_N_Flux_CoagD)+      &
!>   &                          tl_PhyPN*(N_Flux_CoagP+N_Flux_CoagD)
              adfac=PhyPN(ng)*ad_Bio(i,k,iLDeP)
              ad_N_Flux_CoagD=ad_N_Flux_CoagD+adfac
              ad_N_Flux_CoagP=ad_N_Flux_CoagP+adfac
              ad_PhyPN=ad_PhyPN+                                        &
     &                 (N_Flux_CoagP+N_Flux_CoagD)*ad_Bio(i,k,iLDeP)

!>            tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)-                      &
!>   &                          (tl_PhyPN*N_Flux_CoagD+                 &
!>   &                           PhyPN(ng)*tl_N_Flux_CoagD)
              ad_N_Flux_CoagD=ad_N_Flux_CoagD-                          &
     &                        PhyPN(ng)*ad_Bio(i,k,iSDeP)
              ad_PhyPN=ad_PhyPN-N_Flux_CoagD*ad_Bio(i,k,iSDeP)
#endif
!>            tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)+                      &
!>   &                          tl_N_Flux_CoagP+tl_N_Flux_CoagD
              ad_N_Flux_CoagD=ad_N_Flux_CoagD+ad_Bio(i,k,iLDeN)
              ad_N_Flux_CoagP=ad_N_Flux_CoagP+ad_Bio(i,k,iLDeN)

!>            tl_N_Flux_CoagD=tl_Bio(i,k,iSDeN)*cff1+                   &
!>   &                        Bio2(i,k,iSDeN)*tl_cff1
              ad_cff1=ad_cff1+Bio2(i,k,iSDeN)*ad_N_Flux_CoagD
              ad_Bio(i,k,iSDeN)=ad_Bio(i,k,iSDeN)+ad_N_Flux_CoagD*cff1
              ad_N_Flux_CoagD=0.0_r8

!>            tl_N_Flux_CoagP=tl_Bio(i,k,iPhyt)*cff1+                   &
!>   &                        Bio2(i,k,iPhyt)*tl_cff1
              ad_cff1=ad_cff1+Bio2(i,k,iPhyt)*ad_N_Flux_CoagP
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+ad_N_Flux_CoagP*cff1
              ad_N_Flux_CoagP=0.0_r8

!>            tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2+                 &
!>   &                          Bio1(i,k,iSDeN)*tl_cff2
              ad_cff2=ad_cff2+Bio1(i,k,iSDeN)*ad_Bio(i,k,iSDeN)
              ad_Bio(i,k,iSDeN)=ad_Bio(i,k,iSDeN)*cff2

!>            tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)*cff2+                 &
!>   &                          Bio1(i,k,iChlo)*tl_cff2
              ad_cff2=ad_cff2+Bio1(i,k,iChlo)*ad_Bio(i,k,iChlo)
              ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)*cff2

!>            tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)*cff2+                 &
!>   &                          Bio1(i,k,iPhyt)*tl_cff2
              ad_cff2=ad_cff2+Bio1(i,k,iPhyt)*ad_Bio(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)*cff2

!>            tl_cff2=-cff2*cff2*tl_cff1
              ad_cff1=ad_cff1-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8

!>            tl_cff1=tl_fac1*(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))+        &
!>   &                fac1*(tl_Bio(i,k,iSDeN)+tl_Bio(i,k,iPhyt))
              adfac=fac1*ad_cff1
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+adfac
              ad_Bio(i,k,iSDeN)=ad_Bio(i,k,iSDeN)+adfac
              ad_fac1=ad_fac1+(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))*ad_cff1
              ad_cff1=0.0_r8
            END DO
          END DO
!>        tl_fac1=dtdays*tl_CoagR
          ad_CoagR=ad_CoagR+dtdays*ad_fac1
          ad_fac1=0.0_r8
