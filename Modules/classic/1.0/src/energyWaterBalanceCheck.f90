!> \file
!> Checks for energy and water balance closure over modelled
!! area.
!! @author D. Verseghy
!
subroutine energyWaterBalanceCheck (ISTEP, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, & ! Formerly CLASSZ
                                    WTVSTP, WTSSTP, WTGSTP, &
                                    FSGV, FLGV, HFSC, HEVC, HMFC, HTCC, &
                                    FSGS, FLGS, HFSS, HEVS, HMFN, HTCS, &
                                    FSGG, FLGG, HFSG, HEVG, HMFG, HTC, &
                                    PCFC, PCLC, QFCF, QFCL, ROFC, WTRC, &
                                    PCPN, QFN, ROFN, WTRS, PCPG, QFG, &
                                    QFC, ROF, WTRG, CMAI, RCAN, SCAN, &
                                    TCAN, SNO, WSNOW, TSNOW, THLIQ, THICE, &
                                    HCPS, THPOR, DELZW, TBAR, ZPOND, TPOND, &
                                    DELZ, FCS, FGS, FC, FG, &
                                    IL1, IL2, ILG, IG, N)
  !
  !     * JAN 06/09 - D.VERSEGHY. MORE VARIABLES IN PRINT STATEMENTS
  !     *                         SLIGHTLY INCREASED ACCURACY LIMITS.
  !     * NOV 10/06 - D.VERSEGHY. CHECK THAT SUMS OF ENERGY AND WATER
  !     *                         FLUXES FOR CANOPY, SNOW AND SOIL MATCH
  !     *                         CHANGES IN HEAT AND WATER STORAGE OVER
  !     *                         CURRENT TIMESTEP.
  !
  use generalUtilsGAT,   only : abandonCell
  use classicParams,  only : DELT, TFREZ, HCPW, HCPICE, HCPSND, SPHW, &
                             SPHICE, SPHVEG, RHOW, RHOICE

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in)  :: ISTEP    !< Flag indicating position at beginning or end of
  !< time step
  integer, intent(in)  :: IL1, IL2, ILG, IG, N
  integer              :: I, J
  !
  !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
  !     * BALANCES.
  !
  real, intent(inout) :: CTVSTP(ILG)  !< Change in internal energy of vegetation over
  !< current time step \f$[W m^{-2}] \f$
  real, intent(inout) :: CTSSTP(ILG)  !< Change in internal energy of snow pack over
  !< current time step \f$[W m^{-2}] \f$
  real, intent(inout) :: CT1STP(ILG)  !< Change in internal energy of first soil layer
  !< over current time step \f$[W m^{-2}] \f$
  real, intent(inout) :: CT2STP(ILG)  !< Change in internal energy of second soil layer
  !< over current time step \f$[W m^{-2}] \f$
  real, intent(inout) :: CT3STP(ILG)  !< Change in internal energy of third soil layer
  !< over current time step \f$[W m^{-2}] \f$
  real, intent(inout) :: WTVSTP(ILG)  !< Change in vegetation mass over current time step \f$[kg m^{-2}]\f$
  real, intent(inout) :: WTSSTP(ILG)  !< Change in snow mass over current time step \f$[kg m^{-2}]\f$
  real, intent(inout) :: WTGSTP(ILG)  !< Change in soil water storage over current time step \f$[kg m^{-2}]\f$
  !
  real :: QSUMV, QSUMS, QSUM1, QSUM2, QSUM3, WSUMV, WSUMS, WSUMG
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FSGV  (ILG)      !< Diagnosed net shortwave radiation on
  !< vegetation canopy \f$[W m^{-2}] (K_{*,c})\f$
  real, intent(in) :: FLGV  (ILG)      !< Diagnosed net longwave radiation on
  !< vegetation canopy \f$[W m^{-2}] (L_{*,c})\f$
  real, intent(in) :: HFSC  (ILG)      !< Diagnosed sensible heat flux on vegetation
  !< canopy \f$[W m^{-2}] (Q_{H,c})\f$
  real, intent(in) :: HEVC  (ILG)      !< Diagnosed latent heat flux on vegetation
  !< canopy \f$[W m^{-2}] (Q_{E,c})\f$
  real, intent(in) :: HMFC  (ILG)      !< Diagnosed energy associated with phase
  !< change of water on vegetation \f$[W m^{-2}] (Q_{M,c})\f$
  real, intent(in) :: HTCC  (ILG)      !< Diagnosed internal energy change of
  !< vegetation canopy due to conduction and/or
  !< change in mass \f$[W m^{-2}] (Q_{I,c})\f$
  real, intent(in) :: FSGS  (ILG)      !< Diagnosed net shortwave radiation at snow
  !< surface \f$[W m^{-2}] (K_{*,s})\f$
  real, intent(in) :: FLGS  (ILG)      !< Diagnosed net longwave radiation at snow
  !< surface \f$[W m^{-2}] (L_{*,s})\f$
  real, intent(in) :: HFSS  (ILG)      !< Diagnosed sensible heat flux at snow
  !< surface \f$[W m^{-2}] (Q_{H,s})\f$
  real, intent(in) :: HEVS  (ILG)      !< Diagnosed latent heat flux at snow surface
  !< \f$[W m^{-2}] (Q_{E,s})\f$
  real, intent(in) :: HMFN  (ILG)      !< Diagnosed energy associated with phase
  !< change of water in snow pack \f$[W m^{-2}] (Q_{M,s})\f$
  real, intent(in) :: HTCS  (ILG)      !< Diagnosed internal energy change of snow
  !< pack due to conduction and/or change in
  !< mass \f$[W m^{-2}] (Q_{I,s})\f$
  real, intent(in) :: FSGG  (ILG)      !< Diagnosed net shortwave radiation at soil
  !< surface \f$[W m^{-2}] (K_{*,g})\f$
  real, intent(in) :: FLGG  (ILG)      !< Diagnosed net longwave radiation at soil
  !< surface \f$[W m^{-2}] (L_{*,g})\f$
  real, intent(in) :: HFSG  (ILG)      !< Diagnosed sensible heat flux at soil
  !< surface \f$[W m^{-2}] (Q_{H,g})\f$
  real, intent(in) :: HEVG  (ILG)      !< Diagnosed latent heat flux at soil surface
  !< \f$[W m^{-2}] (Q_{E,g})\f$
  real, intent(in) :: HMFG  (ILG,IG)   !< Diagnosed energy associated with phase
  !< change of water in soil layers \f$[W m^{-2}] (Q_{M,g})\f$
  real, intent(in) :: HTC   (ILG,IG)   !< Diagnosed internal energy change of soil
  !< layer due to conduction and/or change in
  !< mass \f$[W m^{-2}] (Q_{I,g})\f$
  real, intent(in) :: PCFC  (ILG)      !< Diagnosed frozen precipitation intercepted
  !< by vegetation \f$[kg m^{-2} s^{-1}] (P_{f,c})\f$
  real, intent(in) :: PCLC  (ILG)      !< Diagnosed liquid precipitation intercepted
  !< by vegetation \f$[kg m^{-2} s^{-1}] (P_{l,c})\f$
  real, intent(in) :: QFCF  (ILG)      !< Diagnosed vapour flux from frozen water on
  !< vegetation \f$[kg m^{-2} s^{-1}] (E_{f,c})\f$
  real, intent(in) :: QFCL  (ILG)      !< Diagnosed vapour flux from liquid water on
  !< vegetation \f$[kg m^{-2} s^{-1}] (E_{l,c})\f$
  real, intent(in) :: ROFC  (ILG)      !< Liquid/frozen water runoff from vegetation
  !< \f$[kg m^{-2} s^{-1}] (R_c)\f$
  real, intent(in) :: WTRC  (ILG)      !< Diagnosed water transferred off the
  !< vegetation canopy \f$[kg m^{-2} s^{-1}] (A_c)\f$
  real, intent(in) :: PCPN  (ILG)      !< Diagnosed precipitation incident on snow
  !< pack \f$[kg m^{-2} s^{-1}] (P_s)\f$
  real, intent(in) :: QFN   (ILG)      !< Diagnosed water vapour flux from snow pack
  !< \f$[kg m^{-2} s^{-1}] (E_s)\f$
  real, intent(in) :: ROFN  (ILG)      !< Liquid water runoff from snow pack
  !< \f$[kg m^{-2} s^{-1}] (R_s)\f$
  real, intent(in) :: WTRS  (ILG)      !< Diagnosed water transferred into or out of
  !< the snow pack \f$[kg m^{-2} s^{-1}] (A_s)\f$
  real, intent(in) :: PCPG  (ILG)      !< Diagnosed precipitation incident on ground
  !< \f$[kg m^{-2} s^{-1}] (P_g)\f$
  real, intent(in) :: QFG   (ILG)      !< Diagnosed water vapour flux from ground
  !< surface \f$[kg m^{-2} s^{-1}] (E_g)\f$
  real, intent(in) :: QFC   (ILG,IG)   !< Diagnosed vapour flux from transpiration
  !< over modelled area \f$[W m^{-2}] (E_c)\f$
  real, intent(in) :: ROF   (ILG)      !< Total runoff from soil \f$[kg m^{-2} s^{-1}] (R_g)\f$
  real, intent(in) :: WTRG  (ILG)      !< Diagnosed water transferred into or out of
  !< the soil \f$[kg m^{-2} s^{-1}] (A_g)\f$
  real, intent(in) :: CMAI  (ILG)      !< Current mass of vegetation canopy \f$[kg m^{-2}] (W_c)\f$
  real, intent(in) :: RCAN  (ILG)      !< Intercepted liquid water stored on canopy
  !< \f$[kg m^{-2}] (W_{l,c})\f$
  real, intent(in) :: SCAN  (ILG)      !< Intercepted frozen water stored on canopy
  !< \f$[kg m^{-2}] (W_{f,c})\f$
  real, intent(in) :: TCAN  (ILG)      !< Vegetation canopy temperature \f$[K] (T_c)\f$
  real, intent(in) :: SNO   (ILG)      !< Mass of snow pack \f$[kg m^{-2}] (W_s)\f$
  real, intent(in) :: WSNOW (ILG)      !< Liquid water content of snow pack \f$[kg m^{-2}] (W_{l,s})\f$
  real, intent(in) :: TSNOW (ILG)      !< Snowpack temperature \f$[K] (T_s)\f$
  real, intent(in) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !< layers \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(in) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !< layers \f$[m^3 m^{-3}] (\theta_f)\f$
  real, intent(in) :: HCPS  (ILG,IG)   !< Volumetric heat capacity of soil particles
  !< \f$[J m^{-3}] (C_g)\f$
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable thickness of soil layer \f$[m] (\Delta_zw)\f$
  real, intent(in) :: DELZ  (IG)       !< Total thickness of soil layer \f$[m] (\Delta_z)\f$
  real, intent(in) :: TBAR  (ILG,IG)   !< Temperature of soil layers \f$[K] (T_g)\f$
  real, intent(in) :: ZPOND (ILG)      !< Depth of ponded water on surface \f$[m] (z_p)\f$
  real, intent(in) :: TPOND (ILG)      !< Total thickness of soil layer \f$[m] (\Delta_z)\f$
  real, intent(in) :: FCS   (ILG)      !< Fractional coverage of vegetation over snow on modelled area [ ]
  real, intent(in) :: FGS   (ILG)      !< Fractional coverage of snow over bare ground on modelled area [ ]
  real, intent(in) :: FC    (ILG)      !< Fractional coverage of vegetation over bare ground on modelled area [ ]
  real, intent(in) :: FG    (ILG)      !< Fractional coverage of bare ground on modelled area [ ]
  !
  ! =================================================================
  !
  !>
  !! In this subroutine, checks are carried out to ensure that the
  !! change in energy storage in each of the components of the
  !! modelled area (canopy, snow and soil) is equal to the sum of the
  !! energy fluxes into and out of them; and that the change in
  !! moisture storage in each of the components is equal to the sum of
  !! the water fluxes into and out of them. The subroutine is called
  !! twice, once at the beginning (ISTEP=0) and once at the end
  !! (ISTEP=1) of each time step. At the beginning, the instantaneous
  !! energy and moisture storage terms are evaluated, and at the end
  !! the differences over the time step are calculated:
  !!
  !! Change in canopy energy storage =
  !!     \f$\Delta [(c_c W_c + c_w W_{l,c} + c_i W_{f,c} )T_c ] / \Delta t\f$
  !!
  !! Change in snow energy storage =
  !!     \f$\Delta [(C_i W_s /\rho_i + C_w W_{l,s} / \rho_w)T_s ]/ \Delta t \f$
  !!
  !! Change in soil layer energy storage =
  !!     \f$\Delta {[(C_w \theta_l + C_i \theta_f + C_g \theta_g) \Delta z_w + C_b (\Delta z – \Delta z_w)]T_j }/ \Delta t\f$
  !! (For the first soil layer, the numerator contains the additional
  !! term \f$C_w z_p T_p\f$.)
  !!
  !! Change in canopy moisture storage = \f$\Delta[W_{l,c} + W_{f,c}]\f$
  !!
  !! Change in snow moisture storage = \f$\Delta[W_s + W_{l,s}]\f$
  !!
  !! Change in soil moisture storage =
  !!     \f$\Delta [(\theta_l \rho_w + \theta_f \rho_i) \Delta z_w + z_p \rho_w ]\f$
  !!
  !! The net energy and moisture fluxes are also evaluated at the end
  !! of the time step:
  !!
  !! Net energy flux for canopy =
  !!     \f$K_{*,c} + L_{*,c} – Q_{H,c} – Q_{E,c} – Q_{M,c} + Q_{I,c}\f$
  !!
  !! Net energy flux for snow =
  !!     \f$K_{*,s} + L_{*,s} – Q_{H,s} – Q_{E,s} – Q_{M,s} + Q_{I,s}\f$
  !!
  !! Net energy flux for first soil layer =
  !!     \f$K_{*,g} + L_{*,g} – Q_{H,g} – Q_{E,g} – Q_{M,1} + Q_{I,1}\f$
  !!
  !! Net energy flux for other soil layers = \f$- Q_{M,j} + Q_{I,j}\f$
  !!
  !! Net moisture flux for canopy =
  !!     \f$P_{l,c} + P_{f,c} – E_{l,c} – E_{f,c} – R_c + A_c\f$
  !!
  !! Net moisture flux for snow = \f$P_s – E_s – R_s + A_s\f$
  !!
  !! Net moisture flux for soil = \f$P_g – E_g – R_g + A_g - E_c\f$
  !!
  !! In these equations the \f$K_*\f$ terms refer to net shortwave radiation,
  !! the \f$L_*\f$ terms to net longwave radiation, the \f$Q_H\f$ terms to sensible
  !! heat flux, the \f$Q_E\f$ terms to latent heat flux, the \f$Q_M\f$ terms to heat
  !! associated with melting or freezing of water, and the \f$Q_I\f$ terms to
  !! changes in heat storage caused by conduction or redistribution of
  !! water. The P terms refer to precipitation, the E terms to
  !! evaporation, the R terms to runoff and the A terms to water
  !! transferred between different components of the landscape.
  !! The subscript 1 refers to the first soil layer, and j to a
  !! generalized other layer.
  !!
  !! Finally, each change in energy or moisture storage is compared in
  !! turn with the corresponding net flux of energy or moisture, and
  !! if the difference is greater than a selected threshold value, an
  !! error message is printed out and the run is stopped.
  !!

  if (ISTEP == 0) then
    !
    !     * SET BALANCE CHECK VARIABLES FOR START OF CURRENT TIME STEP.
    !
    do I = IL1,IL2
      WTGSTP(I) = 0.0
      CTVSTP(I) = - (CMAI(I) * SPHVEG + RCAN(I) * SPHW + &
                  SCAN(I) * SPHICE) * TCAN(I)
      CTSSTP(I) = - TSNOW(I) * (HCPICE * SNO(I) / RHOICE + &
                  HCPW * WSNOW(I) / RHOW)
      CT1STP(I) = - ((HCPW * THLIQ(I,1) + HCPICE * THICE(I,1) &
                  + HCPS(I,1) * (1.0 - THPOR(I,1))) * DELZW(I,1) + &
                  HCPSND * (DELZ(1) - DELZW(I,1))) * TBAR(I,1) - &
                  HCPW * ZPOND(I) * TPOND(I)
      CT2STP(I) = - ((HCPW * THLIQ(I,2) + HCPICE * THICE(I,2) &
                  + HCPS(I,2) * (1.0 - THPOR(I,2))) * DELZW(I,2) + &
                  HCPSND * (DELZ(2) - DELZW(I,2))) * TBAR(I,2)
      CT3STP(I) = - ((HCPW * THLIQ(I,3) + HCPICE * THICE(I,3) &
                  + HCPS(I,3) * (1.0 - THPOR(I,3))) * DELZW(I,3) + &
                  HCPSND * (DELZ(3) - DELZW(I,3))) * TBAR(I,3)
      WTVSTP(I) = - (RCAN(I) + SCAN(I))
      WTSSTP(I) = - SNO(I) - WSNOW(I)
      do J = 1,IG
        WTGSTP(I) = WTGSTP(I) - &
                    (THLIQ(I,J) * RHOW + THICE(I,J) * RHOICE) * &
                    DELZW(I,J)
      end do ! loop 50
      WTGSTP(I) = WTGSTP(I) - ZPOND(I) * RHOW
    end do ! loop 100
    !
  end if
  !
  if (ISTEP == 1) then
    !
    !     * CHECK ENERGY AND WATER BALANCES OVER THE CURRENT TIME STEP.
    !
    do I = IL1,IL2 ! loop 200
      CTVSTP(I) = CTVSTP(I) + (CMAI(I) * SPHVEG + RCAN(I) * SPHW + &
                  SCAN(I) * SPHICE) * TCAN(I)
      CTSSTP(I) = CTSSTP(I) + TSNOW(I) * (HCPICE * SNO(I) / RHOICE + &
                  HCPW * WSNOW(I) / RHOW)
      CT1STP(I) = CT1STP(I) + ((HCPW * THLIQ(I,1) + HCPICE * THICE(I,1) &
                  + HCPS(I,1) * (1.0 - THPOR(I,1))) * DELZW(I,1) + &
                  HCPSND * (DELZ(1) - DELZW(I,1))) * TBAR(I,1) + &
                  HCPW * ZPOND(I) * TPOND(I)
      CT2STP(I) = CT2STP(I) + ((HCPW * THLIQ(I,2) + HCPICE * THICE(I,2) &
                  + HCPS(I,2) * (1.0 - THPOR(I,2))) * DELZW(I,2) + &
                  HCPSND * (DELZ(2) - DELZW(I,2))) * TBAR(I,2)
      CT3STP(I) = CT3STP(I) + ((HCPW * THLIQ(I,3) + HCPICE * THICE(I,3) &
                  + HCPS(I,3) * (1.0 - THPOR(I,3))) * DELZW(I,3) + &
                  HCPSND * (DELZ(3) - DELZW(I,3))) * TBAR(I,3)
      CTVSTP(I) = CTVSTP(I) / DELT
      CTSSTP(I) = CTSSTP(I) / DELT
      CT1STP(I) = CT1STP(I) / DELT
      CT2STP(I) = CT2STP(I) / DELT
      CT3STP(I) = CT3STP(I) / DELT
      WTVSTP(I) = WTVSTP(I) + RCAN(I) + SCAN(I)
      WTSSTP(I) = WTSSTP(I) + SNO(I) + WSNOW(I)
      do J = 1,IG ! loop 150
        WTGSTP(I) = WTGSTP(I) + &
                    (THLIQ(I,J) * RHOW + THICE(I,J) * RHOICE) * &
                    DELZW(I,J)
      end do ! loop 150
      WTGSTP(I) = WTGSTP(I) + ZPOND(I) * RHOW
    end do ! loop 200
    !
    do I = IL1,IL2
      QSUMV = FSGV(I) + FLGV(I) - HFSC(I) - HEVC(I) - &
              HMFC(I) + HTCC(I)
      QSUMS = FSGS(I) + FLGS(I) - HFSS(I) - HEVS(I) - &
              HMFN(I) + HTCS(I)
      QSUM1 = FSGG(I) + FLGG(I) - HFSG(I) - HEVG(I) - &
              HMFG(I,1) + HTC(I,1)
      QSUM2 = - HMFG(I,2) + HTC(I,2)
      QSUM3 = - HMFG(I,3) + HTC(I,3)
      WSUMV = (PCFC(I) + PCLC(I) - &
              QFCF(I) - QFCL(I) - ROFC(I) + &
              WTRC(I)) * DELT
      WSUMS = (PCPN(I) - QFN(I) - &
              ROFN(I) + WTRS(I)) * DELT
      WSUMG = (PCPG(I) - QFG(I) - &
              ROF(I) + WTRG(I)) * DELT
      do J = 1,IG
        WSUMG = WSUMG - QFC(I,J) * DELT
      end do ! loop 250
      !
      if (ABS(CTVSTP(I) - QSUMV) > 1.0) then
        write(6,6441) N,CTVSTP(I),QSUMV
6441    format(2X,'CANOPY ENERGY BALANCE  ',I8,2F20.8)
        write(6,6450) FSGV(I),FLGV(I),HFSC(I), &
              HEVC(I),HMFC(I),HTCC(I)
        write(6,6450) RCAN(I),SCAN(I),TCAN(I)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(CTSSTP(I) - QSUMS) > 7.0) then
        write(6,6442) N,I,CTSSTP(I),QSUMS
6442    format(2X,'SNOW ENERGY BALANCE  ',2I8,2F20.8)
        write(6,6450) FSGS(I),FLGS(I),HFSS(I), &
             HEVS(I),HMFN(I),HTCS(I)
        write(6,6450) TSNOW(I),SNO(I),WSNOW(I)
        write(6,6451) FCS(I),FGS(I),FC(I),FG(I)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(CT1STP(I) - QSUM1) > 5.0) then
        write(6,6443) N,I,CT1STP(I),QSUM1
        write(6,6450) FSGG(I),FLGG(I),HFSG(I), &
             HEVG(I),HMFG(I,1),HTC(I,1)
        write(6,6450) FSGS(I),FLGS(I),HFSS(I), &
             HEVS(I),HMFN(I),HTCS(I)
        write(6,6450) THLIQ(I,1) * RHOW * DELZW(I,1), &
             THLIQ(I,2) * RHOW * DELZW(I,2), &
             THLIQ(I,3) * RHOW * DELZW(I,3), &
             THICE(I,1) * RHOICE * DELZW(I,1), &
             THICE(I,2) * RHOICE * DELZW(I,2), &
             THICE(I,3) * RHOICE * DELZW(I,3), &
             ZPOND(I) * RHOW
        write(6,6451) FCS(I),FGS(I),FC(I),FG(I), &
             DELZW(I,1),DELZW(I,2),DELZW(I,3)
6443    format(2X,'LAYER 1 ENERGY BALANCE  ',2I8,2F20.8)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(CT2STP(I) - QSUM2) > 5.0) then
        write(6,6444) N,I,CT2STP(I),QSUM2
6444    format(2X,'LAYER 2 ENERGY BALANCE  ',2I8,2F20.8)
        write(6,6450) HMFG(I,2),HTC(I,2), &
             THLIQ(I,2),THICE(I,2),THPOR(I,2),TBAR(I,2) - TFREZ
        write(6,6450) HMFG(I,3),HTC(I,3), &
             THLIQ(I,3),THICE(I,3),THPOR(I,3),TBAR(I,3) - TFREZ
        write(6,6450) HMFG(I,1),HTC(I,1), &
             THLIQ(I,1),THICE(I,1),THPOR(I,1),TBAR(I,1) - TFREZ
        write(6,6451) FCS(I),FGS(I),FC(I),FG(I), &
             DELZW(I,2),HCPS(I,2),DELZW(I,3)
6451    format(2X,7E20.6)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(CT3STP(I) - QSUM3) > 10.0) then
        write(6,6445) N,I,CT3STP(I),QSUM3
6445    format(2X,'LAYER 3 ENERGY BALANCE  ',2I8,2F20.8)
        write(6,6450) HMFG(I,3),HTC(I,3), &
             TBAR(I,3)
        write(6,6450) THLIQ(I,3),THICE(I,3),HCPS(I,3), &
                       THPOR(I,3),DELZW(I,3)
        ! CALL EXIT
      end if
      if (ABS(WTVSTP(I) - WSUMV) > 1.0E-3) then
        write(6,6446) N,WTVSTP(I),WSUMV
6446    format(2X,'CANOPY WATER BALANCE  ',I8,2F20.8)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(WTSSTP(I) - WSUMS) > 1.0E-2) then
        write(6,6447) N,I,WTSSTP(I),WSUMS
6447    format(2X,'SNOW WATER BALANCE  ',2I8,2F20.8)
        write(6,6450) PCPN(I) * DELT,QFN(I) * DELT, &
             ROFN(I) * DELT,WTRS(I) * DELT
        write(6,6450) SNO(I),WSNOW(I),TSNOW(I) - TFREZ
        write(6,6451) FCS(I),FGS(I),FC(I),FG(I)
        ! CALL EXIT
        call abandonCell
      end if
      if (ABS(WTGSTP(I) - WSUMG) > 1.0E-1) then
        write(6,6448) N,I,WTGSTP(I),WSUMG
6448    format(2X,'GROUND WATER BALANCE  ',2I8,2F20.8)
        write(6,6450) PCPG(I) * DELT,QFG(I) * DELT, &
             QFC(I,1) * DELT,QFC(I,2) * DELT, &
             QFC(I,3) * DELT,ROF(I) * DELT, &
             WTRG(I) * DELT
        do J = 1,IG
          write(6,6450) THLIQ(I,J) * RHOW * DELZW(I,J), &
                 THICE(I,J) * RHOICE * DELZW(I,J), &
                 DELZW(I,J)
        end do ! loop 390
        write(6,6450) ZPOND(I) * RHOW
6450    format(2X,7F15.6)
        write(6,6451) FCS(I),FGS(I),FC(I),FG(I)
        ! CALL EXIT
        call abandonCell
      end if
    end do ! loop 400

  end if

  return
end subroutine energyWaterBalanceCheck
