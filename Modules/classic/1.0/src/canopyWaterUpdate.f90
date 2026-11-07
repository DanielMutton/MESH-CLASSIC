!> \file
!! Updates liquid and frozen water stores on canopy and in soil in
!! response to calculated sublimation, evaporation and transpiration
!! rates.
!! @author D. Verseghy, M. Lazare
!
subroutine canopyWaterUpdate (EVAP, SUBL, RAICAN, SNOCAN, TCAN, THLIQ, TBAR, ZSNOW, & ! Formerly CANVAP
                              WLOST, CHCAP, QFCF, QFCL, QFN, QFC, HTCC, HTCS, HTC, &
                              FI, CMASS, TSNOW, HCPSNO, RHOSNO, FROOT, THPOR, &
                              THLMIN, DELZW, EVLOST, RLOST, IROOT, &
                              IG, ILG, IL1, IL2, JL, N)

  !     * SEP 15/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3.
  !     * SEP 13/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUN 20/02 - D.VERSEGHY. TIDY UP SUBROUTINE CALL; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE SOIL
  !     *                         PERMEABLE DEPTH.
  !     * DEC 30/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         BUGFIXES IN CALCULATION OF QFN AND
  !     *                         QFC.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         RATIONALIZE CALCULATION OF WLOST
  !     *                         REFINE CALCULATION OF QFCL.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         ADDITIONAL DIAGNOSTIC CALCULATIONS -
  !     *                         HTCC AND HTC.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !                                        NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CALCULATE ACTUAL EVAPORATION,
  !     *                         SUBLIMATION AND TRANSPIRATION FROM
  !     *                         VEGETATION CANOPY.
  !
  use classicParams,       only : DELT, TFREZ, HCPW, SPHW, SPHICE, SPHVEG, &
                                  RHOW, CLHMLT, CLHVAP


  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, ILG, IL1, IL2, JL, N
  integer             :: I, J
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  ! layer \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(inout) :: TBAR  (ILG,IG)   !< Temperature of soil layer \f$[K] (T_g)\f$
  real, intent(inout) :: QFC   (ILG,IG)   !< Transpired water removed from soil layer \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: HTC   (ILG,IG)   !< Internal energy change of soil layer due to
  !! conduction and/or change in mass \f$[W m^{-2}] (I_g)\f$
  !
  real, intent(inout) :: EVAP  (ILG)  !< Evapotranspiration rate from vegetation canopy \f$[m s^{-1}]\f$
  real, intent(inout) :: SUBL  (ILG)  !< Calculated sublimation rate from vegetation canopy \f$[m s^{-1}]\f$
  real, intent(inout) :: RAICAN(ILG)  !< Intercepted liquid water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: SNOCAN(ILG)  !< Intercepted frozen water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: TCAN  (ILG)  !< Temperature of vegetation canopy \f$[K] (T_c)\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_g)\f$
  real, intent(inout) :: WLOST (ILG)  !< Residual amount of water that cannot be
  !! supplied by surface stores \f$[kg m^{-2}]\f$
  real, intent(inout) :: CHCAP (ILG)  !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1}] (C_c)\f$
  real, intent(inout) :: QFCF  (ILG)  !< Sublimation from frozen water in canopy
  !! interception store \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFCL  (ILG)  !< Evaporation from liquid water in canopy
  !! interception store \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(inout) :: QFN   (ILG)  !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: HTCC  (ILG)  !< Internal energy change of canopy due to changes
  !! in temperature and/or mass \f$[W m^{-2}] (I_c)\f$
  real, intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to
  !! conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FROOT (ILG,IG)   !< Fractional contribution of soil layer to
  !! transpiration [ ]
  real, intent(in) :: THPOR(ILG,IG)    !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !! remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer \f$[m] (\Delta z_{g,w})\f$
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !! modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: CMASS (ILG)  !< Mass of vegetation canopy \f$[kg m^{-2}]\f$
  real, intent(in) :: TSNOW (ILG)  !< Temperature of the snow pack [C]
  real, intent(in) :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(in) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}]\f$
  !
  !     * WORK ARRAYS.
  !
  real, intent(inout)    :: EVLOST(ILG), RLOST(ILG)
  !
  integer, intent(inout) :: IROOT(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: SLOST, THTRAN, THLLIM

  ! C-----------------------------------------------------------------------
  !>
  !! The calculated fluxes of liquid and frozen water from the canopy
  !! to the overlying air, obtained as outputs of waterCalcPrep, are applied
  !! to the liquid and frozen intercepted water stores on the
  !! vegetation canopy, and to the liquid and frozen moisture stores
  !! at the surface and in the soil. Since there may not be sufficient
  !! water in one or more of these stores to sustain the calculated
  !! rate over the whole time step, a hierarchy of operations is
  !! followed as described below. The change of internal energy HTC in
  !! the canopy, snow and soil layers as a result of these processes
  !! is calculated as the difference in HTC between the beginning and
  !! end of the subroutine:
  !! \f[
  !! \Delta I_c = X_i \Delta (C_c T_c)/ \Delta t
  !! \f]\f[
  !! \Delta I_s = X_i \Delta (C_s T_s z_s)/ \Delta t
  !! \f]
  !! where the C terms represent volumetric heat capacities and the T
  !! terms temperatures of the canopy and snow pack, \f$\Delta t\f$ is the
  !! length of the time step, \f$z_s\f$ the snow depth, and \f$X_i\f$ the fractional
  !! coverage of the subarea under consideration relative to the
  !! modelled area. For the soil layers, since only the liquid water
  !! content is affected by these calculations, the change in internal
  !! energy of each layer is calculated from the change in liquid
  !! water content \f$\theta_l\f$ as:
  !! \f[ \Delta I_g = X_i C_w \Delta z_{g, w} \Delta (T_g \theta_l)/ \Delta t \f]
  !!


  !     * INITIALIZE ARRAYS.
  !     * (THE WORK ARRAY "IROOT" INDICATES POINTS WHERE TRANSPIRATION
  !     * CAN OCCUR.)
  !
  do I = IL1,IL2 ! loop 50
    if (FI(I) > 0.) then
      RLOST (I) = 0.0
      EVLOST(I) = 0.0
      IROOT (I) = 0
      HTCC  (I) = HTCC(I) - FI(I) * TCAN(I) * CHCAP(I) / DELT
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
    end if
  end do ! loop 50
  !
  do J = 1,IG ! loop 100
    do I = IL1,IL2
      if (FI(I) > 0.) then
        HTC (I,J) = HTC(I,J) - FI(I) * (TBAR(I,J) + TFREZ) * THLIQ(I,J) * &
                    HCPW * DELZW(I,J) / DELT
        if (FROOT(I,J) > 1.0E-5) IROOT(I) = 1
      end if
    end do
  end do ! loop 100
  !>
  !! Sublimation is addressed first. The predicted mass of sublimated
  !! water SLOST is calculated and compared to the frozen water in the
  !! canopy interception store, SNOCAN. If SLOST \f$\leq\f$ SNOCAN, all of the
  !! sublimated water is subtracted from SNOCAN. If not, the excess
  !! sublimation is calculated as SLOST-SNOCAN, QFCF is corrected for
  !! the canopy sublimation difference, and SNOCAN is set to zero.
  !! Next, the new value of SLOST is compared to the snowpack mass,
  !! calculated as ZSNOW*RHOSNO. If SLOST \f$\leq\f$ ZSNOW*RHOSNO, all of the
  !! remaining sublimated water is taken from the snow pack, and QFN
  !! is modified to reflect this loss. Otherwise, the excess
  !! sublimation is calculated as SLOST - ZSNOW*RHOSNO, QFN is
  !! adjusted accordingly, and ZSNOW is set to zero. There now remain
  !! no further frozen moisture stores from which sublimated water can
  !! be taken (frozen water in the soil is assumed to be immobile), so
  !! the remaining energy that had been assigned to sublimation is
  !! assigned to canopy evaporation instead, and QFCL is duly
  !! recalculated. This means, however, that a small imbalance will
  !! arise in the water budget owing to the difference between the
  !! latent heats of sublimation and evaporation. This imbalance is
  !! assigned to the housekeeping variable WLOST.
  !
  !     * SUBLIMATION CASE.  IF SNOW ON CANOPY IS INSUFFICIENT TO SUPPLY
  !     * DEMAND, RESIDUAL IS TAKEN FIRST FROM SNOW UNDERLYING CANOPY AND
  !     * THEN FROM LIQUID WATER ON CANOPY.
  !
  do I = IL1,IL2 ! loop 200
    if (FI(I) > 0. .and. SUBL(I) > 0.) then
      SLOST = SUBL(I) * DELT * RHOW
      if (SLOST <= SNOCAN(I)) then
        SNOCAN(I) = SNOCAN(I) - SLOST
        SUBL(I) = 0.0
      else
        SLOST = SLOST - SNOCAN(I)
        QFCF(I) = QFCF(I) - FI(I) * SLOST / DELT
        SNOCAN(I) = 0.0
        if (SLOST <= ZSNOW(I) * RHOSNO(I)) then
          ZSNOW(I) = ZSNOW(I) - SLOST / RHOSNO(I)
          SUBL(I) = 0.0
          QFN(I) = QFN(I) + FI(I) * SLOST / DELT
        else
          SLOST = SLOST - ZSNOW(I) * RHOSNO(I)
          QFN(I) = QFN(I) + FI(I) * ZSNOW(I) * RHOSNO(I) / DELT
          ZSNOW(I) = 0.0
          WLOST(I) = WLOST(I) - SLOST * CLHMLT / CLHVAP
          EVAP(I) = EVAP(I) + SLOST * (CLHMLT + CLHVAP) / &
                    (CLHVAP * DELT * RHOW)
          QFCL(I) = QFCL(I) + FI(I) * SLOST * (CLHMLT + CLHVAP) / &
                    (CLHVAP * DELT)
        end if
      end if
    end if
  end do ! loop 200
  !>
  !! Now canopy evaporation is addressed. It is assumed that all
  !! intercepted liquid water evaporates before transpiration begins,
  !! since there is a canopy stomatal resistance associated with
  !! transpiration and there is none associated with evaporation. The
  !! predicted mass of evaporated water RLOST is calculated and
  !! compared to the liquid water in the canopy interception store,
  !! RAICAN. If RLOST \f$\leq\f$ RAICAN, all of the evaporated water is
  !! subtracted from RAICAN. If not, the excess evaporation is
  !! calculated as RLOST - RAICAN, and QFCL is corrected for the
  !! canopy evaporation difference. This excess evaporation is now
  !! treated as transpiration. An initial check is done by referring
  !! to the diagnostic flag IROOT, which was set to 1 at the beginning
  !! of the subroutine if there was water available for transpiration
  !! in any of the soil layers. If IROOT is zero, no transpiration can
  !! occur and the excess evaporation is stored in the temporary
  !! variable EVLOST.
  !
  !     * EVAPORATION.  IF WATER ON CANOPY IS INSUFFICIENT TO SUPPLY
  !     * DEMAND, ASSIGN RESIDUAL TO TRANSPIRATION.
  !
  do I = IL1,IL2 ! loop 300
    if (FI(I) > 0. .and. EVAP(I) > 0.) then
      RLOST(I) = EVAP(I) * RHOW * DELT
      if (RLOST(I) <= RAICAN(I)) then
        RAICAN(I) = RAICAN(I) - RLOST(I)
        EVAP  (I) = 0.
        RLOST (I) = 0.
      else
        RLOST(I) = RLOST(I) - RAICAN(I)
        QFCL(I) = QFCL(I) - FI(I) * RLOST(I) / DELT
        if (IROOT(I) == 0) EVLOST(I) = RLOST(I)
        EVAP  (I) = 0.
        RAICAN(I) = 0.
      end if
    end if
  end do ! loop 300
  !>
  !! The next loop is performed if IROOT=1, i.e. if transpiration is
  !! possible. For each soil layer, the volumetric water content that
  !! is removed by transpiration, THTRAN, is calculated from RLOST
  !! (converted to a volumetric content by dividing by the density of
  !! water and the permeable thickness of the soil layer), and the
  !! fractional contribution of the soil layer FROOT. If there is
  !! enough liquid water in the soil layer to supply THTRAN, the
  !! diagnostic transpiration flux QFC for the layer is updated using
  !! THTRAN, and the liquid water content of the layer is updated as
  !! THLIQ-THTRAN. If not, QFC is updated using the available water in
  !! the soil layer, THLIQ is set to THLMIN, and the residual
  !! untranspired water is added to EVLOST.
  !
  !     * TRANSPIRATION.
  !
  do J = 1,IG ! loop 400
    do I = IL1,IL2
      if (FI(I) > 0. .and. IROOT(I) > 0) then
        if (DELZW(I,J) > 0.0) then
          THTRAN = RLOST(I) * FROOT(I,J) / (RHOW * DELZW(I,J))
        else
          THTRAN = 0.0
        end if
        if (THPOR(I,J) < THLMIN(I,J)) then
          THLLIM = THPOR(I,J)
        else
          THLLIM = THLMIN(I,J)
        end if
        if (THTRAN <= (THLIQ(I,J) - THLLIM)) then
          QFC  (I,J) = QFC(I,J) + FI(I) * RLOST(I) * FROOT(I,J) / DELT
          THLIQ(I,J) = THLIQ(I,J) - THTRAN
        else
          QFC  (I,J) = QFC(I,J) + FI(I) * (THLIQ(I,J) - THLLIM) * RHOW * &
                       DELZW(I,J) / DELT
          EVLOST (I) = EVLOST(I) + (THTRAN + THLLIM - THLIQ(I,J)) * RHOW * &
                       DELZW(I,J)
          THLIQ(I,J) = THLLIM
        end if
      end if
    end do
  end do ! loop 400
  !>
  !! In the final cleanup, the canopy heat capacity is recalculated,
  !! the contents of EVLOST are added to WLOST, and the remaining
  !! internal energy calculations are completed.
  !
  !     * CLEANUP.
  !
  do I = IL1,IL2 ! loop 500
    if (FI(I) > 0.) then
      CHCAP(I) = RAICAN(I) * SPHW + SNOCAN(I) * SPHICE + CMASS(I) * SPHVEG
      WLOST(I) = WLOST(I) + EVLOST(I)
      HTCC  (I) = HTCC(I) + FI(I) * TCAN(I) * CHCAP(I) / DELT
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
    end if
  end do ! loop 500
  !
  do J = 1,IG ! loop 550
    do I = IL1,IL2
      if (FI(I) > 0.) then
        HTC (I,J) = HTC(I,J) + FI(I) * (TBAR(I,J) + TFREZ) * THLIQ(I,J) * &
                    HCPW * DELZW(I,J) / DELT
      end if
    end do
  end do ! loop 550
  !
  return
end subroutine canopyWaterUpdate
