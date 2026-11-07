!> \file
!! Checks for freezing or thawing of liquid or frozen water
!! on the vegetation canopy, and adjust canopy temperature and
!! intercepted water stores accordingly.
!! @author D. Verseghy, M. Lazare
!
subroutine canopyPhaseChange (TCAN, RAICAN, SNOCAN, FRAINC, FSNOWC, CHCAP, & ! Formerly CWCALC
                              HMFC, HTCC, FI, CMASS, ILG, IL1, IL2, JL)
  !
  !
  !     * MAR 25/08 - D.VERSEGHY. UPDATE FRAINC AND FSNOWC.
  !     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUN 20/02 - D.VERSEGHY. COSMETIC REARRANGEMENT OF
  !     *                         SUBROUTINE CALL; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * MAR 17/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 13/91 - D.VERSEGHY. ADJUST CANOPY TEMPERATURE AND
  !     *                         INTERCEPTED LIQUID/FROZEN
  !     *                         MOISTURE STORES FOR FREEZING/
  !     *                         THAWING.
  !
  use classicParams,       only : DELT, TFREZ, SPHW, SPHICE, SPHVEG, CLHMLT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL
  integer :: I
  !
  !     * I/O ARRAYS.
  !
  real, intent(inout) :: TCAN  (ILG)  !< Temperature of vegetation canopy \f$[K] (T_c)\f$
  real, intent(inout) :: RAICAN(ILG)  !< Intercepted liquid water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: SNOCAN(ILG)  !< Intercepted frozen water stored on the canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: FRAINC(ILG)  !< Fractional coverage of canopy by liquid water [ ]
  real, intent(inout) :: FSNOWC(ILG)  !< Fractional coverage of canopy by frozen water [ ]
  real, intent(inout) :: CHCAP (ILG)  !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1}] (C_c)\f$
  real, intent(inout) :: HMFC  (ILG)  !< Energy associated with freezing or thawing of water in canopy interception stores \f$[W m^{-2}]\f$
  real, intent(inout) :: HTCC  (ILG)  !< Internal energy change of canopy due to changes in temperature and/or mass \f$[W m^{-2}] (I_c)\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI   (ILG)  !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: CMASS(ILG)  !< Mass of vegetation canopy \f$[kg m^{-2}]\f$

  !
  !     * TEMPORARY VARIABLES.
  !
  real :: HFREZ, HCONV, RCONV, HCOOL, HMELT, SCONV, HWARM
  !
  !---------------------------------------------------------------------
  !>
  !! The change of internal energy \f$I_c\f$ of the vegetation canopy as a
  !! result of the phase change processes treated here is calculated as
  !! the difference in \f$I_c\f$ between the beginning and end of the subroutine:
  !!
  !! \f$\Delta I_c = X_i \Delta [C_c T_c ] / \Delta t\f$
  !!
  !! where \f$C_c\f$ represents the volumetric heat capacity of the canopy, \f$T_c\f$
  !! its temperature, \f$\Delta t\f$ the length of the time step, and \f$X_i\f$ the
  !! fractional coverage of the subarea under consideration relative to
  !! the modelled area.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0.) then
      HTCC  (I) = HTCC(I) - FI(I) * TCAN(I) * CHCAP(I) / DELT
      !>
      !! If there is liquid water stored on the canopy and the
      !! canopy temperature is less than 0 C, the available
      !! energy sink HFREZ is calculated from CHCAP and the
      !! difference between TCAN and 0 C, and
      !! compared with HCONV, calculated as the energy sink
      !! required to freeze all of the liquid water on the
      !! canopy. If HFREZ \f$\leq\f$ HCONV, the amount of water that
      !! can be frozen is calculated using the latent heat
      !! of melting. The fractional coverages of frozen and
      !! liquid water FSNOWC and FRAINC and their masses
      !! SNOCAN and RAICAN are adjusted accordingly, TCAN is
      !! set to 0 C, and the amount of energy
      !! involved is subtracted from the internal energy HTCC
      !! and added to HMFC. Otherwise all of the
      !! intercepted liquid water is converted to frozen
      !! water, and the energy available for cooling the canopy is
      !! calculated as HCOOL = HFREZ – HCONV. This available
      !! energy is applied to decreasing the
      !! temperature of the canopy, using the specific heat of
      !! the canopy elements, and the amount of energy that
      !! was involved in the phase change is subtracted from
      !! HTCC and added to HMFC.
      !!
      if (RAICAN(I) > 0. .and. TCAN(I) < TFREZ) then
        HFREZ = CHCAP(I) * (TFREZ - TCAN(I))
        HCONV = RAICAN(I) * CLHMLT
        if (HFREZ <= HCONV) then
          RCONV = HFREZ / CLHMLT
          FSNOWC(I) = FSNOWC(I) + FRAINC(I) * RCONV / RAICAN(I)
          FRAINC(I) = FRAINC(I) - FRAINC(I) * RCONV / RAICAN(I)
          SNOCAN(I) = SNOCAN(I) + RCONV
          RAICAN(I) = RAICAN(I) - RCONV
          TCAN  (I) = TFREZ
          HMFC  (I) = HMFC(I) - FI(I) * CLHMLT * RCONV / DELT
          HTCC  (I) = HTCC(I) - FI(I) * CLHMLT * RCONV / DELT
        else
          HCOOL = HFREZ - HCONV
          SNOCAN(I) = SNOCAN(I) + RAICAN(I)
          FSNOWC(I) = FSNOWC(I) + FRAINC(I)
          FRAINC(I) = 0.0
          TCAN  (I) = - HCOOL / (SPHVEG * CMASS(I) + SPHICE * &
                      SNOCAN(I)) + TFREZ
          HMFC  (I) = HMFC(I) - FI(I) * CLHMLT * RAICAN(I) / DELT
          HTCC  (I) = HTCC(I) - FI(I) * CLHMLT * RAICAN(I) / DELT
          RAICAN(I) = 0.0
        end if
      end if
      !>
      !! If there is frozen water stored on the canopy and the
      !! canopy temperature is greater than 0 C, the available
      !! energy for melting, HMELT, is calculated from CHCAP and
      !! the difference between TCAN and 0 C, and
      !! compared with HCONV, calculated as the energy required to
      !! melt all of the frozen water on the canopy.
      !! If HMELT \f$\leq\f$ HCONV, the amount of frozen water that can be
      !! melted is calculated using the latent heat
      !! of melting. The fractional coverages of frozen and liquid
      !! water FSNOWC and FRAINC and their masses
      !! SNOCAN and RAICAN are adjusted accordingly, TCAN is set
      !! to 0 C, and the amount of energy
      !! involved is subtracted from HTCC and added to HMFC.
      !! Otherwise, all of the intercepted frozen water is
      !! converted to liquid water, and the energy available for
      !! warming the canopy is calculated as HWARM =
      !! HMELT – HCONV. This available energy is applied to
      !! increasing the temperature of the canopy, using
      !! the specific heats of the canopy elements, and the amount
      !! of energy that was involved in the phase
      !! change is subtracted from HTCC and added to HMFC.
      !!
      if (SNOCAN(I) > 0. .and. TCAN(I) > TFREZ) then
        HMELT = CHCAP(I) * (TCAN(I) - TFREZ)
        HCONV = SNOCAN(I) * CLHMLT
        if (HMELT <= HCONV) then
          SCONV = HMELT / CLHMLT
          FRAINC(I) = FRAINC(I) + FSNOWC(I) * SCONV / SNOCAN(I)
          FSNOWC(I) = FSNOWC(I) - FSNOWC(I) * SCONV / SNOCAN(I)
          SNOCAN(I) = SNOCAN(I) - SCONV
          RAICAN(I) = RAICAN(I) + SCONV
          TCAN(I) = TFREZ
          HMFC  (I) = HMFC(I) + FI(I) * CLHMLT * SCONV / DELT
          HTCC  (I) = HTCC(I) + FI(I) * CLHMLT * SCONV / DELT
        else
          HWARM = HMELT - HCONV
          RAICAN(I) = RAICAN(I) + SNOCAN(I)
          FRAINC(I) = FRAINC(I) + FSNOWC(I)
          FSNOWC(I) = 0.0
          TCAN(I) = HWARM / (SPHVEG * CMASS(I) + SPHW * RAICAN(I)) + &
                    TFREZ
          HMFC  (I) = HMFC(I) + FI(I) * CLHMLT * SNOCAN(I) / DELT
          HTCC  (I) = HTCC(I) + FI(I) * CLHMLT * SNOCAN(I) / DELT
          SNOCAN(I) = 0.0
        end if
      end if
      !>
      !! In the final cleanup, the canopy heat capacity is
      !! recomputed and the remaining internal energy calculations
      !! are completed.
      !!
      CHCAP(I) = SPHVEG * CMASS(I) + SPHW * RAICAN(I) + SPHICE * SNOCAN(I)
      HTCC (I) = HTCC(I) + FI(I) * TCAN(I) * CHCAP(I) / DELT
    end if
  end do ! loop 100
  return
end subroutine canopyPhaseChange
