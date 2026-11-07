!> \file
!! Calculate decrease in snow albedo and increase in density due to aging.
!
subroutine snowAging (ALBSNO, RHOSNO, ZSNOW, HCPSNO, TSNOW, & ! Formerly SNOALBW
                      FI, S, RMELT, WSNOW, RHOMAX, ISAND, &
                      ILG, IG, IL1, IL2, JL)
  !
  !     * APR 17/14 - D.VERSEGHY. MAKE SNOW ALBEDO REFRESHMENT VALUE
  !     *                         CONSISTENT WITH snowAddNew.
  !     * MAR 07/07 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS.
  !     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * SEP 10/04 - R.HARVEY/D.VERSEGHY. INCREASE SNOW ALBEDO
  !     *                         REFRESHMENT AND WETTING THRESHOLDS.
  !     * AUG 04/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE
  !     *                         CALCULATIONS AGAIST ROUNDOFF
  !     *                         ERRORS.
  !     * APR 21/04 - F.SEGLENIEKS/D.VERSEGHY. BUG FIX IN SNOW
  !     *                         TEMPERATURE COMPARISONS.
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * OCT 20/00 - R.BROWN/D.VERSEGHY. MODIFIED SNOW DENSITY
  !     *                                 CALCULATIONS, ACCOUNTING
  !     *                                 FOR SETTLING IN WARM AND
  !     *                                 COLD SNOW.
  !     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         SPECIFY LOCATION OF ICE SHEETS
  !     *                         BY SOIL TEXTURE ARRAY RATHER
  !     *                         THAN BY SOIL COLOUR INDEX.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * MAR 13/92 - M.LAZARE.   CLASS - VERSION 2.1.
  !     *                         CODE FOR MODEL VERSION GCM7 -
  !     *                         DIVIDE PREVIOUS SUBROUTINE
  !     *                         "SNOALB" INTO "snowAlbedoTransmiss" AND
  !     *                         "snowAging" AND VECTORIZE.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. CALCULATE DECREASE IN SNOW ALBEDO
  !     *                         AND INCREASE IN DENSITY DUE TO
  !     *                         AGING. (ASSIGN DIFFERENT LOWER
  !     *                         SNOW ALBEDO LIMITS FOR DRY AND
  !     *                         MELTING SNOW.)
  !
  use classicParams, only : DELT, HCPW, HCPICE, RHOW, RHOICE

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IG, IL1, IL2, JL
  integer :: I, IPTBAD
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(inout) :: ALBSNO(ILG)  !< Albedo of snow \f$[ ] (\alpha_s)\f$
  real, intent(inout) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_s)\f$
  real, intent(out)   :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: TSNOW (ILG)  !< Temperature of the snow pack \f$[C]\f$
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: S     (ILG)  !< Snowfall rate \f$[m s^{-1}] \f$
  real, intent(in) :: RMELT (ILG)  !< Melt rate at top of snow pack \f$[m s^{-1}]\f$
  real, intent(in) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}]\f$
  !
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * WORK ARRAY.
  !
  real, intent(inout) :: RHOMAX(ILG)  !< Maximum density of snow pack \f$[kg m^{-3}] (\rho_{s, max})\f$
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: TIMFAC,RHOOLD
  !
  !----------------------------------------------------------------------

  IPTBAD = 0
  do I = IL1,IL2 ! loop 100
    if (ZSNOW(I) > 0. .and. &
        FI  (I) > 0. .and. S(I) * DELT < 1.0E-4) then
      if (ALBSNO(I) > 0.5001 .and. (RMELT(I) > 1.0E-7 .or. &
          TSNOW(I) >= - 0.01)) then
        ALBSNO(I) = (ALBSNO(I) - 0.50) * EXP( - 0.01 * DELT / 3600.0) + &
                    0.50
      else if (ALBSNO(I) > 0.7001 .and. RMELT(I) <= 1.0E-7) &
                                                    then
        ALBSNO(I) = (ALBSNO(I) - 0.70) * EXP( - 0.01 * DELT / 3600.0) + &
                    0.70
      end if
    end if
    !
    if (FI(I) > 0. .and. ZSNOW(I) > 0.0001) then
      if (TSNOW(I) < - 0.01) then
        RHOMAX(I) = 450.0 - (204.7 / ZSNOW(I)) * &
                    (1.0 - EXP( - ZSNOW(I) / 0.673))
      else
        RHOMAX(I) = 700.0 - (204.7 / ZSNOW(I)) * &
                    (1.0 - EXP( - ZSNOW(I) / 0.673))
      end if
    end if
    !
    if (FI(I) > 0. .and. ZSNOW(I) > 0.0001 .and. &
        RHOSNO(I) < (RHOMAX(I) - 0.01)) then
      RHOOLD = RHOSNO(I)
      RHOSNO(I) = (RHOSNO(I) - RHOMAX(I)) * EXP( - 0.01 * DELT / 3600.0) + &
                  RHOMAX(I)
      ZSNOW(I) = ZSNOW(I) * RHOOLD / RHOSNO(I)
      HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                  (RHOW * ZSNOW(I))
    end if
    if ((ALBSNO(I) < 0.49 .or. ALBSNO(I) > 1.0) .and. &
        ZSNOW (I) > 0. .and. FI(I) > 0.)               IPTBAD = I
  end do ! loop 100
  !
  if (IPTBAD /= 0) then
    write(6,6100) IPTBAD,JL,ALBSNO(IPTBAD)
6100 format('0AT (I,J) = (',I3,',',I3,'),ALBSNO = ',F10.5)
    call errorHandler('snowAging', - 1)
  end if
  !
  return
end subroutine snowAging
!> \file
!!
!! @author D. Verseghy, M. Lazare, R. Brown, F. Seglenieks, Y. Delange, R. Harvey
!!
!! The albedo and density of snow are modelled using empirical
!! exponential decay functions. In the absence
!! of snowfall exceeding the snow albedo refreshment threshold (set to 0.0001 m
!! here and in snowAddNew), the snow albedo \f$\alpha_s\f$ is assumed to decrease
!! exponentially with time from a fresh
!! snow value of 0.84 to a background old snow value \f$\alpha_{s, old}\f$ using an
!! expression based on data given in
!! Aguado (1985) \cite Aguado1985-fv, Robinson and Kukla (1984) and Dirmhirn and Eaton
!! (1975) \cite Dirmhirn1975-vx :
!!
!! \f$\alpha_s (t+1) = [\alpha_s (t) - \alpha_{s, old}] exp [-0.01 \Delta t / 3600] + \alpha_{s, old}\f$
!!
!! where \f$\Delta t\f$ is the length of the time step. If the melt rate RMELT
!! at the top of the snow pack is non-
!! negligible or if the temperature of the snow is close to 0 C,
!! \f$\alpha_{s, old}\f$ is assigned a value of 0.50; otherwise \f$\alpha_{s, old}\f$
!! is assigned a value of 0.70.
!!
!! The maximum snow density \f$\rho_{s, max}\f$ is estimated as a function of
!! snow depth \f$z_s\f$, after Tabler et al. (1990):
!!
!! \f$\rho_{s, max} = A_s - [204.70/ z_s] [1.0 - exp(-z_s /0.673)]\f$
!!
!! The empirical constant \f$A_s\f$ is assigned a value of 450.0 for cold
!! snow packs, and 700.0 for snow packs near
!! the melting point, following Brown et al. (2006) \cite Brown2006-ec.
!!
!! The density of snow \f$\rho_s\f$ increases exponentially with time from its
!! fresh snow value to the background old
!! snow density calculated above, according to an expression
!! analogous to that for albedo, derived from the
!! field measurements of Longley (1960) and Gold (1958) \cite Gold1958-ng:
!!
!! \f$\rho_s (t+1) = [\rho_s (t) - \rho_{s, max} ] exp [-0.01 \Delta t/3600] + \rho{s, max}\f$
!!
!! The snow depth and heat capacity are adjusted (see notes on
!! src/snowSublimation.f90), and a check is
!! performed with a call to abort if for unphysical albedo values
!! are encountered.
!!
