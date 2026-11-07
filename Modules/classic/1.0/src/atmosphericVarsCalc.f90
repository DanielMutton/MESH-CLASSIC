!> \file
!>
!! Evaluates atmospheric variables and rainfall/snowfall rates over modelled area.
!! @author D. Verseghy, M. Lazare, R. Brown, S. Fassnacht, P. Bartlett
!!

subroutine atmosphericVarsCalc (VPD, TADP, PADRY, RHOAIR, RHOSNI, RPCP, TRPCP, & ! Formerly CLASSI
                                SPCP, TSPCP, TA, QA, PCPR, RRATE, SRATE, &
                                PRESSG, IPCP, NL, IL1, IL2)

  !     * NOV 17/11 - M.LAZARE.   REMOVE CALCULATION OF PCPR
  !     *                         FOR IPCP=4 (REDUNDANT SINCE
  !     *                         PCPR MUST BE PASSED IN FOR
  !     *                         IF CONDITION ON LINE 100).
  !     * NOV 22/06 - P.BARTLETT. CALCULATE PCPR IF IPCP=4.
  !     * JUN 06/06 - V.FORTIN.   ADD OPTION FOR PASSING IN
  !     *                         RAINFALL AND SNOWFALL RATES
  !     *                         CALCULATED BY ATMOSPHERIC MODEL.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND
  !     *                         MOVE CLOUD FRACTION CALCULATION
  !     *                         BACK INTO DRIVER.
  !     * SEP 04/03 - D.VERSEGHY. NEW LOWER LIMIT ON PRECIPITATION
  !     *                         RATE.
  !     * AUG 09/02 - D.VERSEGHY. MOVE CALCULATION OF SOME
  !     *                         ATMOSPHERIC VARIABLES HERE
  !     *                         PRIOR TO GATHERING.
  !     * JUL 26/02 - R.BROWN/S.FASSNACHT/D.VERSEGHY. PROVIDE
  !     *                         ALTERNATE METHODS OF ESTIMATING
  !     *                         RAINFALL/SNOWFALL PARTITIONING.
  !     * JUN 27/02 - D.VERSEGHY. ESTIMATE FRACTIONAL CLOUD COVER
  !     *                         AND RAINFALL/SNOWFALL RATES
  !     *                         IF NECESSARY.
  !
  use classicParams,        only : abszero, TFREZ, RGAS, RGASV, RHOW
  use generalutilsGAT,         only : calcEsat

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in)  :: IPCP, NL, IL1, IL2
  integer             :: I
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: VPD   (NL)   !< Vapour pressure deficit of air \f$[mb] (e_d)\f$
  real, intent(out) :: TADP  (NL)   !< Dew point temperature of air [K]
  real, intent(out) :: PADRY (NL)   !< Partial pressure of dry air \f$[Pa] (p_{dry})\f$
  real, intent(out) :: RHOAIR(NL)   !< Density of air \f$[kg m^{-3}] (\rho_a)\f$
  real, intent(out) :: RHOSNI(NL)   !< Density of fresh snow \f$[kg m^{-3}] (\rho_{s, i})\f$
  real, intent(out) :: RPCP  (NL)   !< Calculated rainfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(out) :: TRPCP (NL)   !< Rainfall temperature over modelled area [C]
  real, intent(out) :: SPCP  (NL)   !< Calculated snowfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(out) :: TSPCP (NL)   !< Snowfall temperature over modelled area [C]
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: TA    (NL)   !< Air temperature at reference height \f$[K] (T_a)\f$
  real, intent(in) :: QA    (NL)   !< Specific humidity at reference height \f$[kg kg^{-1}] (q_a)\f$
  real, intent(in) :: PRESSG(NL)   !< Surface atmospheric pressure \f$[Pa] (p)\f$
  real, intent(in) :: PCPR  (NL)   !< Precipitation rate over modelled area
  !! \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: RRATE (NL)   !< Input rainfall rate over modelled area
  !! \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: SRATE (NL)   !< Input snowfall rate over modelled area
  !! \f$[kg m^{-2} s^{-1}]\f$
  !
  !     * WORK ARRAYS.
  !
  real :: PHASE (NL)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: EA, CA, CB, EASAT, CONST
  !
  !----------------------------------------------------------------
  !
  !     * CALCULATION OF ATMOSPHERIC INPUT VARIABLES.
  !
  !>
  !! In the first section, the air vapour pressure deficit, dry air
  !! pressure, air density and dew point temperature are calculated.
  !! The vapour pressure deficit \f$e_d\f$ (in units of mb) is obtained from
  !! the saturated and actual vapour pressures of the air, \f$e_a\f$ and
  !! \f$e_{a, sat}\f$ respectively (in units of Pa), as
  !! \f$e_d = [e_{a, sat} - e_a] /100.0\f$
  !!
  !! The air vapour pressure is obtained from the specific humidity \f$q_a\f$
  !! using the formula
  !!
  !! \f$e_a = q_a p /[0.622 + 0.378 q_a ]\f$
  !!
  !! where p is the surface atmospheric pressure. For the
  !! saturated vapour pressure, following Emanuel (1994) \cite Emanuel1994-dt
  !! \f$e_{a, sat}\f$ is from the temperature \f$T_a\f$ and the freezing
  !! point \f$T_f\f$:
  !!
  !! \f$e_{sat} = exp[53.67957 - 6743.769 / T - 4.8451 * ln(T)]       T \geq T_f\f$
  !!
  !! \f$e_{sat} = exp[23.33086 - 6111.72784 / T + 0.15215 * log(T)]    T < T_f\f$
  !!
  !! The partial pressure of dry air, \f$p_{dry}\f$, is obtained by subtracting
  !! \f$e_a\f$ from p, and the density of the air is calculated as the sum of
  !! the densities of the dry air and the water vapour:
  !! \f$\rho_a = p_{dry} / R_d T_a + e_a / R_v T_a\f$
  !!
  !! where \f$R_d\f$ and \f$R_v\f$ are the gas constants for dry air and water
  !! vapour respectively. The dew point temperature of the air is
  !! evaluated by substituting \f$e_a\f$ for \f$e_{a, sat}\f$ on the left-hand side of
  !! the appropriate equation above, and solving for \f$T_a\f$.
  !!
  !! Finally, if IPCP = 4, this indicates that the partitioning
  !! between rainfall and snowfall has been done outside of CLASS.
  !! The rainfall and snowfall rates RRATE and SRATE that have been
  !! passed into the subroutine are therefore simply assigned to RPCP
  !! and SPCP.
  !!
  do I = IL1,IL2 ! loop 100
    EA = QA(I) * PRESSG(I) / (0.622 + 0.378 * QA(I))
    EASAT = calcEsat(TA(I))
    VPD(I) = MAX(0.0,(EASAT - EA) / 100.0)
    PADRY(I) = PRESSG(I) - EA
    RHOAIR(I) = PADRY(I) / (RGAS * TA(I)) + EA / (RGASV * TA(I))
    CONST = LOG(EA / 611.0)
    TADP(I) = (CB * CONST - CA * TFREZ) / max((CONST - CA), abszero)
    !
    !     * DENSITY OF FRESH SNOW.
    !
    !>
    !! In the next section, the density of fresh snow \f$\rho_{s, i}\f$ is
    !! determined as an empirical function of the air temperature.
    !! For temperatures below 0 C, an equation presented by Hedstrom
    !! and Pomeroy (1998) is used. For temperatures >= 0 C, a
    !! relation following Pomeroy and Gray (1995) is used, with an
    !! upper limit of 200 kg m-3:
    !!
    !! \f$\rho_{s, i} = 67.92 + 51.25 exp[(T_a – T_f)/2.59]\f$ \f$T_a < T_f\f$
    !!
    !! \f$\rho_{s, i} = 119.17 + 20.0 (T_a – T_f)\f$           \f$T_a \geq T_f\f$
    !!
    if (TA(I) <= TFREZ) then
      RHOSNI(I) = 67.92 + 51.25 * EXP((TA(I) - TFREZ) / 2.59) 
    else
      RHOSNI(I) = MIN((119.17 + 20.0 * (TA(I) - TFREZ)),200.0)
    end if
    !
    !     * PRECIPITATION PARTITIONING BETWEEN RAIN AND SNOW.
    !
    !>
    !! In the last section, the partitioning of precipitation
    !! between rainfall RPCP and snowfall SPCP is addressed. Four
    !! options for doing so are provided; the user’s selection is
    !! indicated by the flag IPCP. In each case the rainfall and
    !! snowfall rates are converted to units of \f$ m s^{-1}\f$, by dividing
    !! by the density of water in the case of rain and by \f$\rho_{s, i}\f$ in
    !! the case of snow. The rainfall temperature is set to the
    !! maximum of 0 C and \f$T_a\f$, and the snowfall temperature to the
    !! minimum of 0 C and \f$T_a\f$. If IPCP = 1, the precipitation is
    !! simply diagnosed as rain if the air temperature is greater
    !! than 0 C, and as snow otherwise. If IPCP = 2, an empirical
    !! relation developed by Brown (2001) is used, where the
    !! precipitation is entirely snowfall when \f$T_a \leq 0 C\f$, and
    !! entirely rainfall when \f$T_a \geq 2.0 C\f$, and varies linearly
    !! between the two, with an equal mix of rain and snow at
    !! \f$T_a = 1.0 C\f$. If IPCP = 3, the precipitation is assumed to be
    !! entirely snowfall when \f$T_a \leq 0 C\f$, and entirely rainfall when
    !! \f$T_a \geq 6.0 C\f$, and between the two a polynomial function
    !! presented by Auer (1974) \cite Auer1974-fd is used, relating the fraction of
    !! the precipitation that is snowfall, \f$X_{sf}\f$, to \f$T_a\f$:
    !! \f[
    !! X_{sf} = [0.0202 T_a^6 – 0.3660 T_a^5 + 2.0399 T_a^4 – 1.5089 T_a^3
    !!      – 15.038 T_a^2 + 4.6664 T_a + 100.0]/100.0 \f]
    !!
    RPCP (I) = 0.0
    TRPCP(I) = 0.0
    SPCP (I) = 0.0
    TSPCP(I) = 0.0
    if (PCPR(I) > 1.0E-8) then
      if (IPCP == 1) then
        if (TA(I) > TFREZ) then
          RPCP (I) = PCPR(I) / RHOW
          TRPCP(I) = MAX((TA(I) - TFREZ),0.0)
        else
          SPCP (I) = PCPR(I) / RHOSNI(I)
          TSPCP(I) = MIN((TA(I) - TFREZ),0.0)
        end if
      else if (IPCP == 2) then
        if (TA(I) <= TFREZ) then
          PHASE(I) = 1.0
        else if (TA(I) >= (TFREZ + 2.0)) then
          PHASE(I) = 0.0
        else
          PHASE(I) = 1.0 - 0.5 * (TA(I) - TFREZ)
        end if
        RPCP(I) = (1.0 - PHASE(I)) * PCPR(I) / RHOW
        if (RPCP(I) > 0.0) TRPCP(I) = MAX((TA(I) - TFREZ),0.0)
        SPCP(I) = PHASE(I) * PCPR(I) / RHOSNI(I)
        if (SPCP(I) > 0.0) TSPCP(I) = MIN((TA(I) - TFREZ),0.0)
      else if (IPCP == 3) then
        if (TA(I) <= TFREZ) then
          PHASE(I) = 1.0
        else if (TA(I) >= (TFREZ + 6.0)) then
          PHASE(I) = 0.0
        else
          PHASE(I) = (0.0202 * (TA(I) - TFREZ) ** 6 - 0.3660 * & 
                     (TA(I) - TFREZ) ** 5 + 2.0399 * (TA(I) - TFREZ) ** 4 - &
                     1.5089 * (TA(I) - TFREZ) ** 3 - 15.038 * &
                     (TA(I) - TFREZ) ** 2 + 4.6664 * (TA(I) - TFREZ) + 100.0) / &
                     100.0
          PHASE(I) = MAX(0.0,MIN(1.0,PHASE(I)))
        end if
        RPCP(I) = (1.0 - PHASE(I)) * PCPR(I) / RHOW
        if (RPCP(I) > 0.0) TRPCP(I) = MAX((TA(I) - TFREZ),0.0)
        SPCP(I) = PHASE(I) * PCPR(I) / RHOSNI(I)
        if (SPCP(I) > 0.0) TSPCP(I) = MIN((TA(I) - TFREZ),0.0)
      else if (IPCP == 4) then
        RPCP(I) = RRATE(I) / RHOW
        if (RPCP(I) > 0.0) TRPCP(I) = MAX((TA(I) - TFREZ),0.0)
        SPCP(I) = SRATE(I) / RHOSNI(I)
        if (SPCP(I) > 0.0) TSPCP(I) = MIN((TA(I) - TFREZ),0.0)
      end if
    end if
  end do ! loop 100
  !
  return
end subroutine atmosphericVarsCalc
