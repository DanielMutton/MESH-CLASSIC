!> \file
!> Calculates screen relative humidity based on input screen
!! temperature, screen specific humidity and surface pressure.
!!
!! The formulae used here are consistent with that used elsewhere
!! in the gcm physics.
!! @author M. Lazare, V. Fortin
!
subroutine screenRelativeHumidity (SRH, ST, SQ, PRESSG, FMASK, ILG, IL1, IL2) ! Formerly SCREENRH
  !
  !     * DEC 16, 2014 - V.FORTIN.  REMOVE MERGE IN CALCULATION OF FRACW.
  !     * APR 30, 2009 - M.LAZARE.
  !
  use classicParams, only : EPS1, EPS2, T1S, RW1, RW2, RW3, RI1, RI2, RI3

  implicit none
  !
  !     * OUTPUT FIELD:
  !
  real, intent(out)   :: SRH(ILG)
  !
  !     * INPUT FIELDS.
  !
  real, intent(in)   :: ST(ILG), SQ(ILG), PRESSG(ILG), FMASK(ILG)
  !
  real :: FACTE, EPSLIM, FRACW, ETMP, ESTREF, ESAT, QSW
  !      REAL A, B, EPS1, EPS2, T1S, T2S, AI, BI, AW, BW, SLP
  !      REAL RW1, RW2, RW3, RI1, RI2, RI3
  real :: ESW, ESI, ESTEFF, TTT, UUU
  !
  integer, intent(in) :: ILG, IL1, IL2
  integer             :: IL
  !
  !     * COMPUTES THE SATURATION VAPOUR PRESSURE OVER WATER OR ICE.
  !
  ESW(TTT)        = EXP(RW1 + RW2 / TTT) * TTT ** RW3
  ESI(TTT)        = EXP(RI1 + RI2 / TTT) * TTT ** RI3
  ESTEFF(TTT,UUU) = UUU * ESW(TTT) + (1. - UUU) * ESI(TTT)
  !========================================================================

  EPSLIM = 0.001
  FACTE = 1. / EPS1 - 1.
  do IL = IL1,IL2
    if (FMASK(IL) > 0.) then
      !
      !       * COMPUTE THE FRACTIONAL PROBABILITY OF WATER PHASE
      !       * EXISTING AS A FUNCTION OF TEMPERATURE (FROM ROCKEL,
      !       * RASCHKE AND WEYRES, BEITR. PHYS. ATMOSPH., 1991.)
      !
      if (ST(IL) >= T1S) then
        FRACW = 1.0
      else
        FRACW = 0.0059 + 0.9941 * EXP( - 0.003102 * (T1S - ST(IL)) ** 2)
      end if
      !
      ETMP = ESTEFF(ST(IL),FRACW)

      ESTREF = 0.01 * PRESSG(IL) * (1. - EPSLIM) / (1. - EPSLIM * EPS2)
      if (ETMP < ESTREF) then
        ESAT = ETMP
      else
        ESAT = ESTREF
      end if
      !
      QSW = EPS1 * ESAT / (0.01 * PRESSG(IL) - EPS2 * ESAT)
      SRH(IL) = MIN(MAX((SQ(IL) * (1. + QSW * FACTE)) &
                / (QSW * (1. + SQ(IL) * FACTE)),0.),1.)
    end if
  end do
  !
  return
end subroutine screenRelativeHumidity
