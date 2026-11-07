!> \file
!> Calculates drag coefficients and related variables
!! This work has been described in \cite Abdella1996-em.
!! @author K. Abdella, N. Macfarlane, M. Lazare, D. Verseghy, E. Chan
!
subroutine DRCOEF1 (CDM, CDH, RIB, CFLUX, QG, QA, ZOMIN, ZOHIN, &
                   CRIB, TVIRTG, TVIRTA, VA, FI, ITER, &
                   ILG, IL1, IL2)
  !
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * SEP 10/02 - K.ABDELLA.  BUGFIX IN CALCULATION OF "OLS" (2 PLACES).
  !     * MAY 22/02 - N.MCFARLANE.USE THE ASYMPTOTIC VALUE FOR CDH
  !     *                         GENERALLY TO LIMIT FLUXES AND PREVENT
  !     *                         OVERSHOOT. (COMMMENTED OUT FOR
  !     *                         OFF-LINE RUNS.)
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * OCT 26/99 - E. CHAN.    COMMENT OUT ARTIFICIAL DAMPING OF
  !     *                         TURBULENT FLUXES FOR STAND-ALONE TESTING.
  !     * JUL 18/97 - M. LAZARE, CLASS 2.7. PASS IN ADDITIONAL WORK ARRAYS
  !     *             D.VERSEGHY. ZOMIN AND ZOHIN TO USE INTERNALLY,
  !     *                         SO THAT INPUT ZOMIN AND ZOHIN DO NOT
  !     *                         CHANGE WHEN PASSED BACK TO THE
  !     *                         ITERATION (BUG FROM PREVIOUS
  !     *                         VERSIONS). PREVIOUS ZOM AND ZOH
  !     *                         BECOME WORK ARRAYS.
  !     * MAR 10/97 - M. LAZARE.  PASS IN QG AND QA AND USE TO ONLY
  !     *                         DAMP TURBULENT FLUXES USING CDHMOD
  !     *                         UNDER STABLE CONDITIONS IF MOISTURE
  !     *                         FLUX IS UPWARDS. ALSO, BETTER
  !     *                         DEFINITION OF SURFACE LAYER TOP
  !     *                         FROM K.ABDELLA.
  !     * MAY 21/96 - K. ABDELLA. MODIFICATION FOR FREE-CONVECTIVE
  !     *                         LIMIT ON UNSTABLE SIDE ADDED.
  !     * JAN 10/96 - K. ABDELLA. CORRECT ERROR IN AU1 (UNSTABLE
  !     *                         SIDE) AND PUT IN PRANDTL NUMBER
  !     *                         RANGE (0.74->1)/
  !     *                         "CDHMOD" USED ON BOTH STABLE AND
  !     *                         UNSTABLE SIDE, TO LIMIT FLUXES
  !     *                         OVER LONG TIMESTEP.
  !     * M. LAZARE - FEB 14/95.  USE VARIABLES "ZOLN" AND "ZMLN"
  !     *                         ON UNSTABLE SIDE, FOR OPTIMIZATION.
  !     *                         THIS IS PREVIOUS VERSION "DRCOEFX".
  !     * K. ABDELLA/M. LAZARE. - NOV 30/94.
  !
  use classicParams, only : GRAV, VKC

  implicit none

  !     * INTEGER CONSTANTS.
  integer, intent(in) :: ILG, IL1, IL2
  integer :: I

  !     * OUTPUT ARRAYS.
  real, intent(out)   :: CDM    (ILG) !< STABILITY-DEPENDENT DRAG COEFFICIENT FOR MOMENTUM.
  real, intent(inout) :: CDH    (ILG) !< STABILITY-DEPENDENT DRAG COEFFICIENT FOR HEAT.
  real, intent(inout) :: RIB    (ILG) !< BULK RICHARDSON NUMBER.
  real, intent(out) :: CFLUX  (ILG) !< CD * MOD(V), BOUNDED BY FREE-CONVECTIVE LIMIT.

  !     * INPUT ARRAYS.
  real, intent(in) :: ZOMIN  (ILG) !< ROUGHNESS HEIGHTS FOR MOMENTUM/HEAT NORMALIZED BY REFERENCE HEIGHT.
  real, intent(in) :: ZOHIN  (ILG) !< ROUGHNESS HEIGHTS FOR MOMENTUM/HEAT NORMALIZED BY REFERENCE HEIGHT.
  real, intent(in) :: CRIB   (ILG) !< -RGAS*SLTHKEF/(VA**2), WHERE SLTHKEF=-LOG(MAX(SGJ(ILEV), SHJ(ILEV)))
  real, intent(in) :: TVIRTG (ILG) !< "SURFACE" VIRTUAL TEMPERATURE.
  real, intent(in) :: TVIRTA (ILG) !< LOWEST LEVEL VIRTUAL TEMPERATURE.
  real, intent(in) :: VA     (ILG) !< AMPLITUDE OF LOWEST LEVEL WIND.
  real, intent(in) :: FI     (ILG) !< FRACTION OF SURFACE TYPE BEING STUDIED.
  real, intent(in) :: QG     (ILG) !< SATURATION SPECIFIC HUMIDITY AT GROUND TEMPERATURE.
  real, intent(in) :: QA     (ILG) !< LOWEST LEVEL SPECIFIC HUMIDITY.

  integer, intent(in) :: ITER(ILG) !< INDEX ARRAY INDICATING IF POINT IS UNDERGOING FURTHER ITERATION OR NOT.

  !     * WORK ARRAYS.
  !> ZOM/ZOH: WORK ARRAYS USED FOR SCALING ZOMIN/ZOHIN ON STABLE SIDE, AS PART OF CALCULATION.
  real :: ZOM    (ILG)
  real :: ZOH    (ILG)

  !     * TEMPORARY VARIABLES.
  real :: AA, AA1, BETA, PR, ZLEV, ZS, ZOLN, ZMLN, CPR, ZI, OLSF, OLFACT, &
          ZL, ZMOL, ZHOL, XM, XH, BH1, BH2, BH, WB, WSTAR, RIB0, WSPEED, &
          AU1, OLS, PSIM1, PSIM0, PSIH1, PSIH0, USTAR, TSTAR, WTS, AS1, &
          AS2, AS3, CLIMIT

  !-------------------------------------------------------------
  AA = 9.5285714
  AA1 = 14.285714
  BETA = 1.2
  PR = 1.
  !
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. ITER(I) == 1) then
      RIB(I) = CRIB(I) * (TVIRTG(I) - TVIRTA(I))
      if (RIB(I) >= 0.0) then
        ZLEV = - CRIB(I) * TVIRTA(I) * (VA(I) ** 2) / GRAV
        ZS = MAX(10.,5. * MAX(ZOMIN(I) * ZLEV,ZOHIN(I) * ZLEV))
        ZS = ZLEV * (1. + RIB(I)) / (1. + (ZLEV / ZS) * RIB(I))
        ZOM(I) = ZOMIN(I) * ZLEV / ZS
        ZOH(I) = ZOHIN(I) * ZLEV / ZS
        RIB(I) = RIB(I) * ZS / ZLEV
      else
        ZOM(I) = ZOMIN(I)
        ZOH(I) = ZOHIN(I)
      end if
      ZOLN = LOG(ZOH(I))
      ZMLN = LOG(ZOM(I))
      if (RIB(I) < 0.0) then
        CPR = MAX(ZOLN / ZMLN,0.74)
        CPR = MIN(CPR,1.0)
        ZI = 1000.0
        OLSF = BETA ** 3 * ZI * VKC ** 2 / ZMLN ** 3
        OLFACT = 1.7 * (LOG(1. + ZOM(I) / ZOH(I))) ** 0.5 + 0.9
        OLSF = OLSF * OLFACT
        ZL = - CRIB(I) * TVIRTA(I) * (VA(I) ** 2) / GRAV
        ZMOL = ZOM(I) * ZL / OLSF
        ZHOL = ZOH(I) * ZL / OLSF
        XM = (1.00 - 15.0 * ZMOL) ** (0.250)
        XH = (1.00 - 9.0 * ZHOL) ** 0.25
        BH1 = - LOG( - 2.41 * ZMOL) + LOG(((1. + XM) / 2.) ** 2 * (1. + XM ** 2) / 2.)
        BH1 = BH1 - 2. * ATAN(XM) + ATAN(1.) * 2.
        BH1 = BH1 ** 1.5
        BH2 = - LOG( - 0.25 * ZHOL) + 2. * LOG(((1.00 + XH ** 2) / 2.00))
        BH = VKC ** 3. * BETA ** 1.5 / (BH1 * (BH2) ** 1.5)
        WB = SQRT(GRAV * (TVIRTG(I) - TVIRTA(I)) * ZI / TVIRTG(I))
        WSTAR = BH ** (0.333333) * WB
        RIB0 = RIB(I)

        WSPEED = SQRT(VA(I) ** 2 + (BETA * WSTAR) ** 2)
        RIB(I) = RIB0 * VA(I) ** 2 / WSPEED ** 2
        AU1 = 1. + 5.0 * (ZOLN - ZMLN) * RIB(I) * (ZOH(I) / ZOM(I)) ** 0.25
        OLS = - RIB(I) * ZMLN ** 2 / (CPR * ZOLN) * (1.0 + AU1 / &
              (1.0 - RIB(I) / (ZOM(I) * ZOH(I)) ** 0.25))
        PSIM1 = LOG(((1.00 + (1.00 - 15.0 * OLS) ** 0.250) / 2.00) ** 2 * &
                (1.0 + (1.00 - 15.0 * OLS) ** 0.5) / 2.0) - 2.0 * ATAN( &
                (1.00 - 15.0 * OLS) ** 0.250) + ATAN(1.00) * 2.00
        PSIM0 = LOG(((1.00 + (1.00 - 15.0 * OLS * ZOM(I)) ** 0.250) / 2.00) ** 2 &
                * (1.0 + (1.00 - 15.0 * OLS * ZOM(I)) ** 0.5) / 2.0) - 2.0 * &
                ATAN((1.00 - 15.0 * OLS * ZOM(I)) ** 0.250) + ATAN(1.00) * 2.0
        PSIH1 = LOG(((1.00 + (1.00 - 9.0 * OLS) ** 0.50) / 2.00) ** 2)
        PSIH0 = LOG(((1.00 + (1.00 - 9.0 * OLS * ZOH(I)) ** 0.50) / 2.00) ** 2)

        USTAR = VKC / ( - ZMLN - PSIM1 + PSIM0)
        TSTAR = VKC / ( - ZOLN - PSIH1 + PSIH0)
        CDH(I) = USTAR * TSTAR / PR
        WTS = CDH(I) * WSPEED * (TVIRTG(I) - TVIRTA(I))
        WSTAR = (GRAV * ZI / TVIRTG(I) * WTS) ** (0.333333)

        WSPEED = SQRT(VA(I) ** 2 + (BETA * WSTAR) ** 2)
        RIB(I) = RIB0 * VA(I) ** 2 / WSPEED ** 2
        AU1 = 1. + 5.0 * (ZOLN - ZMLN) * RIB(I) * (ZOH(I) / ZOM(I)) ** 0.25
        OLS = - RIB(I) * ZMLN ** 2 / (CPR * ZOLN) * (1.0 + AU1 / &
              (1.0 - RIB(I) / (ZOM(I) * ZOH(I)) ** 0.25))
        PSIM1 = LOG(((1.00 + (1.00 - 15.0 * OLS) ** 0.250) / 2.00) ** 2 * &
                (1.0 + (1.00 - 15.0 * OLS) ** 0.5) / 2.0) - 2.0 * ATAN( &
                (1.00 - 15.0 * OLS) ** 0.250) + ATAN(1.00) * 2.00
        PSIM0 = LOG(((1.00 + (1.00 - 15.0 * OLS * ZOM(I)) ** 0.250) / 2.00) ** 2 &
                * (1.0 + (1.00 - 15.0 * OLS * ZOM(I)) ** 0.5) / 2.0) - 2.0 * &
                ATAN((1.00 - 15.0 * OLS * ZOM(I)) ** 0.250) + ATAN(1.00) * 2.0
        PSIH1 = LOG(((1.00 + (1.00 - 9.0 * OLS) ** 0.50) / 2.00) ** 2)
        PSIH0 = LOG(((1.00 + (1.00 - 9.0 * OLS * ZOH(I)) ** 0.50) / 2.00) ** 2)

      else

        WSPEED = VA(I)
        AS1 = 10.0 * ZMLN * (ZOM(I) - 1.0)
        AS2 = 5.00 / (2.0 - 8.53 * RIB(I) * EXP( - 3.35 * RIB(I)) + 0.05 * RIB(I) ** 2)
        ! <<<
        AS2 = AS2 * PR * SQRT( - ZMLN) / 2.
        AS3 = 27. / (8. * PR * PR)
        ! >>>
        OLS = RIB(I) * (ZMLN ** 2 + AS3 * AS1 * (RIB(I) ** 2 + AS2 * RIB(I))) &
              / (AS1 * RIB(I) - PR * ZOLN)
        PSIM1 = - 0.667 * (OLS - AA1) * EXP( - 0.35 * OLS) - AA - OLS
        PSIM0 = - 0.667 * (OLS * ZOM(I) - AA1) * EXP( - 0.35 * OLS * ZOM(I)) &
                - AA - OLS * ZOM(I)
        PSIH1 = - (1.0 + 2.0 * OLS / 3.0) ** 1.5 - 0.667 * (OLS - AA1) &
                * EXP( - 0.35 * OLS) - AA + 1.0
        PSIH0 = - (1.0 + 2.0 * OLS * ZOH(I) / 3.0) ** 1.5 - 0.667 * (OLS * ZOH(I) - AA1) &
                * EXP( - 0.35 * OLS * ZOH(I)) - AA + 1.0

      end if

      USTAR = VKC / ( - ZMLN - PSIM1 + PSIM0)
      TSTAR = VKC / ( - ZOLN - PSIH1 + PSIH0)

      CDM(I) = USTAR ** 2.0
      CDH(I) = USTAR * TSTAR / PR
      !
      !         * CALCULATE CD*MOD(V) UNDER FREE-CONVECTIVE LIMIT.
      !
      if (TVIRTG(I) > TVIRTA(I)) then
        CLIMIT = 1.9E-3 * (TVIRTG(I) - TVIRTA(I)) ** 0.333333
      else
        CLIMIT = 0.
      end if
      CFLUX(I) = MAX(CDH(I) * WSPEED,CLIMIT)
    end if
  end do ! loop 100

  return
end subroutine DRCOEF1
