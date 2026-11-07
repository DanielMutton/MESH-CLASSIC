!> \file
!> Calculate the diagnostic values of U, V, T, Q
!! near the surface (ZU and ZT)
!! @author Y. Delage, G. Pellerin, B. Bilodeau, R. Sarrazin, G. Pellerin, D. Verseghy,
!! M. Mackay, F. Seglenieks, P.Bartlett, E.Chan, B.Dugas, J.Melton
!
subroutine DIASURFZ1 (UZ, VZ, TZ, QZ, NI, U, V, TG, QG, Z0, Z0T, ILMO, ZA, &
                     H, UE, FTEMP, FVAP, ZU, ZT, LAT, F, IL1, IL2, JL)

  use classicParams, only : AS, ASX, CI, BETA, FACTN, HMIN, ANGMAX, &
                            GRAV, SPHAIR, VKC

  implicit none

  integer, intent(in) :: NI, JL
  real, intent(in)    :: ZT(NI), ZU(NI)
  real, intent(inout) :: UZ(NI), VZ(NI), TZ(NI), QZ(NI)
  real, intent(in)    :: ZA(NI), U(NI), V(NI)
  real, intent(in)    :: TG(NI), QG(NI), UE(NI), FTEMP(NI), FVAP(NI)
  real, intent(in)    :: ILMO(NI), Z0T(NI), Z0(NI), H(NI), F(NI)
  real :: LAT(NI)
  ! Author
  !          Yves Delage  (Aug1990)
  !
  ! Revision
  ! 001      G. Pellerin(JUN94)
  !          Adaptation to new surface formulation
  ! 002      B. Bilodeau (Nov 95) - Replace VK by KARMAN
  ! 003      R. Sarrazin (Jan 96) - Prevent problems if zu < za
  ! 004      G. Pellerin (Feb 96) - Rewrite stable formulation
  ! 005      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
  ! 006      Y. Delage (Feb 98) - Addition of HMIN
  ! 007      Y. Delage (Sept 00) - Change UE2 by UE
  !                              - Introduce log-linear profile for near-
  !                                 neutral cases
  ! 008      D. Verseghy (Nov 02) - Remove unused constant CLM
  !                                 from common block CLASSD2
  ! 009      M. Mackay (Nov 04) - Change all occurrences of ALOG
  !                               to LOG for greater portability.
  ! 010      F. SeglenieKs (Mar 05) - Declare LAT as REAL(8) for
  !                                   consistency
  ! 011      P.Bartlett (Mar 06) - Set HI to zero for unstable case
  ! 012      E.Chan (Nov 06) - Bracket entire subroutine loop with
  !                            IF (F(J)>0.)
  ! 013      D.Verseghy (Nov 06) - Convert LAT to regular precision
  ! 014      B.Dugas (Jan 09) - "Synchronization" with diasurf2
  ! 015      J.Melton (Sep 12) - Made some numbers explictly reals, they
  !                              were written as ints.
  !
  ! Object
  !          to calculate the diagnostic values of U, V, T, Q
  !          near the surface (ZU and ZT)
  !
  ! Arguments
  !
  !          - Output -
  ! UZ       U component of the wind at Z=ZU
  ! VZ       V component of the wind at Z=ZU
  ! TZ       temperature in kelvins at Z=ZT
  ! QZ       specific humidity at Z=ZT
  !
  !          - Input -
  ! NI       number of points to process
  ! U        U component of wind at Z=ZA
  ! V        V component of wind at Z=ZA
  ! TG       temperature at the surface (Z=0) in Kelvins
  ! QG       specific humidity
  ! PS       surface pressure at the surface
  ! ILMO     inverse of MONIN-OBUKHOV lenth
  ! H        height of boundary layer
  ! UE       friction velocity
  ! Z0       roughness lenth for winds
  ! Z0T      roughness lenth for temperature and moisture
  ! FTEMP    temperature flux at surface
  ! FVAP     vapor flux at surface
  ! ZA       heights of first model level above ground
  ! ZU       heights for computation of wind components
  ! ZT       heights for computation of temperature and moisture
  ! LAT      LATITUDE
  ! F        Fraction of surface type being studied

  real :: ANG, ANGI, VITS, LZZ0, LZZ0T
  real :: CT, DANG, CM
  real :: XX, XX0, YY, YY0, fh, fm, hi
  real :: RAC3
  integer :: J, IL1, IL2

  RAC3 = SQRT(3.)

  do J = IL1,IL2
    if (F(J) > 0.0) then

      LZZ0T = LOG((ZT(J) + Z0(J)) / Z0T(J))
      LZZ0 = LOG(ZU(J) / Z0(J) + 1.0)
      if (ILMO(J) <= 0.) then
        !---------------------------------------------------------------------
        !                      UNSTABLE CASE
        !
        hi = 0.
        ! CDIR IEXPAND
        fh = fhi(ZT(J) + Z0(J),Z0T(j),LZZ0T,ILMO(J),YY,YY0)
        ! CDIR IEXPAND
        fm = fmi(ZU(J) + Z0(J),Z0 (J),LZZ0,ILMO(J),XX,XX0)
      else
        !---------------------------------------------------------------------
        !                        STABLE CASE
        hi = 1.0 / MAX(HMIN,H(J),(ZA(J) + 10.0 * Z0(J)) * factn,factn / &
             (4.0 * AS * BETA * ilmo(j)))
        ! CDIR IEXPAND
        fh = BETA * (LZZ0T + min( psi(ZT(J) + Z0(J),HI,ILMO(J)) &
             - psi(Z0T(J),HI,ILMO(J)), &
             ASX * ILMO(J) * (ZT(J) + Z0(J) - Z0T(J))))
        ! CDIR IEXPAND
        fm = LZZ0 + min( psi(zu(J) + Z0(J),HI,ILMO(J)) &
             - psi(Z0(J),HI,ILMO(J)), &
             ASX * ILMO(J) * ZU(J))
      end if
      !---------------------------------------------------------------------
      ! CT=KARMAN/FH
      ! CM=KARMAN/FM
      CT = VKC / FH
      CM = VKC / FM

      TZ(J) = TZ(J) + F(J) * (TG(J) - FTEMP(J) / (CT * UE(J)) - GRAV / SPHAIR * ZT(J))
      QZ(J) = QZ(J) + F(J) * (QG(J) - FVAP(J) / (CT * UE(J)))
      VITS = UE(J) / CM

      ! CALCULATE WIND DIRECTION CHANGE FROM TOP OF SURFACE LAYER
      DANG = (ZA(J) - ZU(J)) * HI * ANGMAX * SIN(LAT(J))
      ANGI = ATAN2(V(J),SIGN(ABS(U(J)) + 1.e-05,U(J)))
      if (ILMO(J) > 0.) then
        ANG = ANGI + DANG
      else
        ANG = ANGI
      end if

      UZ(J) = UZ(J) + F(J) * VITS * COS(ANG)
      VZ(J) = VZ(J) + F(J) * VITS * SIN(ANG)

    end if
  end do

  return
contains

  !   The following code is taken from the RPN/CMC physics library file
  !   /usr/local/env/armnlib/modeles/PHY_shared/ops/v_4.5/RCS/stabfunc2.cdk, v

  !   Internal function FMI
  !   Stability function for momentum in the unstable regime (ilmo<0)
  !   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 19
  !
  real function FMI(Z2, Z02, LZZ02, ILMO2, X, X0)
    implicit none
    !
    real, intent(in) :: Z2, Z02, LZZ02, ILMO2
    real, intent(out) :: X, X0
    !
    X = (1.0 - CI * Z2 * BETA * ILMO2) ** (0.16666666)
    X0 = (1.0 - CI * Z02 * BETA * ILMO2) ** (0.16666666)
    FMI = LZZ02 + LOG((X0 + 1.0) ** 2 * SQRT(X0 ** 2 - X0 + 1.0) * (X0 ** 2 + X0 + 1.0) ** 1.5 &
          / ((X + 1.0) ** 2 * SQRT(X ** 2 - X + 1.0) * (X ** 2 + X + 1.0) ** 1.5)) &
          + RAC3 * ATAN(RAC3 * ((X ** 2 - 1.0) * X0 - (X0 ** 2 - 1.0) * X) / &
          ((X0 ** 2 - 1.0) * (X ** 2 - 1.0) + 3.0 * X * X0))
    !
    return
  end function FMI
  !
  !   Internal function FHI
  !   Stability function for heat and moisture in the unstable regime (ilmo<0)
  !   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 17
  !
  real function FHI(Z2, Z0T2, LZZ0T2, ILMO2, Y, Y0)
    implicit none
    !
    real, intent(in) :: Z2, Z0T2, LZZ0T2, ILMO2
    real, intent(out) :: Y, Y0
    !
    Y = (1.0 - CI * Z2  * BETA * ILMO2) ** (0.33333333)
    Y0 = (1.0 - CI * Z0T2 * BETA * ILMO2) ** (0.33333333)
    FHI = BETA * (LZZ0T2 + 1.5 * LOG((Y0 ** 2 + Y0 + 1.0) / (Y ** 2 + Y + 1.0)) + RAC3 * &
          ATAN(RAC3 * 2.0 * (Y - Y0) / ((2.0 * Y0 + 1.0) * (2.0 * Y + 1.0) + 3.0)))
    !
    return
  end function FHI
  !
  !   Internal function psi
  !   Stability function for momentum in the stable regime (unsl>0)
  !   Reference :  Y. Delage, BLM, 82 (p23-48) (Eqs.33-37)
  !
  real function PSI(Z2, HI2, ILMO2)
    implicit none
    !
    real :: a, b, c, d
    real, intent(in) :: ILMO2, Z2, HI2
    !
    d = 4.0 * AS * BETA * ILMO2
    c = d * hi2 - hi2 ** 2
    b = d - 2.0 * hi2
    a = sqrt(1.0 + b * z2 - c * z2 ** 2)
    psi = 0.5 * (a - z2 * hi2 - log(1.0 + b * z2 * 0.5 + a) - &
          b / (2.0 * sqrt(c)) * asin((b - 2.0 * c * z2) / d))
    !
    return
  end function PSI
end subroutine DIASURFZ1
