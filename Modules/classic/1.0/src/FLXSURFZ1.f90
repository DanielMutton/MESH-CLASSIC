!> \file
!> Calculates surface layer transfer coefficients and fluxes
!! FLXSURFZ is a variant of FLXSURF3 that permits to input
!! wind and temperature (humidity) at different levels
!! @author Y. Delage, G. Pellerin, B. Bilodeau, M. Desgagne, R. Sarrazin, C. Girard, D. Verseghy, V. Fortin, J. P. Paquin, R. Harvey, J. Melton
subroutine FLXSURFZ1 (CDM, CDH, CTU, RIB, FTEMP, FVAP, ILMO, &
                     UE, FCOR, TA, QA, ZU, ZT, VA, &
                     TG, QG, H, Z0, Z0T, &
                     LZZ0, LZZ0T, FM, FH, N, IL1, IL2, FI, ITER, JL)

  use classicParams,  only : AS, CI, BS, BETA, FACTN, HMIN, DELTA, GRAV, &
                             VKC, ASX

  implicit none

  !
  ! Author
  !          Y.Delage (Jul 1990)
  ! Revision
  ! 001      G. Pellerin (Jun 94) New function for unstable case
  ! 002      G. Pellerin (Jui 94) New formulation for stable case
  ! 003      B. Bilodeau (Nov 95) Replace VK by KARMAN
  ! 004      M. Desgagne (Dec 95) Add safety code in function ff
  !                               and ensures that RIB is non zero
  ! 005      R. Sarrazin (Jan 96) Correction for H
  ! 006      C. Girard (Nov 95) - Diffuse T instead of Tv
  ! 007      G. Pellerin (Feb 96) Revised calculation for H (stable)
  ! 008      G. Pellerin (Feb 96) Remove corrective terms to CTU
  ! 009      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
  ! 010      Y. Delage (Feb 98) - Addition of HMIN
  ! 011      Y. Delage (Sept 00) - Set top of surface layer at ZU +Z0
  !                              - Output UE instead of UE**2
  !                              - Initialise ILMO and H
  !                              - Change iteration scheme for stable case
  !                              - Introduce log-linear profile for near-
  !                                 neutral stable cases
  !                              - set VAMIN inside flxsurf
  ! 012      Y. Delage (Oct 00) - Input of wind and temperatur/humidity
  !                                at different levels
  ! 013      D. Verseghy (Nov 02) - Add calculation of CDH
  ! 014      D. Verseghy (Nov 04) - Pass in JL for troubleshooting
  ! 015      V. Fortin (Jun 06) - Modify calculations of LZZ0 and LZZ0T
  !                               to avoid numerical problems with log
  ! 016      D. Verseghy (Jun 06) - Loop over IL1-IL2 instead of 1-N
  ! 017      J.P. Paquin (Aug 08) - "Synchronization" with flxsurf3
  !                               - Insert code from stabfunc2.cdk (v4.5)
  !                               - VAMIN=2.5 m/s to be coherent with ISBA
  !                                 (temporary measure VAMIN should be added
  !                                  to physics constants)
  !                                 (changes implemented by L.Duarte on Oct 08)
  ! 018      R. Harvey   (Oct 08) - Add fractional subregion cover FI
  !                                 (in addition to ITER) to control
  !                                 calculations over relevant points.
  ! 019      B. Dugas (after L. Spacek for flxsurf3 v_5.0.2) (Jan 09)
  !                               - Correction of the log-linear profile
  !                               - Double precision for rib calculations
  !                               - VAMIN is now retreived from CLASSD3
  ! 020      J. Melton            - Made some reals coded as ints, explicitly real
  !
  ! Object
  !          to calculate surface layer transfer coefficients and fluxes
  !          FLXSURFZ is a variant of FLXSURF3 that permits to input
  !          wind and temperature (humidity) at different levels
  !
  ! Arguments
  !
  !          - Output -
  real, intent(out) :: CDM(N) !< transfer coefficient of momentum squared
  real, intent(inout) :: CTU(N) !< transfer coefficient of temperature times UE
  real, intent(inout) :: RIB(N) !< bulk Richardson number
  real, intent(out) :: FTEMP(N) !< temperature flux
  real, intent(out) :: FVAP(N) !< vapor flux
  real, intent(inout) :: ILMO(N) !< (1/length of Monin-Obukov)
  real, intent(inout) :: UE(N) !< friction velocity
  real, intent(inout) :: H(N) !< height of the boundary layer
  real, intent(inout) :: FM(N) !< momentum stability function
  real, intent(inout) :: FH(N) !< heat stability function
  real, intent(inout) :: LZZ0(N) !< log ((zu+z0)/z0)
  real, intent(inout) :: LZZ0T(N) !< log ((zt+z0)/z0t)
  !          - Input -
  real, intent(in) :: FCOR(N) !< Coriolis factor
  real, intent(in) :: ZU(N) !< height of wind input
  real, intent(in) :: ZT(N) !< height of temperature and humidity input
  real, intent(in) :: TA(N) !< potential temperature at first predictive level above surface
  real, intent(in) :: QA(N) !< specific humidity     "    "      "        "      "     "
  real, intent(in) :: VA(N) !< wind speed            "    "      "        "      "     "
  real, intent(in) :: TG(N) !< surface temperature
  real, intent(in) :: QG(N) !< specific humidity at the surface
  real, intent(in) :: Z0(N) !< roughness length for momentum      flux calculations
  real, intent(in) :: Z0T(N) !< roughness length for heat/moisture flux calculations
  integer, intent(in) :: N !< horizontal dimension
  !
  ! Notes
  !          SEE DELAGE AND GIRARD BLM 58 (19-31)
  !                "       BLM 82 (23-48)
  !
  !     DIVERSES CONSTANTES PHYSIQUES
  !
  integer, intent(in) :: IL1,IL2,ITER(N),JL
  real, intent(out) :: CDH(N)
  real, intent(in) :: FI(N)


  integer :: J
  integer :: IT,ITMAX
  real :: HMAX,CORMIN,EPSLN
  real :: RAC3,CM,CT,ZP
  real :: F,G,DG
  real :: HI,HE,HS,unsl
  real * 8 :: DTHV,TVA,TVS
  real :: HL,U
  real :: CS,XX,XX0,YY,YY0
  real :: ZB,DD,ILMOX
  real :: DF,ZZ,betsasx
  real :: aa,bb,cc
  save HMAX,CORMIN,EPSLN
  save ITMAX
  real, save :: VMODMIN = - 1.0
  !
  real            :: VAMIN
  !     COMMON /CLASSD3/VAMIN
  data VAMIN / 1.0 /
  !
  data CORMIN,HMAX / 0.7E-4,1500.0 /
  data ITMAX / 3 /
  data EPSLN / 1.0e-05 /
  !
  DF(ZZ) = (1.0 - ZZ * HI) * sqrt(1.0 + (4.0 * AS * BETA * ilmo(j)) * ZZ / (1.0 - ZZ * HI))
  !
  RAC3 = sqrt(3.0)
  CS = AS * 2.5
  if (VMODMIN < 0.0) VMODMIN = SQRT(VAMIN)
  betsasx = 1. / asx
  !
  do J = IL1,IL2
    if (FI(J) > 0. .and. ITER(J) == 1) then
      !
      !  CALCULATE THE RICHARDSON NUMBER
      ZP = ZU(J) ** 2 / (ZT(J) + Z0(J) - Z0T(J))
      u = max(vmodmin,va(j))
      tva = (1.d0 + DELTA * QA(J)) * TA(J)
      tvs = (1.d0 + DELTA * QG(J)) * TG(J)
      dthv = tva - tvs
      RIB(J) = GRAV / (tvs + 0.5 * dthv) * ZP * dthv / (u * u)
      if (rib(j) >= 0.0) rib(j) = max(rib(j),EPSLN)
      if (rib(j) < 0.0) rib(j) = min(rib(j), - EPSLN)
      !
      !  FIRST APPROXIMATION TO ILMO
      LZZ0(J) = LOG(Z0(J) + ZU(J)) - LOG(Z0(J))
      LZZ0T(J) = LOG(ZT(J) + Z0(J)) - LOG(Z0T(J))
      if (RIB(J) > 0.) then
        FM(J) = LZZ0(J) + CS * RIB(J) / max(2.0 * z0(j),1.0)
        FH(J) = BETA * (LZZ0T(J) + CS * RIB(J)) / &
                max(sqrt(z0(j) * z0t(j)),1.0)
        ILMO(J) = RIB(J) * FM(J) * FM(J) / (ZP * FH(J))
        F = MAX(ABS(FCOR(J)),CORMIN)
        H(J) = BS * sqrt(VKC * u / (ILMO(J) * F * fm(j)))
      else
        FM(J) = LZZ0(J) - min(0.7 + log(1.0 - rib(j)),LZZ0(J) - 1.0)
        FH(J) = BETA * (LZZ0T(J) - min(0.7 + log(1.0 - rib(j)),LZZ0T(J) - 1.0))
      end if
      ILMO(J) = RIB(J) * FM(J) * FM(J) / (ZP * FH(J))
    end if
  end do

  ! - - - - - - - - -  BEGINNING OF ITERATION LOOP - - - - - - - - - - -
  do IT = 1,ITMAX ! loop 35
    do J = IL1,IL2
      if (FI(J) > 0. .and. ITER(J) == 1) then
        u = max(vmodmin,va(j))
        ZP = ZU(J) ** 2 / (ZT(J) + Z0(J) - Z0T(J))
        if (RIB(J) > 0.) then
          !----------------------------------------------------------------------
          !  STABLE CASE
          ILMO(J) = max(EPSLN,ILMO(J))
          hl = (ZU(J) + 10.0 * Z0(J)) * FACTN
          F = MAX(ABS(FCOR(J)),CORMIN)
          hs = BS * sqrt(VKC * u / (ILMO(J) * F * fm(j)))
          H(J) = MAX(HMIN,hs,hl,factn / (4.0 * AS * BETA * ILMO(J)))
          HI = 1.0 / H(J)
          unsl = ILMO(J)
          ! CDIR IEXPAND
          fm(J) = LZZ0(J) + psi(ZU(J) + Z0(J),hi,unsl) - psi(Z0(J),hi,unsl)
          ! CDIR IEXPAND
          fh(J) = BETA * (LZZ0T(J) + psi(ZT(J) + Z0(J),hi,unsl) &
                  - psi(Z0T(J),hi,unsl))
          DG = - ZP * FH(J) / (FM(J) * FM(J)) * (1.0 + beta * (DF(ZT(J) + Z0(J)) &
               - DF(Z0T(J))) / (2.0 * FH(J)) - (DF(ZU(J) + Z0(J)) &
               - DF(Z0(J))) / FM(J))
          !----------------------------------------------------------------------
          !  UNSTABLE CASE
        else
          ILMO(J) = MIN(0.,ILMO(J))
          ! CDIR IEXPAND
          FM(J) = fmi(zu(j) + z0(j),z0 (j),lzz0 (j),ilmo(j),xx,xx0)
          ! CDIR IEXPAND
          FH(J) = fhi(zt(j) + z0(j),z0t(j),lzz0t(j),ilmo(j),yy,yy0)
          DG = - ZP * FH(J) / (FM(J) * FM(J)) * (1.0 + beta / FH(J) * (1.0 / YY - 1.0 / YY0) &
               - 2.0 / FM(J) * (1.0 / XX - 1.0 / XX0))
        end if
        !----------------------------------------------------------------------
        if (IT < ITMAX) then
          G = RIB(J) - FH(J) / (FM(J) * FM(J)) * ZP * ILMO(J)
          ILMO(J) = ILMO(J) - G / DG
        end if
      end if
    end do
  end do ! loop 35
  ! - - - - - -  - - - END OF ITERATION LOOP - - - - - - - - - - - - - -
  do J = IL1,IL2
    if (FI(J) > 0. .and. ITER(J) == 1) then
      u = max(vmodmin,va(j))
      if (asx < as) then
        !----------------------------------------------------------------------
        !  CALCULATE ILMO AND STABILITY FUNCTIONS FROM LOG-LINEAR PROFILE
        !     (SOLUTION OF A QUADRATIC EQATION)
        !
        zb = zu(j) / (zt(j) + z0(j) - z0t(j))
        !  DISCRIMINANT
        dd = (beta * lzz0t(j) * zb) ** 2 - 4.0 * rib(j) * asx * lzz0(j) * &
             (beta * lzz0t(j) * zb - lzz0(j))
        if (rib(j) > 0. .and. rib(j) < betsasx .and. dd >= 0.) then
          !  COEFFICIENTS
          aa = asx * asx * rib(j) - asx
          bb = - beta * lzz0t(j) * zb + 2.0 * rib(j) * asx * lzz0(j)
          cc = rib(j) * lzz0(j) ** 2
          !  SOLUTION
          if (bb >= 0.) then
            ilmox = ( - bb - sqrt(dd)) &
                    / (2.0 * zu(j) * aa)
          else
            ilmox = 2.0 * cc / (zu(j) * ( - bb + sqrt(dd)))
          end if
          if (ilmox < ilmo(j)) then
            ilmo(j) = ilmox
            fm(j) = lzz0(j) + asx * zu(j) * ilmox
            fh(j) = beta * lzz0t(j) + asx * (zt(j) + z0(j) - z0t(j)) * ilmox
          end if
        end if
      end if
      !----------------------------------------------------------------------
      CM = VKC / FM(J)
      CT = VKC / FH(J)
      UE(J) = u * CM
      CDM(J) = CM ** 2
      CTU(J) = CT * UE(J)
      CDH(J) = CM * CT
      if (rib(j) > 0.0) then
        !             cas stable
        H(J) = MIN(H(J),hmax)
      else
        !             cas instable
        F = MAX(ABS(FCOR(J)),CORMIN)
        he = max(HMIN,0.3 * UE(J) / F)
        H(J) = MIN(he,hmax)
      end if
      FTEMP(J) = - CTU(J) * (TA(J) - TG(J))
      FVAP(J) = - CTU(J) * (QA(J) - QG(J))
    end if
  end do ! loop 80
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
end subroutine FLXSURFZ1
