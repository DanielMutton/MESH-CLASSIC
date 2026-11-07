!> \file
!! Checks for closure of surface water budget, and for
!! unphysical values of certain variables.
!! @author D. Verseghy, M. Lazare, B. Dugas
!
subroutine checkWaterBudget (ISFC, PCPR, EVAP, RUNOFF, WLOST, RAICAN, SNOCAN, & ! Formerly CHKWAT
                             RAICNI, SNOCNI, ZPOND, ZPONDI, THLIQ, THICE, &
                             THLIQI, THICEI, ZSNOW, RHOSNO, XSNOW, SNOWI, &
                             WSNOW, WSNOWI, FCS, FGS, FI, BAL, THPOR, THLMIN, &
                             DELZW, ISAND, IG, ILG, IL1, IL2, JL, N)

  !     * APR 28/10 - B.DUGAS.    INTRODUCE SEPARATE ACCURACY LIMITS
  !     *                         FOR BAL AND FOR THE OTHER CHECKS.
  !     * JUN 06/06 - D.VERSEGHY. MODIFY CHECK ON RUNOFF.
  !     * APR 03/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
  !     * SEP 26/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 300 LOOP.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 04/04 - D.VERSEGHY. RELAX ACCURACY LIMIT FOR WATER
  !     *                         BALANCE CHECK, CONSISTENT WITH
  !     *                         ROUNDOFF ERROR CONSTRAINTS.
  !     * JUN 25/02 - D.VERSEGHY. RENAME VARIABLES FOR CLARITY; UPDATES
  !     *                         CAUSED BY CONVERSION OF PONDING DEPTH
  !     *                         TO A PROGNOSTIC VARIABLE; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * MAY 23/02 - D.VERSEGHY. MOVE CALCULATION OF "XSNOW" INTO
  !     *                         THIS ROUTINE.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         RATIONALIZE USE OF WLOST.
  !     *                         ALSO INTRODUCE NEW VALUE OF ACCLMT
  !     *                         CORRESPONDING TO 3 MM/YR AS USED
  !     *                         BY THE PILPS COMMUNITY.
  !     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS.
  !     * JAN 31/94 - D.VERSEGHY. LOCAL VERSION FOR CLIMATE RESEARCH
  !     *                         NETWORK: CHECK RAICAN AND SNOCAN
  !     *                         AGAINST -ACCLMT INSTEAD OF AGAINST 0.
  !     * AUG 16/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         RETURN WATER BALANCE CHECK TO ROUTINE
  !     *                         USE (COMMENTED OUT IN PREVIOUS VERSION)
  !     *                         AND RENAME SUBROUTINE FROM "CHKVAL"
  !     *                         TO "CHKWAT". (now "checkWaterBudget" JUL 2019)
  !     * MAY 15/92 - M.LAZARE.   CLASS - VERSION 2.1.
  !     *                         MOISTURE BALANCE CHECKS EXTRACTED FROM
  !     *                         "waterBudgetDriver" AND VECTORIZED.
  !     * APR 11/89 - D.VERSEGHY. THE FOLLOWING MOISTURE BALANCE CHECKS
  !     *                         ARE CARRIED OUT: INTERCEPTED MOISTURE
  !     *                         STORES AND LOCAL RUNOFF MUST BE >=0
  !     *                         LIQUID SOIL LAYER MOISTURE STORES MUST
  !     *                         BE LESS THAN THE PORE VOLUME AND GREATER
  !     *                         THAN THE LIMITING VALUE "THLMIN"; FROZEN
  !     *                         SOIL LAYER MOISTURE STORES MUST BE LESS
  !     *                         THAN THE MAXIMUM AVAILABLE VOLUME (THE
  !     *                         PORE VOLUME - THLMIN) AND GE.0; AND THE
  !     *                         MOISTURE BALANCE OF THE TOTAL CANOPY/
  !     *                         SNOW/SOIL COLUMN MUST BE WITHIN A
  !     *                         SPECIFIED TOLERANCE.  THE TOLERANCE
  !     *                         LEVEL ADOPTED IS DESIGNATED BY "ACCLMT".
  !
  use classicParams,       only : DELT, RHOW, RHOICE

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ISFC  !< Type of surface (1 = canopy over snow, 2 = snow covered ground, 3 = canopy over bare ground, 4 = bare ground)
  integer, intent(in) :: IG, ILG, IL1, IL2, JL, N
  integer :: I, J, K
  !
  integer :: IPTBAD, JPTBAD, KPTBAD, IPTBDI, JPTBDI, KPTBDI, LPTBDI, &
             IPTBDJ, JPTBDJ, KPTBDJ, LPTBDJ
  !
  !     * INPUT FIELDS.

  real, intent(in) :: PCPR  (ILG)  !< Precipitation rate over modelled subarea \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: EVAP  (ILG)  !< Evapotranspiration rate over modelled subarea \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: RUNOFF(ILG)  !< Total runoff over modelled subarea [m]
  real, intent(in) :: WLOST (ILG)  !< Residual amount of water that cannot be
  !! supplied by surface stores \f$[kg m^{-2}]\f$
  real, intent(in) :: RAICAN(ILG)  !< Intercepted liquid water on canopy at end of
  !! time step \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCAN(ILG)  !< Intercepted frozen water on canopy at end of
  !! time step \f$[kg m^{-2}]\f$
  real, intent(in) :: RAICNI(ILG)  !< Intercepted liquid water on canopy at beginning
  !! of time step \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCNI(ILG)  !< Intercepted frozen water on canopy at beginning
  !! of time step \f$[kg m^{-2}]\f$
  real, intent(in) :: ZPOND (ILG)  !< Depth of ponded water on ground at end of time
  !! step [m]
  real, intent(in) :: ZPONDI(ILG)  !< Depth of ponded water on ground at beginning of
  !! time step [m]
  real, intent(in) :: ZSNOW (ILG)  !< Depth of snow pack [m]
  real, intent(in) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}]\f$
  real, intent(inout) :: XSNOW (ILG)  !< Switch to indicate presence of snow cover [ ]
  real, intent(in) :: SNOWI (ILG)  !< Snow pack mass at beginning of time step \f$[kg m^{-2}]\f$
  real, intent(in) :: WSNOW (ILG)  !< Liquid water content of snow pack at end of
  !! time step \f$[kg m^{-2}]\f$
  real, intent(in) :: WSNOWI(ILG)  !< Liquid water content of snow pack at beginning
  !! of time step \f$[kg m^{-2}]\f$
  real, intent(in) :: FCS   (ILG)  !< Fractional coverage of canopy over snow on
  !! modelled area [ ]
  real, intent(in) :: FGS   (ILG)  !< Fractional coverage of snow over bare ground on
  !! modelled area [ ]
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !! modelled area [ ]

  !
  real, intent(in) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !! layers at end of time step \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !! layers at end of time step \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLIQI(ILG,IG)   !< Volumetric frozen water content of soil
  !! layers at beginning of time step \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICEI(ILG,IG)   !< Volumetric frozen water content of soil
  !! layers at beginning of time step \f$[m^3 m^{-3}]\f$
  !
  !     * WORK ARRAYS.
  !
  real, intent(inout) :: BAL(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: ACCLMT, BALLMT, CANFAC, SNOFAC
  !
  !     * SOIL INFORMATION ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !! remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer [m]
  !
  integer, intent(in) :: ISAND (ILG,IG)
  !
  !      ACCLMT=3.0*DELT/3.1536E7
  ACCLMT = 1.0E-3
  BALLMT = 1.0E-1
  !-----------------------------------------------------------------------
  !>
  !! This subroutine is called from waterBudgetDriver to perform water balance
  !! checks for each of the four subareas. The flag ISFC indicates
  !! which subarea is being addressed: ISFC=1 for vegetation over
  !! snow, ISFC=2 for snow over bare ground, ISFC=3 for vegetation
  !! over bare ground, and ISFC=4 for bare ground. If a problem is
  !! discovered, a flag is set to the index of the modelled area, and
  !! a call to errorHandler is performed with an error message. Checks for
  !! unphysical values of certain water balance variables are
  !! performed against an accuracy limit ACCLMT, currently set to
  !! \f$1x10^{-3} kg m^{-2}\f$ or \f$m^3 m^{-3}\f$. The overall water balance of the
  !! subarea is checked against an accuracy limit BALLMT, currently
  !! set to \f$1x10^{-1} kg m^{-2}\f$. (These values reflect expected roundoff
  !! errors associated with 32-bit computation.)
  !!
  if (ISFC == 1 .or. ISFC == 3) then
    IPTBAD = 0
    JPTBAD = 0
  end if
  KPTBAD = 0
  !>
  !! In loop 100, for canopy-covered subareas, the intercepted rain
  !! RAICAN and snow SNOCAN are checked to ensure that if they are
  !! negative, they are vanishingly small. A similar check is done for
  !! the runoff.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. ISAND(I,1) > - 4) then
      if (ISFC == 1 .or. ISFC == 3) then
        if (RAICAN(I) < ( - 1.0 * ACCLMT)) IPTBAD = I
        if (SNOCAN(I) < ( - 1.0 * ACCLMT)) JPTBAD = I
      end if
      if (RUNOFF(I) < ( - 1.0 * ACCLMT)) KPTBAD = I
    end if
  end do ! loop 100
  !
  if (ISFC == 1 .or. ISFC == 3) then
    if (IPTBAD /= 0) then
      write(6,6100) IPTBAD,JL,ISFC,RAICAN(IPTBAD)
6100  format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' RAICAN = ', &
                E13.5)
      call errorHandler('checkWaterBudget', - 1)
    end if
    if (JPTBAD /= 0) then
      write(6,6150) JPTBAD,JL,ISFC,SNOCAN(JPTBAD)
6150  format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' SNOCAN = ', &
                E13.5)
      call errorHandler('checkWaterBudget', - 2)
    end if
  end if
  if (KPTBAD /= 0) then
    write(6,6200) KPTBAD,JL,ISFC,RUNOFF(KPTBAD)
6200 format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' RUNOFF = ', &
            E13.5)
    call errorHandler('checkWaterBudget', - 3)
  end if
  !
  IPTBDI = 0
  JPTBDI = 0
  KPTBDI = 0
  LPTBDI = 0
  !>
  !! In the 150 loop, for all areas that are not continental ice
  !! sheets (ISAND=-4), the liquid water content in each soil layer is
  !! checked to ensure that it is not larger than the pore volume and
  !! that it is not smaller than the minimum liquid water content
  !! (except for rock layers). The ice content is similarly checked to
  !! ensure that the sum of it, converted to an equivalent liquid
  !! water content, plus the minimum water content, is not greater
  !! than the pore volume (except for rock layers). It is also checked
  !! to ensure that if it is negative, it is vanishingly small.
  !!
  do J = 1,IG ! loop 150
    do I = IL1,IL2
      if (FI(I) > 0. .and. ISAND(I,1) > - 4) then
        if ((THLIQ(I,J) - THPOR(I,J)) > ACCLMT) then
          write(6,6009) I,J,N,THLIQ(I,J),THPOR(I,J),DELZW(I,J), &
                         THICE(I,J),THLIQI(I,J),THICEI(I,J), &
                         ZPOND(I),ZPONDI(I),ISAND(I,J)
6009      format(2X,3I6,8F16.8,I6)
          do K = 1,IG ! loop 145
            write(6,6008) K,THLIQ(I,K),THLIQI(I,K),THICE(I,K), &
                          THICEI(I,K),DELZW(I,K),THPOR(I,K),ISAND(I,K)
6008        format(2X,I6,6F14.8,I6)
          end do ! loop 145
          IPTBDI = I
          IPTBDJ = J
        end if
        if (THLIQ(I,J) < (THLMIN(I,J) - ACCLMT) .and. &
            ISAND(I,J) /= - 3) then
          JPTBDI = I
          JPTBDJ = J
        end if
        if ((THICE(I,J) * RHOICE / RHOW - THPOR(I,J) + THLMIN(I,J)) &
            > ACCLMT .and. ISAND(I,J) /= - 3) then
          KPTBDI = I
          KPTBDJ = J
        end if
        if (THICE(I,J) < - 1. * ACCLMT) then
          LPTBDI = I
          LPTBDJ = J
        end if
      end if
    end do
  end do ! loop 150
  !
  if (IPTBDI /= 0) then
    write(6,6250) IPTBDI,JL,ISFC,N,THLIQ(IPTBDI,IPTBDJ), &
                   THPOR(IPTBDI,IPTBDJ),IPTBDJ
6250 format('0AT (I,JL) = (',I6,',',I6,'),ISFC = ',I2,' STEP = ',I8, &
             ' THLIQ = ',E13.5,' THPOR = ',E13.5,' FOR J = ',I2)
    call errorHandler('checkWaterBudget', - 4)
  end if
  if (JPTBDI /= 0) then
    write(6,6300) JPTBDI,JL,ISFC,THLIQ(JPTBDI,JPTBDJ),JPTBDJ
6300 format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' THLIQ = ', &
             E13.5,' FOR J = ',I2)
    call errorHandler('checkWaterBudget', - 5)
  end if
  if (KPTBDI /= 0) then
    write(6,6350) KPTBDI,JL,ISFC,THICE(KPTBDI,KPTBDJ), &
                   THPOR(KPTBDI,KPTBDJ),KPTBDJ
6350 format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' THICE = ', &
             E13.5,' THPOR = ',E13.5,' FOR J = ',I2)
    write(6,6460) PCPR(KPTBDI) * DELT,EVAP(KPTBDI) * RHOW * DELT, &
         RUNOFF(KPTBDI) * RHOW,WLOST(KPTBDI), &
         RAICNI(KPTBDI) - RAICAN(KPTBDI),SNOCNI(KPTBDI) - &
         SNOCAN(KPTBDI),(ZPOND(KPTBDI) - ZPONDI(KPTBDI)) * RHOW
    write(6,6460) RAICAN(KPTBDI),RAICNI(KPTBDI), &
          SNOCAN(KPTBDI),SNOCNI(KPTBDI),ZPOND(KPTBDI), &
          ZPONDI(KPTBDI)
    write(6,6460) ZSNOW(KPTBDI) * RHOSNO(KPTBDI), &
         SNOFAC * SNOWI(KPTBDI),WSNOW(KPTBDI),WSNOWI(KPTBDI), &
         SNOFAC * SNOWI(KPTBDI) - ZSNOW(KPTBDI) * RHOSNO(KPTBDI), &
         WSNOWI(KPTBDI) - WSNOW(KPTBDI)
    write(6,6460) ZSNOW(KPTBDI),RHOSNO(KPTBDI),SNOFAC, &
         SNOWI(KPTBDI)
    do J = 1,IG ! loop 250
      write(6,6460) &
         THLIQ(KPTBDI,J),THLIQI(KPTBDI,J), &
         THICE(KPTBDI,J),THICEI(KPTBDI,J), &
         DELZW(KPTBDI,J),THPOR(KPTBDI,J), &
         (THLIQ(KPTBDI,J) - THLIQI(KPTBDI,J)) * RHOW * DELZW(KPTBDI,J), &
         (THICE(KPTBDI,J) - THICEI(KPTBDI,J)) * RHOICE * DELZW(KPTBDI,J)
    end do ! loop 250
    write(6,6470) FCS(KPTBDI),FGS(KPTBDI)
    call errorHandler('checkWaterBudget', - 6)
  end if
  if (LPTBDI /= 0) then
    write(6,6400) LPTBDI,JL,ISFC,THICE(LPTBDI,LPTBDJ),LPTBDJ
6400 format('0AT (I,JL) = (',I3,',',I3,'),ISFC = ',I2,' THICE = ', &
             E13.5,' FOR J = ',I2)
    call errorHandler('checkWaterBudget', - 7)
  end if
  !
  IPTBAD = 0
  if (ISFC == 1 .or. ISFC == 3) then
    CANFAC = 1.0
  else
    CANFAC = 0.0
  end if
  !>
  !! Finally, in loop 300, the overall water balance BAL is calculated
  !! and compared to BALLMT. BAL is evaluated as the residual of the
  !! precipitation, the evaporation, the runoff, the water loss term
  !! WLOST, the change in canopy intercepted liquid and frozen water
  !! (for vegetation-covered areas), the change in surface ponded
  !! water, the change in snow pack and snow liquid water content (for
  !! snow-covered areas), and the changes in the soil layer liquid and
  !! frozen water contents. If the absolute value of BAL is greater
  !! than BALLMT, a flag is set, all of the terms entering BAL are
  !! printed out, and a call to errorHandler is performed.
  !!
  do I = IL1,IL2 ! loop 300
    if (FI(I) > 0. .and. ZSNOW(I) > 0.) XSNOW(I) = 1.0
    if (FI(I) > 0. .and. ISAND(I,1) > - 4) then
      if (ISFC == 1 .or. ISFC == 2) then
        SNOFAC = 1.0 / (FCS(I) + FGS(I))
      else
        SNOFAC = 0.0
      end if
      BAL(I) = PCPR(I) * DELT - &    ! precip
               EVAP(I) * RHOW * DELT - &    ! evap
               RUNOFF(I) * RHOW + WLOST(I) - & ! runoff + wlost
               CANFAC * (RAICAN(I) - RAICNI(I) + SNOCAN(I) - SNOCNI(I)) - & ! canopy snow and water
               (ZPOND(I) - ZPONDI(I)) * RHOW - & ! ponded
               ZSNOW(I) * RHOSNO(I) + SNOFAC * SNOWI(I) - & ! snow
               WSNOW(I) + WSNOWI(I)

      do J = 1,IG ! loop 275
        BAL(I) = BAL(I) - &
                 (THLIQ(I,J) - THLIQI(I,J)) * RHOW * DELZW(I,J) - &   ! change in soil liquid content
                 (THICE(I,J) - THICEI(I,J)) * RHOICE * DELZW(I,J)   ! change in soil ice content
      end do ! loop 275
      if (ABS(BAL(I)) > BALLMT) then
        IPTBAD = I
      end if
    end if
  end do ! loop 300

  if (IPTBAD /= 0) then
    write(6,6450) IPTBAD,JL,N,ISFC,BAL(IPTBAD)
    write(6,6460) PCPR(IPTBAD) * DELT,EVAP(IPTBAD) * RHOW * DELT, &
         RUNOFF(IPTBAD) * RHOW,WLOST(IPTBAD), &
         RAICNI(IPTBAD) - RAICAN(IPTBAD),SNOCNI(IPTBAD) - &
         SNOCAN(IPTBAD),(ZPOND(IPTBAD) - ZPONDI(IPTBAD)) * RHOW
    write(6,6460) RAICAN(IPTBAD),RAICNI(IPTBAD), &
          SNOCAN(IPTBAD),SNOCNI(IPTBAD),ZPOND(IPTBAD), &
          ZPONDI(IPTBAD)
    write(6,6460) ZSNOW(IPTBAD) * RHOSNO(IPTBAD), &
         SNOFAC * SNOWI(IPTBAD),WSNOW(IPTBAD),WSNOWI(IPTBAD), &
         SNOFAC * SNOWI(IPTBAD) - ZSNOW(IPTBAD) * RHOSNO(IPTBAD), &
         WSNOWI(IPTBAD) - WSNOW(IPTBAD)
    write(6,6460) ZSNOW(IPTBAD),RHOSNO(IPTBAD),SNOFAC, &
         SNOWI(IPTBAD)
    do J = 1,IG ! loop 350
      write(6,6460) &
         THLIQ(IPTBAD,J),THLIQI(IPTBAD,J), &
         THICE(IPTBAD,J),THICEI(IPTBAD,J), &
         DELZW(IPTBAD,J),THPOR(IPTBAD,J), &
         (THLIQ(IPTBAD,J) - THLIQI(IPTBAD,J)) * RHOW * DELZW(IPTBAD,J), &
         (THICE(IPTBAD,J) - THICEI(IPTBAD,J)) * RHOICE * DELZW(IPTBAD,J)
    end do ! loop 350
    write(6,6470) FCS(IPTBAD),FGS(IPTBAD)
6450 format('0AT (I,JL) = (',I8,',',I8,'),TIME = ',I8,' ISFC = ',I2, &
         ' BAL = ',E13.5)
6460 format(2X,8F15.8)
6470 format(2X,4E20.6)
    call errorHandler('checkWaterBudget', - 8)
  end if

  return
end subroutine checkWaterBudget
