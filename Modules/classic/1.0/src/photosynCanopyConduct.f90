!> \file
!> Net Photosynthesis and canopy conductance
!> @author V. Arora, J. Melton, M. Lazare
subroutine photosynCanopyConduct (AILCG, FCANC, TCAN, CO2CONC, PRESSG, FC, & !In ! Formerly PHTSYN3
                                  CFLUX, QA, QSWV, IC, THLIQ, ISAND, & !In
                                  TA, RMAT, COSZS, XDIFFUS, ILG, & !In
                                  IL1, IL2, IG, ICC, ISNOW, SLAI, & !In
                                  THFC, THLW, FCANCMX, L2MAX, NOL2PFTS, & !In
                                  CO2I1, CO2I2, & ! In/Out
                                  RC, AN_VEG, RML_VEG, DAYL, DAYL_MAX) ! Out

  !     HISTORY:
  !
  !     * JAN 5/15  - J. Melton      Use TCAN again. We instituted a PAI minimum for trees (grass and crops already had it) of
  !                                  of 1. This removed the wild TCAN values that could occur.
  !     * JUL 22/15 - J. Melton      SM_FUNC2 was not correctly set up for soil layers > 3 (Noted by Yuanqiao Wu). Now fixed.
  !     * FEB 27/15 - J. Melton      Increase REQUITR to 10, it was not converging correctly at 4. Also rename
  !                                  WILTSM and FIELDSM to THLW and THFC, respectively, for consistency with CLASS.
  !     * SEP 15/14 - J. Melton      Since SN is converted to INT, I made the assignment explicitly INT, rather
  !     *                            than a real :: that gets cast as INT later.
  !     * APR 11/14 - J. Melton      When at wilting point SM_FUNC should be 0. Was incorrectly set to 0.01
  !     * JAN 15/14 - J. Melton      Fixed bug in SM_FUNC calculation that prevented
  !     *                            higher number PFTs from feeling water stress.
  !     *                            Also have new parameter values in Vmax, SN, and RMLcoef
  !     * AUG 16/12 - J. MELTON      Fixed GB constraint, removed incorrect
  !     *                            scaling of RML, made reals being treated
  !     *                            as ints more explicit, changed declaration
  !     *                            of usebb to real. Also changed use of TA to
  !     *                            TCAN for RH and GB calculations to ensure
  !     *                            consistency with rest of subroutine.
  !     * NOV 15, 2011 - M.LAZARE.   New version for gcm15j+ to support
  !     *                            class_v3.5:
  !     *                            - SAND and CLAY no longer passed in.
  !     *                            - ISAND passed in so automatic array
  !     *                              removed.
  !     *                            - FIELDSM and WILTSM now passed in
  !     *                              (calculated elsewhere in model)
  !     *                              so no need for other soil propery
  !     *                              automatic arrays.
  !     * NOV 12, 2010 - M.LAZARE.   Revised (cosmetic) to put if condition
  !     *                            on execution of loop 460 based on "fc"
  !     *                            (known in tsolvc6 as "fct") since
  !     *                            "qswv" (known in tsolvc6 as "qswnvc")
  !     *                            is only defined for fct>0. This modification
  !     *                            has been tested to give bit-for-bit the same
  !     *                            answer. The model runs ok without this for
  !     *                            the usual "noxlfqinitauto=on" but otherwise
  !     *                            crashes without the change.
  !     * APR 23, 2010 - V.ARORA/    PREVIOUS VERSION PHTSYN2 FOR GCM15I:
  !     *                M.LAZARE.   - GAMMA_W NOW 0.25 INSTEAD OF 0.425.
  !     *                            - QA IS NOW INPUT WITH RH AS AN
  !     *                              ALLOCATABLE ARRAY, INSTEAD OF RH
  !     *                              AS INPUT AND EA, EASAT AS ALLOCATABLE
  !     *                              ARRAYS. THE LATTER TWO DON'T HAVE
  !     *                              TO BE ARRAYS TO CALCULATE VPD AND
  !     *                              DOING IT THIS WAY IS MORE CONSISTENT
  !     *                              WITH THE REST OF THE MODEL.
  !     * 7 SEP.  2001 - PROGRAM TO CALCULATE STOMATAL CONDUCTANCE, TO BE
  !     * V. ARORA       USED BY CLASS, BY ESTIMATING PHOTOSYNTHESIS AND THEN
  !     *                RELATING PHOTOSYNTHESIS TO STOMATAL CONDUCTANCE
  !     *                PREVIOUS VERSION PHTSYN).
  !
  !     ------------------------------------------------------------------
  !
  use classicParams, only : KN, TUP, TLOW, alpha_phtsyn, omega_phtsyn, &
                            ISC4, MM, BB, VPD0, SN, SMSCALE, VMAX, REQITER, &
                            CO2IMAX, BETA1, BETA2, INICO2I, CHI, RMLCOEFF, &
                            GAMMA_W, GAMMA_M, TFREZ, ZERO, STD_PRESS
  use generalutilsGAT,  only : calcEsat
  
  implicit none
  !
  integer, DIMENSION(:,:), ALLOCATABLE  :: USESLAI
  integer, DIMENSION(:), ALLOCATABLE    :: SORT

  real, DIMENSION(:), ALLOCATABLE       :: FC_TEST, SIGMA, TGAMMA, &
                                           KC, KO, IPAR, GB, RH, VPD, O2_CONC, CO2A, USEBB

  real, DIMENSION(:,:), ALLOCATABLE     :: USEAILCG, SM_FUNC, &
                                           AVE_SM_FUNC, VMAXC, JE3, SM_FUNC2, TOT_RMAT, &
                                           VMUNS1, VMUNS2, VMUNS3, VMUNS, VM, CO2I, PREV_CO2I, &
                                           FPAR, JC, JC1, JC2, JC3, JE, JE1, JE2, JS, A_VEG, &
                                           RC_VEG, GCTU, GCMIN, GCMAX, VPD_TERM, CO2LS, GC
  !    -----------------------------------------------------------------
  !                 VARIABLES USED ONLY FOR THE TWO LEAF MODEL
  !
  real, DIMENSION(:), ALLOCATABLE :: IPAR_SUN, IPAR_SHA

  real, DIMENSION(:,:), ALLOCATABLE     :: GDIR, KB, FPAR_SUN, &
                                           FPAR_SHA, VMAXC_SUN, VMAXC_SHA, VMUNS1_SUN, VMUNS1_SHA, &
                                           VMUNS_SUN, VMUNS_SHA, VM_SUN, VM_SHA, CO2I_SUN, PREV_CO2I_SUN, &
                                           CO2I_SHA, PREV_CO2I_SHA, JC1_SUN, JC1_SHA, JC3_SUN, JC3_SHA, &
                                           JC_SUN, JC_SHA, JE1_SUN, JE1_SHA, &
                                           JE2_SUN, JE2_SHA, JE_SUN, JE_SHA, JS_SUN, JS_SHA, &
                                           A_VEG_SUN, A_VEG_SHA, RML_SUN, RML_SHA, AN_SUN, AN_SHA, &
                                           CO2LS_SUN, CO2LS_SHA, AILCG_SUN, AILCG_SHA, GC_SUN, GC_SHA, &
                                           GCTU_SUN, GCTU_SHA

  integer, intent(in) :: ILG      !< NO. OF GRID CELLS IN LATITUDE CIRCLE
  integer, intent(in) :: IC       !< NO. OF CLASS VEGETATION TYPES, 4
  integer, intent(in) :: IL1      !< IL1=1
  integer, intent(in) :: IL2      !< IL2=ILG
  integer, intent(in) :: IG       !< NO. OF SOIL LAYERS, 3
  integer, intent(in) :: ICC      !< NO. OF CTEM's PFTs, CURRENTLY 9
  integer, intent(in) :: ISNOW    !< integer, intent(in) (0 or 1) TELLING IF PHTSYN IS TO BE RUN OVER CANOPY OVER SNOW OR CANOPY OVER GROUND SUBAREA
  integer, intent(in) :: NOL2PFTS(IC) !<
  integer, intent(in) :: L2MAX    !< MAX. NUMBER OF LEVEL 2 PFTs
  integer :: PS_COUP, LEAFOPT, K1, K2, ICOUNT, IT_COUNT, I, J, M, N  !<
  !
  real, intent(in) :: FCANC(ILG,ICC) !< FRACTIONAL COVERAGE OF CTEM's 9 PFTs
  real, intent(in) :: AILCG(ILG,ICC) !< GREEN LEAF AREA INDEX FOR USE BY PHOTOSYNTHESIS, \f$M^2/M^2\f$
  real, intent(in) :: TCAN(ILG)      !< CANOPY TEMPERATURE, KELVIN
  real, intent(in) :: FC(ILG)        !< SUM OF ALL FCANC OVER A GIVEN SUB-AREA
  real, intent(in) :: CFLUX(ILG)     !< AERODYNAMIC CONDUCTANCE, M/S
  real, intent(in) :: SLAI(ILG,ICC)  !< SCREEN LEVEL HUMIDITY IN KG/KG - STORAGE LAI. THIS LAI IS USED FOR PHTSYN EVEN IF ACTUAL LAI IS ZERO. ESTIMATE OF NET PHOTOSYNTHESIS BASED ON SLAI IS USED FOR INITIATING LEAF ONSET. SEE PHENOLGY SUBROUTINE FOR MORE DETAILS.
  real, intent(in) :: QA(ILG)        !<
  !
  real, intent(in) :: CO2CONC(ILG) !< ATMOS. \f$CO_2\f$ IN PPM, AND THEN CONVERT IT TO PARTIAL PRESSURE, PASCALS, CO2A, FOR USE IN THIS SUBROUTINE
  real :: Q10_FUNCN        !<
  real :: Q10_FUNC         !<
  real, intent(in) :: PRESSG(ILG)      !< ATMOS. PRESSURE, PASCALS
  real, intent(out) :: RML_VEG(ILG,ICC) !< LEAF RESPIRATION RATE, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) FOR EACH PFT
  real, intent(out) :: AN_VEG(ILG,ICC)  !< NET PHOTOSYNTHESIS RATE, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) FOR EACH PFT
  real, intent(in) :: QSWV(ILG)        !< ABSORBED VISIBLE PART OF SHORTWAVE RADIATION, \f$W/M^2\f$
  real, intent(in) :: TA(ILG)          !< AIR TEMPERATURE IN KELVINS
  real, intent(in) :: RMAT(ILG, ICC,IG) !< FRACTION OF ROOTS IN EACH LAYER (grid cell, vegetation, layer)
  real, intent(inout) :: CO2I1(ILG,ICC)   !< INTERCELLULAR \f$CO_2\f$ CONCENTRATION FROM THE PREVIOUS TIME STEP WHICH GETS UPDATED FOR THE SINGLE LEAF OR THE SUNLIT PART OF THE TWO LEAF MODEL
  real, intent(inout) :: CO2I2(ILG,ICC)   !< INTERCELLULAR \f$CO_2\f$ CONCENTRATION FOR THE SHADED PART OF THE TWO LEAF MODEL FROM THE PREVIOUS TIME STEP
  real :: CA               !<
  real :: CB               !<
  real, intent(in) :: THLIQ(ILG,IG)    !< LIQUID MOIS. CONTENT OF 3 SOIL LAYERS
  real, intent(in) :: THFC(ILG,IG)     !< SOIL FIELD CAPACITY.
  real, intent(in) :: THLW(ILG,IG)     !< SOIL WILT CAPACITY.
  real, intent(in) :: FCANCMX(ILG,ICC) !< MAX. FRACTIONAL COVERAGES OF CTEM's 8 PFTs. THIS IS DIFFERENT FROM FCANC AND FCANCS (WHICH MAY VARY WITH SNOW DEPTH).
  !< FCANCMX DOESN'T CHANGE, UNLESS OF COURSE ITS CHANGED BY LAND USE CHANGE OR DYNAMIC VEGETATION.
  real :: Q10_FUNCD        !<
  !
  real :: Q10
  real, intent(out) :: RC(ILG)      !< GRID-AVERAGED STOMATAL RESISTANCE, S/M
  real, intent(in) :: COSZS(ILG)   !< COS OF ZENITH ANGLE
  real, intent(in) :: XDIFFUS(ILG) !< FRACTION OF DIFFUSED PAR
  !
  real :: TEMP_B
  real :: TEMP_C
  real :: TEMP_R
  real :: TEMP_Q1
  real :: TEMP_Q2
  real :: TEMP_JP
  real ::  TEMP_AN
  !
  integer, intent(in) :: ISAND(ILG,IG) !< SAND INDEX.
  real, intent(in) :: DAYL_MAX(ILG)      !< MAXIMUM DAYLENGTH FOR THAT LOCATION
  real, intent(in) :: DAYL(ILG)          !< DAYLENGTH FOR THAT LOCATION
  real :: use_vmax
  !
  !     FOR LIMITING CO2 UPTAKE
  !
  real :: DELTA_CO2(ILG)
  real :: N_EFFECT(ILG)
  logical :: SMOOTH, MIN2, MIN3
  !
  !
  real :: TEMP_PHI1, TEMP_PHI2, EA, EASAT, T_TEMP(ILG)
  !     -----------------------------------------------------------------
  !     DECIDE HERE IF WE WANT TO USE SINGLE LEAF OR TWO-LEAF MODEL
  !     CHOOSE 1 FOR SINGLE-LEAF MODEL, AND 2 FOR TWO-LEAF MODEL
  !
  LEAFOPT = 1
  !
  !     DECIDE IF WE WANT TO USE BWB (1) OR LEUNING TYPE (2) PHOTOSYNTHESIS
  !     STOMATAL CONDUCTANCE COUPLING
  !
  PS_COUP = 2
  !
  !     DECIDE IF WE WANT TO ESTIMATE PHOTOSYNTHETIC RATE AS A SMOOTHED
  !     AVERAGE OF THE 3 LIMITING RATES, OR AS A MINIMUM OF THE 3 LIMITING
  !     RATES, OR AS A MINIMUM OF THE RUBSICO AND LIGHT LIMITING RATES
  !
  SMOOTH = .true.
  MIN2 = .false.
  MIN3 = .false.
  !
  !     --------------------------------------------------------------
  !
  !
  ALLOCATE(USESLAI(ILG,ICC))
  ALLOCATE(SORT(ICC))
  ALLOCATE(USEBB(ICC))

  ALLOCATE(FC_TEST(ILG))
  ALLOCATE(USEAILCG(ILG,ICC))
  ALLOCATE(SM_FUNC(ILG,IG))
  ALLOCATE(SM_FUNC2(ILG,IG))
  ALLOCATE(AVE_SM_FUNC(ILG,ICC))
  ALLOCATE(TOT_RMAT(ILG,ICC))

  ALLOCATE(VMAXC(ILG,ICC))
  ALLOCATE(JE3(ILG,ICC))
  ALLOCATE(VMUNS1(ILG,ICC),VMUNS2(ILG,ICC),VMUNS3(ILG,ICC))
  ALLOCATE(VMUNS(ILG,ICC))
  ALLOCATE(VM(ILG,ICC))
  ALLOCATE(SIGMA(ILG))
  ALLOCATE(TGAMMA(ILG))
  ALLOCATE(KC(ILG))
  ALLOCATE(KO(ILG))
  ALLOCATE(CO2I(ILG,ICC),PREV_CO2I(ILG,ICC))
  ALLOCATE(FPAR(ILG,ICC))
  ALLOCATE(JC(ILG,ICC))
  ALLOCATE(JC1(ILG,ICC),JC2(ILG,ICC),JC3(ILG,ICC))
  ALLOCATE(JE(ILG,ICC),JE1(ILG,ICC),JE2(ILG,ICC))
  ALLOCATE(IPAR(ILG))
  ALLOCATE(JS(ILG,ICC))
  ALLOCATE(A_VEG(ILG,ICC))
  ALLOCATE(GB(ILG))
  ALLOCATE(RC_VEG(ILG,ICC))
  ALLOCATE(GCTU(ILG,ICC))
  ALLOCATE(GCMIN(ILG,ICC))
  ALLOCATE(GCMAX(ILG,ICC))
  ALLOCATE(RH(ILG))
  ALLOCATE(VPD(ILG),VPD_TERM(ILG,ICC))

  ALLOCATE(CO2LS(ILG,ICC))
  ALLOCATE(GC(ILG,ICC))
  ALLOCATE(O2_CONC(ILG))
  ALLOCATE(CO2A(ILG))

  ALLOCATE(GDIR(ILG,ICC),KB(ILG,ICC))
  ALLOCATE(FPAR_SUN(ILG,ICC),FPAR_SHA(ILG,ICC))
  ALLOCATE(VMAXC_SUN(ILG,ICC),VMAXC_SHA(ILG,ICC))
  ALLOCATE(VMUNS1_SUN(ILG,ICC),VMUNS1_SHA(ILG,ICC))
  ALLOCATE(VMUNS_SUN(ILG,ICC),VMUNS_SHA(ILG,ICC))
  ALLOCATE(VM_SUN(ILG,ICC),VM_SHA(ILG,ICC))
  ALLOCATE(CO2I_SUN(ILG,ICC),PREV_CO2I_SUN(ILG,ICC))
  ALLOCATE(CO2I_SHA(ILG,ICC),PREV_CO2I_SHA(ILG,ICC))
  ALLOCATE(JC1_SUN(ILG,ICC),JC1_SHA(ILG,ICC))
  ALLOCATE(JC3_SUN(ILG,ICC),JC3_SHA(ILG,ICC))
  ALLOCATE(JC_SUN(ILG,ICC),JC_SHA(ILG,ICC))
  ALLOCATE(IPAR_SUN(ILG),IPAR_SHA(ILG))
  ALLOCATE(JE1_SUN(ILG,ICC),JE1_SHA(ILG,ICC))
  ALLOCATE(JE2_SUN(ILG,ICC),JE2_SHA(ILG,ICC))
  ALLOCATE(JE_SUN(ILG,ICC),JE_SHA(ILG,ICC))
  ALLOCATE(JS_SUN(ILG,ICC),JS_SHA(ILG,ICC))
  ALLOCATE(A_VEG_SUN(ILG,ICC),A_VEG_SHA(ILG,ICC))
  ALLOCATE(RML_SUN(ILG,ICC),RML_SHA(ILG,ICC))
  ALLOCATE(AN_SUN(ILG,ICC),AN_SHA(ILG,ICC))
  ALLOCATE(CO2LS_SUN(ILG,ICC),CO2LS_SHA(ILG,ICC))
  ALLOCATE(AILCG_SUN(ILG,ICC),AILCG_SHA(ILG,ICC))
  ALLOCATE(GC_SUN(ILG,ICC),GC_SHA(ILG,ICC))
  ALLOCATE(GCTU_SUN(ILG,ICC),GCTU_SHA(ILG,ICC))

  !
  !     --------------------------------------------------------------
  !
  !     INITIALIZATION
  !
  if (LEAFOPT == 1) then
    !       INITIALIZE REQUIRED ARRAYS TO ZERO FOR SINGLE LEAF MODEL
    do I = IL1,IL2
      CO2A(I) = 0.0
      IPAR(I) = 0.0
      SIGMA(I) = 0.0
      TGAMMA(I) = 0.0
      KC(I) = 0.0
      KO(I) = 0.0
      GB(I) = 0.0
      RC(I) = 0.0
      FC_TEST(I) = 0.0
    end do ! loop 100
    !
    do J = 1,ICC
      USEBB(J) = 0.0
      do I = IL1,IL2
        FPAR(I,J) = 0.0
        VMAXC(I,J) = 0.0
        VMUNS1(I,J) = 0.0
        VMUNS2(I,J) = 0.0
        VMUNS3(I,J) = 0.0
        VMUNS(I,J) = 0.0
        AVE_SM_FUNC(I,J) = 0.0
        TOT_RMAT(I,J) = 0.0
        VM(I,J) = 0.0
        JC1(I,J) = 0.0
        JC2(I,J) = 0.0
        JC3(I,J) = 0.0
        JC(I,J) = 0.0
        JE(I,J) = 0.0
        JE1(I,J) = 0.0
        JE2(I,J) = 0.0
        JS(I,J) = 0.0
        A_VEG(I,J) = 0.0
        RML_VEG(I,J) = 0.0
        AN_VEG(I,J) = 0.0
        CO2LS(I,J) = 0.0
        GC(I,J) = 0.0
        GCTU(I,J) = 0.0
        RC_VEG(I,J) = 5000.0
        USESLAI(I,J) = 0
        USEAILCG(I,J) = 0.0
      end do ! loop 201
    end do ! loop 200
    !
  else if (LEAFOPT == 2) then
    !       INITIALIZE ARRAYS FOR THE TWO LEAF MODEL
    do I = IL1,IL2
      CO2A(I) = 0.0
      IPAR(I) = 0.0
      SIGMA(I) = 0.0
      TGAMMA(I) = 0.0
      KC(I) = 0.0
      KO(I) = 0.0
      GB(I) = 0.0
      RC(I) = 0.0
      IPAR_SUN(I) = 0.0
      IPAR_SHA(I) = 0.0
      FC_TEST(I) = 0.0
    end do ! loop 210
    !
    do J = 1,ICC
      USEBB(J) = 0.0
      do I = IL1,IL2
        GDIR(I,J) = 0.0
        KB(I,J) = 0.0
        AILCG_SUN(I,J) = 0.0
        AILCG_SHA(I,J) = 0.0
        FPAR_SUN(I,J) = 0.0
        FPAR_SHA(I,J) = 0.0
        VMAXC_SUN(I,J) = 0.0
        VMAXC_SHA(I,J) = 0.0
        VMUNS1_SUN(I,J) = 0.0
        VMUNS1_SHA(I,J) = 0.0
        VMUNS_SUN(I,J) = 0.0
        VMUNS_SHA(I,J) = 0.0
        AVE_SM_FUNC(I,J) = 0.0
        TOT_RMAT(I,J) = 0.0
        VM_SUN(I,J) = 0.0
        VM_SHA(I,J) = 0.0
        JC1_SUN(I,J) = 0.0
        JC1_SHA(I,J) = 0.0
        JC3_SUN(I,J) = 0.0
        JC3_SHA(I,J) = 0.0
        JC_SUN(I,J) = 0.0
        JC_SHA(I,J) = 0.0
        JE_SUN(I,J) = 0.0
        JE_SHA(I,J) = 0.0
        JE1_SUN(I,J) = 0.0
        JE1_SHA(I,J) = 0.0
        JE2_SUN(I,J) = 0.0
        JE2_SHA(I,J) = 0.0
        JS_SUN(I,J) = 0.0
        JS_SHA(I,J) = 0.0
        A_VEG_SUN(I,J) = 0.0
        A_VEG_SHA(I,J) = 0.0
        RML_SUN(I,J) = 0.0
        RML_SHA(I,J) = 0.0
        AN_SUN(I,J) = 0.0
        AN_SHA(I,J) = 0.0
        CO2LS_SUN(I,J) = 0.0
        CO2LS_SHA(I,J) = 0.0
        GC_SUN(I,J) = 0.0
        GC_SHA(I,J) = 0.0
        GCTU_SUN(I,J) = 0.0
        GCTU_SHA(I,J) = 0.0
        RC_VEG(I,J) = 5000.0
        AN_VEG(I,J) = 0.0
        RML_VEG(I,J) = 0.0
        USESLAI(I,J) = 0
        USEAILCG(I,J) = 0.0
      end do ! loop 220
    end do ! loop 230
  end if
  !
  !     FOLLOWING VARIABLES AND CONSTANTS ARE COMMON TO BOTH SINGLE AND TWO-LEAF MODEL
  do J = 1,IG
    do I = IL1,IL2
      SM_FUNC(I,J) = 0.0
      SM_FUNC2(I,J) = 0.0
    end do ! loop 250
  end do ! loop 240
  !
  do J = 1,ICC
    do I = IL1,IL2
      GCMIN(I,J) = 0.0
      GCMAX(I,J) = 0.0
    end do ! loop 270
  end do ! loop 260
  !
  !     GENERATE THE SORT INDEX FOR CORRESPONDENCE BETWEEN 9 PFTs AND THE
  !     12 VALUES IN THE PARAMETER VECTORS
  !
  ICOUNT = 0
  do J = 1,IC
    do M = 1,NOL2PFTS(J)
      N = (J - 1) * L2MAX + M
      ICOUNT = ICOUNT + 1
      SORT(ICOUNT) = N
    end do ! loop 281
  end do ! loop 280
  !
  !     INITIALIZATION ENDS
  !
  !     -------------------------------------------------------------------
  !>
  !! IF LAI IS LESS THAN SLAI THAN WE USE STORAGE LAI TO PHOTOSYNTHESIZE. HOWEVER, WE DO NOT USE THE STOMATAL
  !! RESISTANCE ESTIMATED IN THIS CASE, BECAUSE STORAGE LAI IS AN IMAGINARY LAI, AND WE SET STOMATAL RESISTANCE
  !! TO ITS MAX. NOTE THAT THE CONCEPT OF STORAGE/IMAGINARY LAI IS USED FOR PHENOLOGY PURPOSES AND THIS
  !! IMAGINARY LAI ACTS AS MODEL SEEDS.
  !!
  do J = 1,ICC
    do I = IL1,IL2
      if (AILCG(I,J) < SLAI(I,J)) then
        USESLAI(I,J) = 1
        USEAILCG(I,J) = SLAI(I,J)
      else
        USEAILCG(I,J) = AILCG(I,J)
      end if
    end do ! loop 350
  end do ! loop 340
  !>
  !! SET MIN. AND MAX. VALUES FOR STOMATAL CONDUCTANCE. WE MAKE SURE THAT MAX. STOMATAL RESISTANCE IS AROUND
  !! 5000 S/M AND MIN. STOMATAL RESISTANCE IS 51 S/M.
  !!
  do J = 1,ICC
    do I = IL1,IL2
      GCMIN(I,J) = 0.0002 * (TFREZ / TCAN(I)) * (1. / 0.0224) * &
                   (PRESSG(I) / STD_PRESS)
      !
      if (LEAFOPT == 1) then
        !           GCMAX(I,J)=0.0196 * (TFREZ/TCAN(I)) * (1./0.0224) *
        GCMAX(I,J) = 0.1    * (TFREZ / TCAN(I)) * (1. / 0.0224) * &
                     (PRESSG(I) / STD_PRESS)
      else if (LEAFOPT == 2) then
        !           GCMAX(I,J)=0.0196 * (TFREZ/TCAN(I)) * (1./0.0224) *
        GCMAX(I,J) = 0.1    * (TFREZ / TCAN(I)) * (1. / 0.0224) * &
                     (PRESSG(I) / STD_PRESS) * 0.5
      end if
      !
    end do ! loop 370
  end do ! loop 360
  !>
  !! IF WE ARE USING LEUNING TYPE PHOTOSYNTHESIS-STOMATAL CONDUCTANCE COUPLING WE NEED VAPOR PRESSURE DEFICIT
  !! AS WELL. CALCULATE THIS FROM THE RH AND AIR TEMPERATURE WE HAVE. WE FIND E_SAT, E, AND VPD IN PASCALS.
  !!
  if (PS_COUP == 2) then
    do J = 1,ICC
      do I = IL1,IL2
        VPD_TERM(I,J) = 0.0
      end do ! loop 400
    end do ! loop 390
    !
    do I = IL1,IL2
      VPD(I) = 0.0
      EASAT  = calcEsat(TCAN(I))
      EA     = QA(I) * PRESSG(I) / (0.622 + 0.378 * QA(I))
      RH(I)  = EA / EASAT
      VPD(I) = EASAT - EA
      VPD(I) = MAX(0.0,VPD(I))

    end do ! loop 420
    !
    K1 = 0
    do J = 1,IC
      if (J == 1) then
        K1 = K1 + 1
      else
        K1 = K1 + NOL2PFTS(J - 1)
      end if
      K2 = K1 + NOL2PFTS(J) - 1
      do M = K1,K2
        do I = IL1,IL2
          VPD_TERM(I,M) = 1.0 / ( 1.0 + ( VPD(I) / VPD0(SORT(M)) ) )
        end do ! loop 450
      end do ! loop 445
    end do ! loop 440
  end if
  !>
  !! ESTIMATE PARTIAL PRESSURE OF \f$CO_2\f$ AND IPAR
  !!
  do I = IL1,IL2
    if (FC(I) > 0.) then
      !> CONVERT CO2CONC FROM PPM TO PASCALS
      CO2A(I) = CO2CONC(I) * PRESSG(I) / 1E+06
      !> CHANGE PAR FROM W/M^2 TO MOL/M^2.S
      IPAR(I) = QSWV(I) * 4.6 * 1E-06
      !>
      !> SUNLIT PART GETS BOTH DIRECT AND DIFFUSED, WHILE THE SHADED PART GETS ONLY DIFFUSED
      !>
      IPAR_SUN(I) = QSWV(I) * 4.6 * 1E-06
      IPAR_SHA(I) = QSWV(I) * 4.6 * 1E-06 * XDIFFUS(I)
    end if
  end do ! loop 460
  !
  K1 = 0
  do J = 1,IC
    if (J == 1) then
      K1 = K1 + 1
    else
      K1 = K1 + NOL2PFTS(J - 1)
    end if
    K2 = K1 + NOL2PFTS(J) - 1
    do M = K1,K2
      do I = IL1,IL2
        if (FCANC(I,M) > ZERO) then
          !>
          !> FOR TWO-LEAF MODEL FIND Kb AS A FUNCTION OF COSZS AND LEAF ANGLE DISTRIBUTION (VEGETATION DEPENDENT)
          !>
          if (LEAFOPT == 2) then
            if (COSZS(I) > 0.0) then
              !>
              !> MAKE SURE -0.4 < CHI < 0.6
              CHI(SORT(M)) = MIN (MAX (CHI(SORT(M)), - 0.4),0.6)
              !> MAKE VALUES CLOSE TO ZERO EQUAL TO 0.01
              if (ABS(CHI(SORT(M))) <= 0.01) CHI(SORT(M)) = 0.01
              !
              TEMP_PHI1 = 0.5 - 0.633 * &
                          CHI(SORT(M)) - 0.33 * CHI(SORT(M)) * CHI(SORT(M))
              TEMP_PHI2 = 0.877 * (1. - 2. * TEMP_PHI1)
              GDIR(I,M) = TEMP_PHI1 + TEMP_PHI2 * COSZS(I)
              KB(I,M) = (GDIR(I,M) / COSZS(I))
              KB(I,M) = KB(I,M) * (SQRT(1. - omega_phtsyn(SORT(M)) ))
              !>
              !> ALSO FIND SUNLIT AND SHADED LAI
              AILCG_SUN(I,M) = (1.0 / KB(I,M)) * &
                               (1.0 - EXP( - 1.0 * KB(I,M) * USEAILCG(I,M) ) )
              AILCG_SHA(I,M) = USEAILCG(I,M) - AILCG_SUN(I,M)
              !>
              !> FOLLOWING FEW LINES TO MAKE SURE THAT ALL LEAVES ARE SHADED WHEN XDIFFUS EQUALS 1. NOT DOING
              !! SO GIVES ERRATIC RESULTS WHEN TWO LEAF OPTION IS USED
              !!
              if (XDifFUS(I) > 0.99) then
                AILCG_SUN(I,M) = 0.0
                AILCG_SHA(I,M) = USEAILCG(I,M)
              end if
              !
            end if
          end if
          !>
          !> FIND FPAR - FACTOR FOR SCALING PHOTOSYNTHESIS TO CANOPY BASED ON ASSUMPTION THAT NITROGEN IS
          !! OPTIMALLY DISTRIBUTED. THE TWO-LEAF MODEL IS NOT THAT DIFFERENT FROM THE SINGLE-LEAF MODEL.
          !! ALL WE DO IS USE TWO SCALING FACTORS (I.E. SCALING FROM LEAF TO CANOPY) INSTEAD OF ONE, AND
          !! THUS PERFORM CALCULATIONS TWICE, AND IN THE END ADD CONDUCTANCE AND NET PHOTOSYNTHESIS FROM
          !! THE TWO LEAVES TO GET THE TOTAL.
          !!
          FPAR(I,M) = (1.0 / KN(SORT(M))) * (1.0 - EXP( - KN(SORT(M)) &
                      * USEAILCG(I,M)))

          if (LEAFOPT == 2) then
            FPAR_SUN(I,M) = ( 1.0 / (KN(SORT(M)) + KB(I,M)) ) * &
                            (1.0 - EXP( - 1. * (KN(SORT(M)) + KB(I,M)) * USEAILCG(I,M) ) )
            FPAR_SHA(I,M) = FPAR(I,M) - FPAR_SUN(I,M)
            !>
            !! IF ALL RADIATION IS DIFFUSED, THEN ALL LEAVES ARE SHADED, AND WE ADJUST FPARs ACCORDINGLY.
            !! WITHOUT THIS THE TWO LEAF MODELS MAY BEHAVE ERRATICALLY
            !!
            if (XDifFUS(I) > 0.99) then
              FPAR_SHA(I,M) = FPAR(I,M)
              FPAR_SUN(I,M) = 0.0
            end if
          end if
          !>
          !> FIND Vmax, canopy, THAT IS Vmax SCALED BY LAI FOR THE SINGLE LEAF MODEL
          !>
          !> ------------- Changing Vcmax seasonally -----------------------
          !!
          !! Based on \cite Bauerle2012-c29  and \cite Alton2017-pd
          !! there is good evidence for the Vcmax varying throughout the season for deciduous tree
          !! species. We are adopting a parameterization based upon their paper with some differences.
          !! We don't apply it to evergreens like they suggest. Their paper had only one evergreen species
          !! and other papers (\cite Miyazawa2006-so) don't seem to back that up. Grasses and crops are also
          !! not affected by the dayl. \cite Alton2017-pd seems to indicate that all PFTs except
          !! BDL-EVG tropical should vary intra-annually (see their figure 8).
          !!

          if ((m == 2 .or. m == 4)) then
            use_vmax = vmax(sort(m)) * (dayl(i) / dayl_max(i)) ** 2
          else ! evergreens,crops,and grasses
            use_vmax = vmax(sort(m))
          end if

          vmaxc(i,m) = use_vmax * fpar(i,m)

          if (LEAFOPT == 2) then
            !> The two leaf is assumed to be affect by the insolation seasonal cycle the
            !> same for each sun/shade leaf
            VMAXC_SUN(I,M) = use_vmax * FPAR_SUN(I,M)
            VMAXC_SHA(I,M) = use_vmax * FPAR_SHA(I,M)
          end if
          !>
          !> ------------- Changing Vcmax seasonally -----------------------///


          !> FIND Vm, unstressed (DUE TO WATER) BUT STRESSED DUE TO TEMPERATURE
          !>
          Q10 = 2.00
          Q10_FUNC = Q10 ** (0.1 * (TCAN(I) - 298.16))
          !
          if (COSZS(I) > 0.0) then
            if (LEAFOPT == 1) then
              VMUNS1(I,M) = VMAXC(I,M) * Q10_FUNC
            else if (LEAFOPT == 2) then
              VMUNS1_SUN(I,M) = VMAXC_SUN(I,M) * Q10_FUNC
              VMUNS1_SHA(I,M) = VMAXC_SHA(I,M) * Q10_FUNC
            end if
          end if
          !>
          !> ASSUMING THAT SUNLIT AND SHADED TEMPERATURES ARE SAME
          !>
          VMUNS2(I,M) = (1. + EXP(0.3 * (TCAN(I) - TUP(SORT(M)))))
          VMUNS3(I,M) = (1. + EXP(0.3 * (TLOW(SORT(M)) - TCAN(I))))
          !
          if (LEAFOPT == 1) then
            VMUNS(I,M) = VMUNS1(I,M) / &
                         (VMUNS2(I,M) * VMUNS3(I,M))
          else if (LEAFOPT == 2) then
            VMUNS_SUN(I,M) = VMUNS1_SUN(I,M) / &
                             (VMUNS2(I,M) * VMUNS3(I,M))
            VMUNS_SHA(I,M) = VMUNS1_SHA(I,M) / &
                             (VMUNS2(I,M) * VMUNS3(I,M))
          end if
          !
        end if
      end do ! loop 490
    end do ! loop 485
  end do ! loop 480
  !>
  !> CALCULATE SOIL MOIS STRESS TO ACCOUNT FOR REDUCTION IN PHOTOSYN
  !> DUE TO LOW SOIL MOISTURE, THREE STEPS HERE
  !> 1. FIND WILTING
  !> POINT AND FIELD CAPACITY SOIL MOIS. CONTENT FOR ALL THREE LAYERS.
  !> 2. USING THESE FIND THE SOIL MOISTURE STRESS TERM FOR ALL
  !> THREE LAYERS
  !> 3. AVERAGE THIS SOIL MOISTURE STRESS TERM
  !> OVER THE 3 LAYERS USING FRACTION OF ROOTS PRESENT IN EACH LAYER
  !> FOR EACH PFT.
  !!
  !! NOTE THAT WHILE SOIL MOISTURE IS UNIFORM OVER
  !> AN ENTIRE GCM GRID CELL, THE SOIL MOISTURE STRESS FOR EACH
  !> PFT IS NOT BECAUSE OF DIFFERENCES IN ROOT DISTRIBUTION.
  !>
  !> WILTING POINT CORRESPONDS TO MATRIC POTENTIAL OF 150 M
  !> FIELD CAPACITY CORRESPONDS TO HYDARULIC CONDUCTIVITY OF
  !> 0.10 MM/DAY -> 1.157x1E-09 M/S
  !>
  do J = 1,IG
    do I = IL1,IL2
      !
      if (ISAND(I,J) == - 3 .or. ISAND(I,J) == - 4) then
        SM_FUNC(I,J) = 0.0
      else ! I.E., ISAND/=-3 OR -4
        if (THLIQ(I,J) <= THLW(I,J)) then
          SM_FUNC(I,J) = 0.0
        else if (THLIQ(I,J) > THLW(I,J) .and. &
       THLIQ(I,J) < THFC(I,J)) then
          SM_FUNC(I,J) = (THLIQ(I,J) - THLW(I,J)) / &
                         (THFC(I,J) - THLW(I,J))
        else if (THLIQ(I,J) >= THFC(I,J)) then
          SM_FUNC(I,J) = 1.0
        end if
        !
      end if ! ISAND==-3 OR -4
    end do ! loop 510
  end do ! loop 500
  !
  K1 = 0
  do J = 1,IC
    if (J == 1) then
      K1 = K1 + 1
    else
      K1 = K1 + NOL2PFTS(J - 1)
    end if
    K2 = K1 + NOL2PFTS(J) - 1
    do M = K1,K2
      do I = IL1,IL2
        do N = 1,IG
          if (ISAND(I,N) /= - 3) then ! ONLY FOR NON-BEDROCK LAYERS
            SM_FUNC2(I,N) = (1.0 - (1.0 - SM_FUNC(I,N)) ** SN(SORT(M)))
            SM_FUNC2(I,N) = SM_FUNC2(I,N) + (1.0 - SM_FUNC2(I,N)) &
                            * SMSCALE(SORT(M))
            AVE_SM_FUNC(I,M) = AVE_SM_FUNC(I,M) + SM_FUNC2(I,N) * RMAT(I,M,N)
            TOT_RMAT(I,M) = TOT_RMAT(I,M) + RMAT(I,M,N)

          end if
        end do ! loop 535

        AVE_SM_FUNC(I,M) = AVE_SM_FUNC(I,M) / TOT_RMAT(I,M)

        if (TOT_RMAT(I,M) < 0.9) then
          write(6, * )'PFT = ',M,' I = ',I
          write(6, * )'RMAT ADD = ',TOT_RMAT(I,M)
          call errorHandler('PHTSYN', - 99)
        end if

      end do ! loop 530
    end do ! loop 525
  end do ! loop 520
  !>
  !> USE SOIL MOISTURE FUNCTION TO MAKE Vm, unstressed -> Vm STRESSED
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (COSZS(I) > 0.0) then
          if (LEAFOPT == 1) then
            VM(I,J) = VMUNS(I,J) * AVE_SM_FUNC(I,J)
          else if (LEAFOPT == 2) then
            VM_SUN(I,J) = VMUNS_SUN(I,J) * AVE_SM_FUNC(I,J)
            VM_SHA(I,J) = VMUNS_SHA(I,J) * AVE_SM_FUNC(I,J)
          end if
        end if
        !
      end if
    end do ! loop 550
  end do ! loop 540
  !>
  !> FIND TEMPERATURE DEPENDENT PARAMETER VALUES
  !>
  do I = IL1,IL2
    !>
    !> FIND RUBISCO SPECIFICITY FOR \f$CO_2\f$ RELATIVE TO \f$O_2\f$ - SIGMA
    !>
    Q10 = 0.57
    Q10_FUNC = Q10 ** (0.1 * (TCAN(I) - 298.16))
    SIGMA(I) = 2600.0 * Q10_FUNC
    !>
    !> FIND \f$CO_2\f$ COMPENSATION POINT USING RUBISCO SPECIFICITY - TGAMMA.
    !> KEEP IN MIND THAT \f$CO_2\f$ COMPENSATION POINT FOR C4 PLANTS IS ZERO,
    !> SO THE FOLLOWING VALUE IS RELEVANT FOR C3 PLANTS ONLY
    !>
    O2_CONC(I) = .2095 * PRESSG(I)
    TGAMMA(I) = O2_CONC(I) / (2.0 * SIGMA(I))
    !>
    !> ESTIMATE MICHELIS-MENTON CONSTANTS FOR \f$CO_2\f$ (Kc) and \f$O_2\f$ (Ko) TO
    !> BE USED LATER FOR ESTIMATING RUBISCO LIMITED PHOTOSYNTHETIC RATE
    !>
    Q10 = 2.10
    Q10_FUNC = Q10 ** (0.1 * (TCAN(I) - 298.16))
    KC(I) = 30.0 * Q10_FUNC
    !
    Q10 = 1.20
    Q10_FUNC = Q10 ** (0.1 * (TCAN(I) - 298.16))
    KO(I) = 30000.0 * Q10_FUNC
    !
  end do ! loop 570
  !>
  !> CHOOSE A VALUE OF INTERCELLULAR \f$CO_2\f$ CONCENTRATION \f$(CO_2i)\f$ IF STARTING
  !> FOR THE FIRST TIME, OR USE VALUE FROM THE PREVIOUS TIME STEP
  !>
  IT_COUNT = 0
  !ignoreLint(1)
999 CONTINUE ! Anchor for a GO TO statement that needs to be removed
  !
  do J = 1,ICC
    do I = IL1,IL2
      !
      if (LEAFOPT == 1) then
        CO2I(I,J) = CO2I1(I,J)
        if (CO2I(I,J) <= ZERO) then
          CO2I(I,J) = INICO2I(SORT(J)) * CO2A(I)
        end if
      else if (LEAFOPT == 2) then
        CO2I_SUN(I,J) = CO2I1(I,J)
        if (CO2I_SUN(I,J) <= ZERO) then
          CO2I_SUN(I,J) = INICO2I(SORT(J)) * CO2A(I)
        end if
        CO2I_SHA(I,J) = CO2I2(I,J)
        if (CO2I_SHA(I,J) <= ZERO) then
          CO2I_SHA(I,J) = INICO2I(SORT(J)) * CO2A(I)
        end if
      end if
      !
    end do ! loop 590
  end do ! loop 580
  !>
  !> ESTIMATE RUBISCO LIMITED PHOTOSYNTHETIC RATE
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (COSZS(I) > 0.0) then
          if (LEAFOPT == 1) then
            JC1(I,J) = CO2I(I,J) - TGAMMA(I)
            JC2(I,J) = KC(I) * (1.0 + (O2_CONC(I) / KO(I)) )
            JC3(I,J) = CO2I(I,J) + JC2(I,J)
            !
            if (ISC4(SORT(J))) then
              JC(I,J)  = VM(I,J)
            else
              JC(I,J)  = VM(I,J) * (JC1(I,J) / JC3(I,J))
            end if
          else if (LEAFOPT == 2) then
            JC1_SUN(I,J) = CO2I_SUN(I,J) - TGAMMA(I)
            JC1_SHA(I,J) = CO2I_SHA(I,J) - TGAMMA(I)
            JC2(I,J) = KC(I) * (1.0 + (O2_CONC(I) / KO(I)) )
            JC3_SUN(I,J) = CO2I_SUN(I,J) + JC2(I,J)
            JC3_SHA(I,J) = CO2I_SHA(I,J) + JC2(I,J)
            !
            if (ISC4(SORT(J))) then
              JC_SUN(I,J) = VM_SUN(I,J)
              JC_SHA(I,J) = VM_SHA(I,J)
            else
              JC_SUN(I,J) = VM_SUN(I,J) * (JC1_SUN(I,J) / JC3_SUN(I,J))
              JC_SHA(I,J) = VM_SHA(I,J) * (JC1_SHA(I,J) / JC3_SHA(I,J))
            end if
          end if
        end if
        !
      end if
    end do ! loop 610
  end do ! loop 600
  !>
  !> ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY AVAILABLE LIGHT
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (COSZS(I) > 0.0) then
          if (LEAFOPT == 1) then
            JE1(I,J) = FPAR(I,J) * alpha_phtsyn(SORT(J)) &
                       * (1.0 - omega_phtsyn(SORT(J)))
            JE2(I,J) = ( CO2I(I,J) - TGAMMA(I) ) / &
                       (CO2I(I,J) + (2.0 * TGAMMA(I)) )
            !
            if (ISC4(SORT(J))) then
              JE(I,J) = IPAR(I) * JE1(I,J)
            else
              JE(I,J) = IPAR(I) * JE1(I,J) * JE2(I,J)
            end if
          else if (LEAFOPT == 2) then
            JE1_SUN(I,J) = FPAR_SUN(I,J) * alpha_phtsyn(SORT(J))
            JE1_SHA(I,J) = FPAR_SHA(I,J) * alpha_phtsyn(SORT(J))
            JE2_SUN(I,J) = ( CO2I_SUN(I,J) - TGAMMA(I) ) / &
                           (CO2I_SUN(I,J) + (2.0 * TGAMMA(I)) )
            JE2_SHA(I,J) = ( CO2I_SHA(I,J) - TGAMMA(I) ) / &
                           (CO2I_SHA(I,J) + (2.0 * TGAMMA(I)) )

            if (ISC4(SORT(J))) then
              JE_SUN(I,J) = (IPAR_SUN(I) + IPAR_SHA(I)) * JE1_SUN(I,J)
              JE_SHA(I,J) = IPAR_SHA(I) * JE1_SHA(I,J)
            else
              JE_SUN(I,J) = (IPAR_SUN(I) + IPAR_SHA(I)) * JE1_SUN(I,J) * &
                            JE2_SUN(I,J)
              JE_SHA(I,J) = IPAR_SHA(I) * JE1_SHA(I,J) * JE2_SHA(I,J)
            end if
          end if
        end if
        !
      end if
    end do ! loop 630
  end do ! loop 620
  !>
  !> ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY TRANSPORT CAPACITY
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (COSZS(I) > 0.0) then
          if (LEAFOPT == 1) then
            if (ISC4(SORT(J))) then
              JS(I,J) = 20000.0 * VM(I,J) * (CO2I(I,J) / PRESSG(I))
            else
              JS(I,J) = 0.5 * VM(I,J)
            end if
          else if (LEAFOPT == 2) then
            if (ISC4(SORT(J))) then
              JS_SUN(I,J) = 20000.0 * VM_SUN(I,J) * (CO2I_SUN(I,J) / PRESSG(I))
              JS_SHA(I,J) = 20000.0 * VM_SHA(I,J) * (CO2I_SHA(I,J) / PRESSG(I))
            else
              JS_SUN(I,J) = 0.5 * VM_SUN(I,J)
              JS_SHA(I,J) = 0.5 * VM_SHA(I,J)
            end if
          end if
        end if
        !
      end if
    end do ! loop 650
  end do ! loop 640
  !>
  !> INCLUDE NUTRIENT LIMITATION EFFECT BY DOWN-REGULATING PHOTOSYNTHESIS
  !> N_EFFECT DECREASES FROM 1.0 AS \f$CO_2\f$ INCREASES ABOVE 288 PPM.
  !>
  do I = IL1,IL2
    DELTA_CO2(I) = CO2CONC(I) - 288.0
    TEMP_R = LOG(1.0 + (DELTA_CO2(I) / 288.0))
    TEMP_B = 1.0 + TEMP_R * (GAMMA_W)
    TEMP_C = 1.0 + TEMP_R * (GAMMA_M)
    N_EFFECT(I) = TEMP_B / TEMP_C
    !> LIMIT N_EFFECT TO MAX OF 1.0 SO THAT NO UP-REGULATION OCCURS
    N_EFFECT(I) = MIN(1.0,N_EFFECT(I))
  end do ! loop 641
  !>
  !> FIND THE SMOOTHED AVERAGE OF THREE PHOTOSYNTHETIC RATES JC, JE,
  !> AND JS USING COLLATZ'S TWO QUADRATIC EQUATIONS, OR FIND THE MIN.
  !> OF THIS TWO RATES OR FIND MIN. OF JC AND JE.
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (COSZS(I) > 0.0) then
          if (LEAFOPT == 1) then
            if (SMOOTH) then
              TEMP_B  = 0.0
              TEMP_C  = 0.0
              TEMP_R  = 0.0
              TEMP_Q1 = 0.0
              TEMP_Q2 = 0.0
              TEMP_JP = 0.0
              !
              TEMP_B  = JC(I,J) + JE(I,J)
              TEMP_C  = JC(I,J) * JE(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA1 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_JP = MIN(TEMP_Q1,TEMP_Q2)
              !
              TEMP_B  = TEMP_JP + JS(I,J)
              TEMP_C  = TEMP_JP * JS(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA2 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
              A_VEG(I,J) = MIN(TEMP_Q1,TEMP_Q2)
            else if (MIN2) then
              A_VEG(I,J) = MIN(JC(I,J),JE(I,J))
            else if (MIN3) then
              A_VEG(I,J) = MIN(JC(I,J),JE(I,J),JS(I,J))
            else
              call errorHandler('PHTSYN', - 1)
            end if
            !> DOWN-REGULATE PHOTOSYNTHESIS FOR C3 PLANTS
            if (.not. ISC4(SORT(J))) then
              A_VEG(I,J) = A_VEG(I,J) * N_EFFECT(I)
            end if
            A_VEG(I,J) = MAX(0.0,A_VEG(I,J))
          else if (LEAFOPT == 2) then
            if (SMOOTH) then
              TEMP_B  = 0.0
              TEMP_C  = 0.0
              TEMP_R  = 0.0
              TEMP_Q1 = 0.0
              TEMP_Q2 = 0.0
              TEMP_JP = 0.0
              !
              TEMP_B  = JC_SUN(I,J) + JE_SUN(I,J)
              TEMP_C  = JC_SUN(I,J) * JE_SUN(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA1 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_JP = MIN(TEMP_Q1,TEMP_Q2)

              TEMP_B  = TEMP_JP + JS_SUN(I,J)
              TEMP_C  = TEMP_JP * JS_SUN(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA2 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
              A_VEG_SUN(I,J) = MIN(TEMP_Q1,TEMP_Q2)
            else if (MIN2) then
              A_VEG_SUN(I,J) = MIN(JC_SUN(I,J),JE_SUN(I,J))
            else if (MIN3) then
              A_VEG_SUN(I,J) = MIN(JC_SUN(I,J),JE_SUN(I,J),JS_SUN(I,J))
            else
              call errorHandler('PHTSYN', - 2)
            end if
            if (.not. ISC4(SORT(J))) then
              A_VEG_SUN(I,J) = A_VEG_SUN(I,J) * N_EFFECT(I)
            end if
            A_VEG_SUN(I,J) = MAX(0.0,A_VEG_SUN(I,J))
            !
            if (SMOOTH) then
              TEMP_B  = 0.0
              TEMP_C  = 0.0
              TEMP_R  = 0.0
              TEMP_Q1 = 0.0
              TEMP_Q2 = 0.0
              TEMP_JP = 0.0
              !
              TEMP_B  = JC_SHA(I,J) + JE_SHA(I,J)
              TEMP_C  = JC_SHA(I,J) * JE_SHA(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA1 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
              TEMP_JP = MIN(TEMP_Q1,TEMP_Q2)
              !
              TEMP_B  = TEMP_JP + JS_SHA(I,J)
              TEMP_C  = TEMP_JP * JS_SHA(I,J)
              TEMP_R  = MAX( (TEMP_B ** 2 - 4. * BETA2 * TEMP_C),0.0)
              TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
              TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
              A_VEG_SHA(I,J) = MIN(TEMP_Q1,TEMP_Q2)
            else if (MIN2) then
              A_VEG_SHA(I,J) = MIN(JC_SHA(I,J),JE_SHA(I,J))
            else if (MIN3) then
              A_VEG_SHA(I,J) = MIN(JC_SHA(I,J),JE_SHA(I,J),JS_SHA(I,J))
            end if
            if (.not. ISC4(SORT(J))) then
              A_VEG_SHA(I,J) = A_VEG_SHA(I,J) * N_EFFECT(I)
            end if
            A_VEG_SHA(I,J) = MAX(0.0,A_VEG_SHA(I,J))
          end if
        end if
        !
      end if
    end do ! loop 670
  end do ! loop 660
  !>
  !> ESTIMATE LEAF MAINTENANCE RESPIRATION RATES AND NET PHOTOSYNTHETIC
  !> RATE. THIS NET PHOSYNTHETIC RATE IS /M^2 OF VEGETATED LAND.
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !>
        !> RECENT STUDIES SHOW RmL IS LESS TEMPERATURE SENSITIVE THAN
        !> PHOTOSYNTHESIS DURING DAY, THAT'S WHY A SMALL Q10 VALUE IS
        !> USED DURING DAY.
        !>
        Q10_FUNCN = 2.00 ** (0.1 * (TCAN(I) - 298.16))
        Q10_FUNCD = 1.30 ** (0.1 * (TCAN(I) - 298.16))
        !
        if (LEAFOPT == 1) then
          if (COSZS(I) > 0.0) then
            RML_VEG(I,J) = RMLCOEFF(SORT(J)) * VMAXC(I,J) * Q10_FUNCD
          else
            RML_VEG(I,J) = RMLCOEFF(SORT(J)) * VMAXC(I,J) * Q10_FUNCN
          end if
          AN_VEG(I,J) = A_VEG(I,J) - RML_VEG(I,J)
        else if (LEAFOPT == 2) then
          if (COSZS(I) > 0.0) then
            RML_SUN(I,J) = RMLCOEFF(SORT(J)) * VMAXC_SUN(I,J) * Q10_FUNCD
            RML_SHA(I,J) = RMLCOEFF(SORT(J)) * VMAXC_SHA(I,J) * Q10_FUNCD
          else
            RML_SUN(I,J) = RMLCOEFF(SORT(J)) * VMAXC_SUN(I,J) * Q10_FUNCN
            RML_SHA(I,J) = RMLCOEFF(SORT(J)) * VMAXC_SHA(I,J) * Q10_FUNCN
          end if
          AN_SUN(I,J) = A_VEG_SUN(I,J) - RML_SUN(I,J)
          AN_SHA(I,J) = A_VEG_SHA(I,J) - RML_SHA(I,J)
        end if
        !
      end if
    end do ! loop 690
  end do ! loop 680
  !>
  !> FIND \f$CO_2\f$ CONCENTRATION AT LEAF SURFACE FOR ALL VEGETATION TYPES.
  !> ALTHOUGH WE ARE FINDING \f$CO_2\f$ CONC AT THE LEAF SURFACE SEPARATELY
  !> FOR ALL VEGETATION TYPES, THE BIG ASSUMPTION HERE IS THAT THE
  !> AERODYNAMIC CONDUCTANCE IS SAME OVER ALL VEGETATION TYPES. CLASS
  !> FINDS AERODYNAMIC RESISTANCE OVER ALL THE 4 SUB-AREAS, BUT NOT
  !> FOR DIFFERENT VEGETATION TYPES WITHIN A SUB-AREA.
  !> ALSO CHANGE AERODYNAMIC CONDUCTANCE, CFLUX, FROM M/S TO \f$MOL/M^2/S\f$
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        GB(I) = CFLUX(I) * (TFREZ / TCAN(I)) * (PRESSG(I) / STD_PRESS) &
                * (1. / 0.0224)
        !
        GB(I) = MIN(10.0,MAX(0.1,GB(I)) )
        !
        if (LEAFOPT == 1) then
          TEMP_AN = AN_VEG(I,J)
          CO2LS(I,J) = 0.5 * (CO2LS(I,J) + &
                       (CO2A(I) - ( (TEMP_AN * 1.37 * PRESSG(I)) / GB(I)))  )
          CO2LS(I,J) = MAX (1.05 * TGAMMA(I),CO2LS(I,J))
        else if (LEAFOPT == 2) then
          TEMP_AN = AN_SUN(I,J)
          CO2LS_SUN(I,J) = 0.5 * (CO2LS_SUN(I,J) + &
                           (CO2A(I) - ( (TEMP_AN * 1.37 * PRESSG(I)) / GB(I))) )
          CO2LS_SUN(I,J) = MAX (1.05 * TGAMMA(I),CO2LS_SUN(I,J))
          !
          TEMP_AN = AN_SHA(I,J)
          CO2LS_SHA(I,J) = 0.5 * (CO2LS_SHA(I,J) + &
                           (CO2A(I) - ( (TEMP_AN * 1.37 * PRESSG(I)) / GB(I))) )
          CO2LS_SHA(I,J) = MAX (1.05 * TGAMMA(I),CO2LS_SHA(I,J))
        end if
        !
      end if
    end do ! loop 710
  end do ! loop 700
  !>
  !> FIND STOMATAL CONDUCTANCE AS PER BALL-WOODROW-BERRY FORMULATION
  !> USED BY COLLATZ ET AL. OR USE THE LEUNING TYPE FORMULATION WHICH
  !> USES VPD INSTEAD OF RH
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !>
        !> IF LIGHT IS TOO LESS MAKE PARAMETER BB VERY SMALL
        if (QSWV(I) < 2.0) then
          USEBB(J) = 0.001
        else
          USEBB(J) = BB(SORT(J))
        end if
        RH(I) = MAX(0.3,RH(I))          ! FOLLOWING IBIS
        !
        if (LEAFOPT == 1) then
          TEMP_AN = AN_VEG(I,J)
          if (PS_COUP == 1) then
            GC(I,J) = ( (MM(SORT(J)) * RH(I) * PRESSG(I) * TEMP_AN) / &
                      CO2LS(I,J) ) + USEBB(J) * USEAILCG(I,J) * AVE_SM_FUNC(I,J)
          else if (PS_COUP == 2) then
            GC(I,J) = ( (MM(SORT(J)) * VPD_TERM(I,J) * PRESSG(I) * TEMP_AN) / &
                      (CO2LS(I,J) - TGAMMA(I)) ) + &
                      USEBB(J) * USEAILCG(I,J) * AVE_SM_FUNC(I,J)
          end if
          !
          GC(I,J) = MAX(GCMIN(I,J), &
                    USEBB(J) * USEAILCG(I,J) * AVE_SM_FUNC(I,J),GC(I,J))
          GC(I,J) = MIN(GCMAX(I,J),GC(I,J))
        else if (LEAFOPT == 2) then
          TEMP_AN = AN_SUN(I,J)
          if (PS_COUP == 1) then
            GC_SUN(I,J) = ( (MM(SORT(J)) * RH(I) * PRESSG(I) * TEMP_AN) / &
                          CO2LS_SUN(I,J) ) + &
                          USEBB(J) * AILCG_SUN(I,J) * AVE_SM_FUNC(I,J)
          else if (PS_COUP == 2) then
            GC_SUN(I,J) = ((MM(SORT(J)) * VPD_TERM(I,J) * PRESSG(I) * TEMP_AN) / &
                          (CO2LS_SUN(I,J) - TGAMMA(I) ) ) + &
                          USEBB(J) * AILCG_SUN(I,J) * AVE_SM_FUNC(I,J)
          end if
          !
          GC_SUN(I,J) = MAX(GCMIN(I,J), &
                        USEBB(J) * AILCG_SUN(I,J) * AVE_SM_FUNC(I,J),GC_SUN(I,J))
          GC_SUN(I,J) = MIN(GCMAX(I,J),GC_SUN(I,J))
          !
          TEMP_AN = AN_SHA(I,J)
          if (PS_COUP == 1) then
            GC_SHA(I,J) = ( (MM(SORT(J)) * RH(I) * PRESSG(I) * TEMP_AN) / &
                          CO2LS_SHA(I,J) ) + &
                          USEBB(J) * AILCG_SHA(I,J) * AVE_SM_FUNC(I,J)
          else if (PS_COUP == 2) then
            GC_SHA(I,J) = ((MM(SORT(J)) * VPD_TERM(I,J) * PRESSG(I) * TEMP_AN) / &
                          (CO2LS_SHA(I,J) - TGAMMA(I)) ) + &
                          USEBB(J) * AILCG_SHA(I,J) * AVE_SM_FUNC(I,J)
          end if
          !
          GC_SHA(I,J) = MAX(GCMIN(I,J), &
                        USEBB(J) * AILCG_SHA(I,J) * AVE_SM_FUNC(I,J),GC_SHA(I,J))
          GC_SHA(I,J) = MIN(GCMAX(I,J),GC_SHA(I,J))
        end if
        !
      end if
    end do ! loop 730
  end do ! loop 720
  !>
  !> FIND THE INTERCELLULAR \f$CO_2\f$ CONCENTRATION BASED ON ESTIMATED VALUE OF GC
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (LEAFOPT == 1) then
          TEMP_AN = AN_VEG(I,J)
          CO2I(I,J) = 0.5 * (CO2I(I,J) + (CO2LS(I,J) - &
                      ( (TEMP_AN * 1.65 * PRESSG(I)) / GC(I,J) ) ) )
          CO2I(I,J) = MAX(1.05 * TGAMMA(I),MIN(CO2IMAX,CO2I(I,J)))
          PREV_CO2I(I,J) = CO2I(I,J)
        else if (LEAFOPT == 2) then
          TEMP_AN = AN_SUN(I,J)
          CO2I_SUN(I,J) = 0.5 * (CO2I_SUN(I,J) + (CO2LS_SUN(I,J) - &
                          ( (TEMP_AN * 1.65 * PRESSG(I)) / GC_SUN(I,J))) )
          CO2I_SUN(I,J) = MAX(1.05 * TGAMMA(I),MIN(CO2IMAX, &
                          CO2I_SUN(I,J)))
          PREV_CO2I_SUN(I,J) = CO2I_SUN(I,J)
          !
          TEMP_AN = AN_SHA(I,J)
          CO2I_SHA(I,J) = 0.5 * (CO2I_SHA(I,J) + (CO2LS_SHA(I,J) - &
                          ( (TEMP_AN * 1.65 * PRESSG(I)) / GC_SHA(I,J))) )
          CO2I_SHA(I,J) = MAX(1.05 * TGAMMA(I),MIN(CO2IMAX, &
                          CO2I_SHA(I,J)))
          PREV_CO2I_SHA(I,J) = CO2I_SHA(I,J)
        end if
        !
      end if
    end do ! loop 750
  end do ! loop 740
  !
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then

        if (LEAFOPT == 1) then
          CO2I1(I,J) = PREV_CO2I(I,J)
          CO2I2(I,J) = 0.0
        else if (LEAFOPT == 2) then
          CO2I1(I,J) = PREV_CO2I_SUN(I,J)
          CO2I2(I,J) = PREV_CO2I_SHA(I,J)
        end if

      end if
    end do ! loop 770
  end do ! loop 760
  !
  IT_COUNT = IT_COUNT + 1
  !>
  !> SEE IF WE HAVE PERFORMED THE REQUIRED NO. OF ITERATIONS, IF NOT
  !> THEN WE GO BACK AND DO ANOTHER ITERATION
  !>
  if (IT_COUNT < REQITER) then
    GO TO 999
  end if
  !>
  !> WHEN REQUIRED NO. OF ITERATIONS HAVE BEEN PERFORMED THEN FIND
  !> STOMATAL CONDUCTANCES FOR ALL VEGETATION TYPES IN M/S AND THEN
  !> USE CONDUCTANCES TO FIND RESISTANCES. GCTU IMPLIES GC IN TRADITIONAL
  !> UNITS OF M/S
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (FCANC(I,J) > ZERO) then
        !
        if (LEAFOPT == 1) then
          GCTU(I,J) = GC(I,J) * (TCAN(I) / TFREZ) * &
                      (STD_PRESS / PRESSG(I)) * 0.0224
          RC_VEG(I,J) = 1. / GCTU(I,J)
        else if (LEAFOPT == 2) then
          GCTU_SUN(I,J) = GC_SUN(I,J) * (TCAN(I) / TFREZ) * &
                          (STD_PRESS / PRESSG(I)) * 0.0224
          GCTU_SHA(I,J) = GC_SHA(I,J) * (TCAN(I) / TFREZ) * &
                          (STD_PRESS / PRESSG(I)) * 0.0224
          !
          if (COSZS(I) < 0.0 .or. QSWV(I) < 2.0) then
            !> DON'T WANT TO REDUCE RESISTANCE AT NIGHT TO LESS THAN
            !> OUR MAX. VALUE OF AROUND 5000 S/M
            GCTU(I,J) = 0.5 * (GCTU_SUN(I,J) + GCTU_SHA(I,J))
          else
            GCTU(I,J) = GCTU_SUN(I,J) + GCTU_SHA(I,J)
          end if
          !
          RC_VEG(I,J) = 1. / GCTU(I,J)
          AN_VEG(I,J) = AN_SUN(I,J) + AN_SHA(I,J)
          RML_VEG(I,J) = RML_SUN(I,J) + RML_SHA(I,J)
        end if
        !
      end if
    end do ! loop 790
  end do ! loop 780
  !>
  !> IF USING STORAGE LAI THEN WE SET STOMATAL RESISTANCE TO ITS MAXIMUM VALUE.
  !>
  do J = 1,ICC
    do I = IL1,IL2
      if (USESLAI(I,J) == 1 .and. AILCG(I,J) < 0.2) then
        RC_VEG(I,J) = 5000.0
      end if
    end do ! loop 810
  end do ! loop 800
  !>
  !> AND FINALLY TAKE WEIGHTED AVERAGE OF RC_VEG BASED ON FRACTIONAL
  !> COVERAGE OF OUR 4 VEGETATION TYPES
  !>
  do J = 1,ICC
    do I = IL1,IL2
      RC(I) = RC(I) + FCANC(I,J) * RC_VEG(I,J)
    end do ! loop 830
  end do ! loop 820
  !
  do I = IL1,IL2
    FC_TEST(I) = SUM(FCANC(I,:))
    if (FC_TEST(I) > ZERO) then
      RC(I) = RC(I) / FC_TEST(I)
    else
      RC(I) = 5000.0
    end if
  end do ! loop 840
  !>
  !> CONVERT AN_VEG AND RML_VEG TO u-MOL CO2/M2.SEC
  !>
  do J = 1,ICC
    do I = IL1,IL2
      AN_VEG(I,J) = AN_VEG(I,J) * 1.0E+06
      RML_VEG(I,J) = RML_VEG(I,J) * 1.0E+06
    end do ! loop 880
  end do ! loop 870
  !

  DEALLOCATE(USESLAI)
  DEALLOCATE(SORT)
  DEALLOCATE(USEBB)

  DEALLOCATE(FC_TEST)
  DEALLOCATE(USEAILCG)
  DEALLOCATE(SM_FUNC)
  DEALLOCATE(SM_FUNC2)
  DEALLOCATE(AVE_SM_FUNC)
  DEALLOCATE(TOT_RMAT)

  DEALLOCATE(VMAXC)
  DEALLOCATE(JE3)
  DEALLOCATE(VMUNS1,VMUNS2,VMUNS3)
  DEALLOCATE(VMUNS)
  DEALLOCATE(VM)
  DEALLOCATE(SIGMA)
  DEALLOCATE(TGAMMA)
  DEALLOCATE(KC)
  DEALLOCATE(KO)
  DEALLOCATE(CO2I,PREV_CO2I)
  DEALLOCATE(FPAR)
  DEALLOCATE(JC)
  DEALLOCATE(JC1,JC2,JC3)
  DEALLOCATE(JE,JE1,JE2)
  DEALLOCATE(IPAR)
  DEALLOCATE(JS)
  DEALLOCATE(A_VEG)
  DEALLOCATE(GB)
  DEALLOCATE(RC_VEG)
  DEALLOCATE(GCTU)
  DEALLOCATE(GCMIN)
  DEALLOCATE(GCMAX)
  DEALLOCATE(RH)

  DEALLOCATE(VPD,VPD_TERM)

  DEALLOCATE(CO2LS)
  DEALLOCATE(GC)
  DEALLOCATE(O2_CONC)
  DEALLOCATE(CO2A)

  DEALLOCATE(GDIR,KB)
  DEALLOCATE(FPAR_SUN,FPAR_SHA)
  DEALLOCATE(VMAXC_SUN,VMAXC_SHA)
  DEALLOCATE(VMUNS1_SUN,VMUNS1_SHA)
  DEALLOCATE(VMUNS_SUN,VMUNS_SHA)
  DEALLOCATE(VM_SUN,VM_SHA)
  DEALLOCATE(CO2I_SUN,PREV_CO2I_SUN)
  DEALLOCATE(CO2I_SHA,PREV_CO2I_SHA)
  DEALLOCATE(JC1_SUN,JC1_SHA)
  DEALLOCATE(JC3_SUN,JC3_SHA)
  DEALLOCATE(JC_SUN,JC_SHA)
  DEALLOCATE(IPAR_SUN,IPAR_SHA)
  DEALLOCATE(JE1_SUN,JE1_SHA)
  DEALLOCATE(JE2_SUN,JE2_SHA)
  DEALLOCATE(JE_SUN,JE_SHA)
  DEALLOCATE(JS_SUN,JS_SHA)
  DEALLOCATE(A_VEG_SUN,A_VEG_SHA)
  DEALLOCATE(RML_SUN,RML_SHA)
  DEALLOCATE(AN_SUN,AN_SHA)
  DEALLOCATE(CO2LS_SUN,CO2LS_SHA)
  DEALLOCATE(AILCG_SUN,AILCG_SHA)
  DEALLOCATE(GC_SUN,GC_SHA)
  DEALLOCATE(GCTU_SUN,GCTU_SHA)


  return
end subroutine photosynCanopyConduct
!> \file
!> All biogeochemical processes in CLASSIC are simulated at a daily time step except gross photosynthetic uptake and
!! associated calculation of canopy conductance, which are simulated on a half hour time step with CLASS (physics).
!! The photosynthesis module of CLASSIC calculates the net canopy photosynthesis rate, which, together with
!! atmospheric \f$CO_2\f$ concentration and vapour pressure or relative humidity, is used to calculate canopy
!! conductance. This canopy conductance is then used by CLASSIC in its energy and water balance calculations.
!!
!! The photosynthesis parametrization is based upon the approach of Farquhar et al. (1980) \cite Farquhar1980-96e and
!! Collatz et al. (1991, 1992) \cite Collatz1991-5bc \cite Collatz1992-jf as implemented in SiB2 (Sellers et al. 1996)
!! \cite Sellers1996-bh and MOSES (Cox et al. 1999) \cite Cox1999-ia with some minor modifications as described in
!! Arora (2003)\cite Arora2003-3b7. Arora (2003) \cite Arora2003-3b7 outlines four possible configurations for the
!! model based on choice of a \f$\textit{big-leaf}\f$ or \f$\textit{two-leaf}\f$ (sunlight and shaded leaves) mode
!! and stomatal conductance formulations based on either Ball et a. (1987) \cite Ball1987-ou or Leuning (1995)
!! \cite Leuning1995-ab. The Ball et al. (1987) \cite Ball1987-ou formulation uses relative humidity while
!! Leuning (1995) \cite Leuning1995-ab uses vapour pressure deficit in calculation of canopy conductance.
!! While the model remains capable of all four possible configurations, in practice, the model is usually run using the
!! big-leaf parametrization with the stomatal conductance formulation of \cite Leuning1995-ab, which is the
!! configuration described here. The original description of the CLASSIC photosynthesis
!! parametrization in Arora (2003) \cite Arora2003-3b7 did not include discussion of all the PFTs
!! simulated by CLASSIC, which we expand upon here.
!!
!! The gross leaf photosynthesis rate, \f$G_\mathrm{o}\f$, depends upon the maximum assimilation rate allowed by the light
!! (\f$J_\mathrm{e}\f$), Rubisco (\f$J_\mathrm{c}\f$) and transport capacity (\f$J_\mathrm{s}\f$). The limitation placed
!! on \f$G_\mathrm{o}\f$ by the amount of available light is calculated as (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$)
!!
!! \f[
!!  J_\mathrm{e} = \left\{\begin{array}{l l}\varepsilon\, (1-{\nu})I \left[\frac{c_{i} - \Gamma}{c_{i} + 2\Gamma}\right], \qquad C_3 plants\\
!! J_\mathrm{e} = \varepsilon\, (1-{\nu})I, \qquad C_4 plants \end{array} \right. \qquad (Eqn 1)
!! \f]
!! where \f$I\f$ is the incident photosynthetically active radiation (\f$PAR\f$;\f$mol\, photons\, m^{-2}\, s^{-1}\f$), \f${\nu}\f$
!! is the leaf scattering coefficient, with values of 0.15 and 0.17 for \f$C_3\f$ and \f$C_4\f$ plants, respectively,
!! and \f$\varepsilon\f$ is the quantum efficiency (\f$mol\, {CO_2}\, (mol\, photons)^{-1}\f$; values of 0.08 and 0.04 are
!! used for \f$C_3\f$ and \f$C_4\f$ plants, respectively). \f$c_\mathrm{i}\f$ is the partial pressure of \f$CO_2\f$
!! in the leaf interior (\f$Pa\f$) and \f$\Gamma\f$ is the \f$CO_2\f$ compensation point (\f$Pa\f$) (described below).
!!
!! The Rubisco enzyme limited photosynthesis rate, \f$J_\mathrm{c}\f$, is given by
!! \f[
!! J_\mathrm{c} = \left\{\begin{array}{l l} V_\mathrm{m} \left[\frac{c_\mathrm{i} - \Gamma}{c_\mathrm{i} + K_\mathrm{c}(1
!! + O_\mathrm{a}/K_\mathrm{o})}\right], \qquad C_3 plants\\
!! J_\mathrm{c} = V_\mathrm{m}, \qquad C_4 plants \end{array} \right. \qquad(Eqn 2)
!! \f]
!!
!! where \f$V_\mathrm{m}\f$ is the maximum catalytic capacity of Rubisco (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), adjusted
!! for temperature and soil moisture, as described below. \f$K_\mathrm{o}\f$ and \f$K_\mathrm{c}\f$ are the
!! Michaelis--Menten constants for \f$O_2\f$ and \f$CO_2\f$, respectively. \f$O_\mathrm{a}\f$ is the partial pressure (\f$Pa\f$) of oxygen.
!!
!! The transport capacity (\f$J_\mathrm{s}\f$) limitation determines the maximum capacity to transport the
!! products of photosynthesis for \f$C_3\f$ plants, while for \f$C_4\f$ plants it represents \f$CO_2\f$ limitation
!! \f[
!! J_\mathrm{s} = \left\{\begin{array}{l l} 0.5 V_\mathrm{m}, \qquad C_3 plants\\ 2
!! \times 10^4\, V_\mathrm{m} \frac{c_\mathrm{i}}{p}, \qquad C_4 plants \end{array} \right.
!! \qquad (Eqn 3)\f]
!!
!! where \f$p\f$ is surface atmospheric pressure (\f$Pa\f$).
!!
!! \f$V_\mathrm{m}\f$ is calculated as
!! \f[
!! V_\mathrm{m} = \frac{V_{max}f_{25}(2.0)S_{root}(\theta) \times 10^{-6}}
!! {[1+ \exp{0.3(T_\mathrm{c} - T_{high})}][1 + \exp{0.3(T_{low} - T_\mathrm{c})}]}, \label{V_m}
!! \qquad (Eqn 4)\f]
!!
!! where \f$T_\mathrm{c}\f$ is the canopy temperature (\f$C\f$) and \f$T_{low}\f$ and \f$T_{high}\f$ are PFT-dependent
!! lower and upper temperature limits for photosynthesis (see also classicParams.f90). \f$f_{25}\f$ is the standard
!! \f$Q_{10}\f$ function at \f$25\, C\f$ (\f$(f_{25}(Q_{10}) = Q^{(0.1(T_\mathrm{c}-25))}_{10}\f$) and \f$V_{max}\f$
!! is the PFT-dependent maximum rate of carboxylation by the enzyme Rubisco (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$; see
!! also classicParams.f90). The constant \f$10^{-6}\f$ converts \f$V_{max}\f$ from units of
!! \f${\mu}mol\, CO_2\, m^{-2}\, s^{-1}\f$ to \f$mol\, CO_2\, m^{-2}\, s^{-1}\f$.
!!
!! The influence of soil moisture stress is simulated via \f$S_{root}(\theta)\f$, which represents a soil moisture stress term formulated as
!! \f[
!! S_{root}(\theta) = \sum_{i=1}^g S(\theta_i) r_{i}, \qquad (Eqn 5)\f]
!!
!! \f[
!! \label{soilmoist_str} S(\theta_i) = \left[1 - \left\{1 - \phi_i \right\}\right]^\varrho,
!! \qquad (Eqn 6)\f]
!!
!! where \f$S_{root}(\theta)\f$ is calculated by weighting \f$S(\theta_i)\f$ with the fraction of roots, \f$r_{i}\f$, in
!! each soil layer \f$i\f$ and \f$\varrho\f$ is a PFT-specific sensitivity to soil moisture stress (unitless; see also
!! classicParams.f90).  \f$\phi_i\f$ is the degree of soil saturation (soil wetness) given by
!! \f[
!! \label{phitheta} \phi_{i}(\theta_{i}) = \max \left[0, \min \left(1, \frac{\theta_{i} - \theta_{i, wilt}}{\theta_{i, field} - \theta_{i, wilt}} \right) \right],
!! \qquad (Eqn 7)\f]
!!
!! where \f$\theta_{i}\f$ is the volumetric soil moisture (\f$m^{3} water\, (m^{3} soil)^{-1}\f$) of the \f$i\f$th soil
!! layer and \f$\theta_{i, field}\f$ and \f$\theta_{i, wilt}\f$ the soil moisture at field capacity and wilting point, respectively.
!!
!! The \f$CO_2\f$ compensation point (\f$\Gamma\f$) is the \f$CO_2\f$ partial pressure where photosynthetic uptake
!! equals the leaf respiratory losses (used in Eqs. 1 and 2). \f$\Gamma\f$ is zero for \f$C_4\f$ plants
!! but is sensitive to oxygen partial pressure for \f$C_3\f$ plants as
!! \f[
!! \label{co2comp} \Gamma = \left\{\begin{array}{l l} \frac{O_\mathrm{a}}{2 \sigma}, C_3 plants\\ 0, C_4 plants, \end{array} \right.\qquad (Eqn 8)
!! \f]
!!
!! where \f$\sigma\f$ is the selectivity of Rubisco for \f$CO_2\f$ over \f$O_2\f$ (unitless), estimated by
!! \f$\sigma = 2600f_{25}(0.57)\f$. The \f$CO_2\f$ (\f$K_\mathrm{c}\f$) and \f$O_2\f$ (\f$K_\mathrm{o}\f$)
!! Michaelis--Menten constants used in Eq. 2 are determined via
!! \f[
!! \label{K_c} K_\mathrm{c} = 30f_{25}(2.1), \qquad (Eqn 9)
!! \f]
!! \f[
!! \label{K_o} K_o = 3 \times 10^4 f_{25}(1.2).\qquad (Eqn 9)
!! \f]
!!
!! Given the light (\f$J_\mathrm{e}\f$), Rubsico (\f$J_\mathrm{c}\f$) and transportation capacity
!! (\f$J_\mathrm{s}\f$) limiting rates, the leaf-level gross photosynthesis rate,
!! \f$G_\mathrm{o}\f$ (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), is then determined following a
!! minimization based upon smallest roots of the following two quadratic equations
!! \f[
!! J_\mathrm{p} = \frac{(J_\mathrm{c} + J_\mathrm{e}) \pm \sqrt{(J_\mathrm{c} + J_\mathrm{e})^2 - 4\beta_1
!! (J_\mathrm{c} + J_\mathrm{e})}}{2\beta_1} , \qquad (Eqn 10)
!! \f]
!! \f[
!! \label{G_o}G_\mathrm{o} = \frac{(J_\mathrm{p} + J_\mathrm{s}) \pm \sqrt{(J_\mathrm{p} + J_\mathrm{s})^2
!! - 4\beta_2 (J_\mathrm{p} + J_\mathrm{s})}}{2\beta_2}, \qquad (Eqn 10)
!! \f]
!!
!! where \f$\beta_1\f$ is 0.95 and \f$\beta_2\f$ is 0.99. When soil moisture stress is occurring, both
!! the \f$J_\mathrm{s}\f$ and \f$J_\mathrm{c}\f$ terms are reduced since the \f$V_\mathrm{m}\f$ term
!! (Eq. 4) includes the effect of soil moisture stress through the \f$S(\theta)\f$ term and this reduces the
!! leaf-level gross photosynthesis rate.
!!
!! The current version of CLASSIC does not include nutrient constraints on photosynthesis and, as a result,
!! increasing atmospheric \f$CO_2\f$ concentration leads to unconstrained increase in photosynthesis.
!! In natural ecosystems, however, down regulation of photosynthesis occurs due to constraints imposed
!! by availability of nitrogen, as well as phosphorus. To capture this effect,
!! CTEM uses a nutrient limitation term, based on experimental plant growth studies, to down regulate the photosynthetic
!! response to elevated \f$CO_2\f$ concentrations (Arora et al., 2009)\cite Arora2009-9bc. The parametrization, and its
!! rationale, are fully described in Arora et al. (2009) \cite Arora2009-9bc but the basic relations are summarized here.
!! The leaf-level gross photosynthetic rate is scaled by the down-regulation term, \f$\Xi_\mathrm{N}\f$,
!! to yield the nutrient limited leaf level gross photosynthetic rate as
!! \f[
!! G_{\mathrm{o}, N-limited} = \Xi_\mathrm{N} G_\mathrm{o}, \\ \Xi_\mathrm{N}
!! = \frac{1 + \gamma_{gd} \ln(c_\mathrm{a}/c_{0})}{1 + \gamma_g \ln(c_\mathrm{a}/c_{0})},
!! \qquad (Eqn 11)\f]
!!
!! where \f$c_\mathrm{a}\f$ is the atmospheric \f$CO_2\f$ concentration in ppm, \f$c_{0}\f$ is the pre-industrial
!! \f$CO_2\f$ concentration (\f$285.0\, ppm\f$), \f$\gamma_g\f$ is 0.95 (Arora et al. 2009) \cite Arora2009-9bc.
!! A value of \f$\gamma_{gd}\f$ lower than \f$\gamma_g\f$ ensures that \f$\Xi_\mathrm{N}\f$ gradually decreases
!! from its pre-industrial value of one as \f$c_\mathrm{a}\f$ increases to constrain the rate of increase
!! of photosynthesis with rising atmospheric \f$CO_2\f$ concentrations.
!!
!! Finally, the leaf-level gross photosynthesis rate, \f$G_{\mathrm{o}, N-limited}\f$ is scaled up to the
!! canopy-level, \f$G_{canopy}\f$, by considering the exponential vertical profile of radiation along the depth of the canopy as
!! \f[
!! G_{canopy} = G_{\mathrm{o}, N-limited} f_{PAR}, \\ \label{fpar} f_{PAR} = \frac{1}{k_\mathrm{n}}(1-\exp^{-k_\mathrm{n}LAI}), \qquad (Eqn 12)
!! \f]
!! which yields the gross primary productivity (\f$G_{canopy}\f$, GPP). \f$k_\mathrm{n}\f$ is the extinction
!! coefficient that describes the nitrogen and time-mean photosynthetically absorbed radiation (\f$PAR\f$)
!! profile along the depth of the canopy (see also classicParams.f90) (Ingestad and Lund, 1986; \cite Ingestad1986-td
!! Field and Mooney, 1986 \cite Field1986-kd), and \f$LAI\f$ (\f$m^{2}\, leaf\, (m^{2}\, ground)^{-1}\f$) is the leaf area index.
!!
!! The net canopy photosynthetic rate, \f$G_{canopy, net}\f$ (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), is calculated
!! by subtracting canopy leaf maintenance respiration costs (\f$R_{mL}\f$; see mainres.f) as
!! \f[
!! \label{Gnet} G_{canopy, net} = G_{canopy} - R_{mL}.\qquad (Eqn 13)
!! \f]
!!
!! **Coupling of photosynthesis and canopy conductance**
!!
!!
!! When using the Leuning (1995) \cite Leuning1995-ab approach for photosynthesis--canopy conductance coupling,
!! canopy conductance (\f$g_\mathrm{c}\f$; \f$mol\, m^{-2}\, s^{-1}\f$) is expressed as a function of the net canopy photosynthesis rate, \f$G_{canopy, net}\f$, as
!! \f[
!! \label{canopy_cond} g_\mathrm{c} = m \frac{G_{canopy, net} p}{(c_\mathrm{s} - \Gamma)}\frac{1}{(1+V/V_\mathrm{o})} + b {LAI}\qquad (Eqn 14)
!! \f]
!! where \f$p\f$ is the surface atmospheric pressure (\f$Pa\f$), the parameter \f$m\f$ is set to 9.0 for needle-leaved trees,
!! 12.0 for other \f$C_3\f$ plants and 6.0 for \f$C_4\f$ plants, parameter \f$b\f$ is assigned the values of 0.01 and 0.04
!! for \f$C_3\f$ and \f$C_4\f$ plants, respectively. \f$V\f$ is the vapour pressure deficit (\f$Pa\f$) and the parameter
!! \f$V_\mathrm{o}\f$ is set to \f$2000\, Pa\f$ for trees and \f$1500\, Pa\f$ for crops and grasses. The partial
!! pressure of \f$CO_2\f$ at the leaf surface, \f$c_\mathrm{s}\f$, is found via
!! \f[
!! \label{c_s} c_\mathrm{s} = c_{ap} - \frac{1.37 G_{canopy, net} p}{g_b}.\qquad (Eqn 15)
!! \f]
!!
!! Here, \f$c_{ap}\f$ is the atmospheric \f$CO_2\f$ partial pressure (\f$Pa\f$) and \f$g_b\f$ is the aerodynamic conductance
!! estimated by CLASS (\f$mol\, m^{-2}\, s^{-1}\f$). The intra-cellular \f$CO_2\f$ concentration required in Eqs. 1--3 is calculated as
!! \f[
!! \label{c_i} c_\mathrm{i} = c_\mathrm{s} - \frac{1.65 G_{canopy, net} p}{g_\mathrm{c}}.\qquad (Eqn 16)
!! \f]
!!
!! Since calculations of \f$G_{canopy, net}\f$ and \f$c_\mathrm{i}\f$ depend on each other, the photosynthesis-canopy
!! conductance equations need to be solved iteratively. The initial value of \f$c_\mathrm{i}\f$ used in calculation
!! of \f$G_{canopy, net}\f$ is the value from the previous time step or, in its absence, \f$c_\mathrm{i}\f$ is assumed to be \f$0.7c_{ap}\f$.
!!
!! Canopy (\f$g_\mathrm{c}\f$) and aerodynamic (\f$g_b\f$) conductance used in above calculations are expressed in
!! units of \f$mol\, CO_2\, m^{-2}\, s^{-1}\f$ but can be converted to the traditional units of \f$m\, s^{-1}\f$ as follows
!! \f[
!! g_\mathrm{c} (m\, s^{-1}) = 0.0224\, \frac{T_\mathrm{c}}{T_\mathrm{f}}\, \frac{p_0}{p}\, g_\mathrm{c} (mol\, m^{-2}\, s^{-1}), \qquad (Eqn 17)
!! \f]
!!
!! where \f$p_0\f$ is the standard atmospheric pressure (\f$101\, 325\, Pa\f$) and \f$T_\mathrm{f}\f$ is freezing temperature (\f$273.16\, K\f$).
