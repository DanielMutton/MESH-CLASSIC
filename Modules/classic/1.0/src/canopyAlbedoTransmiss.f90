!> \file
!! Calculates vegetation albedos, transmissivities and
!! stomatal resistances.
!! @author D. Verseghy, M. Lazare, P. Bartlett, R. Harvey, J. Melton
!
subroutine canopyAlbedoTransmiss (ALVSCN, ALIRCN, ALVSCS, ALIRCS, TRVSCN, TRIRCN, & ! Formerly CANALB
                                  TRVSCS, TRIRCS, RC, RCS, &
                                  ALVSC, ALIRC, RSMIN, QA50, VPDA, VPDB, PSIGA, PSIGB, &
                                  FC, FCS, FSNOW, FSNOWC, FSNOCS, FCAN, FCANS, PAI, PAIS, &
                                  AIL, PSIGND, FCLOUD, COSZS, QSWINV, VPD, TA, &
                                  ACVDAT, ACIDAT, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                                  ILG, IL1, IL2, JL, IC, ICP1, IG, IALC, &
                                  CXTEFF, TRVS, TRIR, RCACC, RCG, RCV, RCT, GC)

  !     * Mar 28/17 - S.SUN Expand PFTs from 4 to 5 in CLASS
  !     *
  !     * AUG 04/15 - D.VERSEGHY/M.LAZARE. REMOVE FLAG VALUE OF RC FOR
  !     *                         VERY DRY SOILS.
  !     * SEP 05/14 - P.BARTLETT. INCREASED ALBEDO VALUES FOR SNOW-
  !     *                         COVERED CANOPY.
  !     * JUN 27/13 - D.VERSEGHY/ USE LOWER BOUND OF 0.01 INSTEAD OF
  !     *             M.LAZARE.   0. IN LOOP 900 TO AVOID CRASH
  !     *                         IN EXTREME CASE OF LOW VISIBLE
  !     *                         INCIDENT SUN.
  !     * AUG 17/12 - J.MELTON    ADD CONSTRAINT TO TRCLRV AT LOW COSZS
  !                               OTHERWISE IT HAS NUMERICAL PROBLEMS
  !                               WHEN MODEL IS COMPILED WITH GFORTRAN.
  !     * DEC 21/11 - M.LAZARE.   DEFINE CONSTANTS "EXPMAX1", EXPMAX2",
  !     *                         "EXPMAX3" TO AVOID REDUNDANT EXP
  !     *                         CALCULATIONS.
  !     * OCT 16/08 - R.HARVEY.   ADD LARGE LIMIT FOR EFFECTIVE
  !     *                         EXTINCTION COEFFICIENT (CXTEFF) IN
  !     *                         (RARE) CASES WHEN CANOPY TRANSMISSIVITY
  !     *                         IN THE VISIBLE IS ZERO EXACTLY.
  !     * MAR 25/08 - D.VERSEGHY. DISTINGUISH BETWEEN LEAF AREA INDEX
  !     *                         AND PLANT AREA INDEX.
  !     * OCT 19/07 - D.VERSEGHY. SIMPLIFY ALBEDO CALCULATIONS FOR
  !     *                         SNOW-FREE CROPS AND GRASS; CORRECT
  !     *                         BUG IN CALCULATION OF RC.
  !     * APR 13/06 - P.BARTLETT/D.VERSEGHY. CORRECT OVERALL CANOPY
  !     *                                    ALBEDO, INTRODUCE SEPARATE
  !     *                                    GROUND AND SNOW ALBEDOS FOR
  !     *                                    OPEN OR CANOPY-COVERED AREAS.
  !     * MAR 21/06 - P.BARTLETT. PROTECT RC CALCULATION AGAINST DIVISION
  !     *                         BY ZERO.
  !     * SEP 26/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 600 LOOP.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 05/04 - D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS AGAINST
  !     *                         ROUNDOFF ERRORS.
  !     * JAN 24/02 - P.BARTLETT/D.VERSEGHY. REFINE CALCULATION OF NEW
  !     *                                    STOMATAL RESISTANCES.
  !     * JUL 30/02 - P.BARTLETT/D.VERSEGHY. NEW STOMATAL RESISTANCE
  !     *                                    FORMULATION INCORPORATED.
  !     * MAR 18/02 - D.VERSEGHY. ALLOW FOR ASSIGNMENT OF SPECIFIED TIME-
  !     *                         VARYING VALUES OF VEGETATION SNOW-FREE
  !     *                         ALBEDO.
  !     * NOV 29/94 - M.LAZARE. CLASS - VERSION 2.3.
  !     *                       CALL ABORT CHANGED TO CALL errorHandler TO ENABLE
  !     *                       RUNNING ON PC'S.
  !     * MAY 06/93 - D.VERSEGHY. EXTENSIVE MODIFICATIONS TO CANOPY
  !     *                         ALBEDO LOOPS.
  !     * MAR 03/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CANOPY ALBEDOS AND TRANSMISSIVITIES.
  !
  use classicParams,     only : CANEXT, DELT, ALVSWC, ALIRWC, CXTLRG, &
                                classpfts

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IC, ICP1, IG, IALC
  integer :: I, J, IPTBAD, JPTBAD, JPTBDI
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(inout) :: ALVSCN(ILG)  !< Visible albedo of vegetation over bare ground [ ]
  real, intent(inout) :: ALIRCN(ILG)  !< Near-IR albedo of vegetation over bare ground [ ]
  real, intent(inout) :: ALVSCS(ILG)  !< Visible albedo of vegetation over snow [ ]
  real, intent(inout) :: ALIRCS(ILG)  !< Near-IR albedo of vegetation over snow [ ]
  real, intent(inout) :: TRVSCN(ILG)  !< Visible transmissivity of vegetation over bare
  !! ground \f$[ ] (\tau_c)\f$
  real, intent(inout) :: TRIRCN(ILG)  !< Near-IR transmissivity of vegetation over bare
  !! ground \f$[ ] (\tau_c)\f$
  real, intent(inout) :: TRVSCS(ILG)  !< Visible transmissivity of vegetation over snow
  !! \f$[ ] (\tau_c)\f$
  real, intent(inout) :: TRIRCS(ILG)  !< Near-IR transmissivity of vegetation over snow
  !! \f$[ ] (\tau_c)\f$
  real, intent(inout) :: RC(ILG)  !< Stomatal resistance of vegetation over bare
  !! ground \f$[s m^{-1}] (r_c)\f$
  real, intent(out) :: RCS(ILG)  !< Stomatal resistance of vegetation over snow \f$[s m^{-1}]\f$
  !
  !     * 2-D INPUT ARRAYS.
  !
  real, intent(in) :: ALVSC(ILG,ICP1) !< Background average visible albedo of
  !! vegetation category [ ]
  real, intent(in) :: ALIRC(ILG,ICP1) !< Background average near-infrared albedo of
  !! vegetation category [ ]
  real, intent(in) :: RSMIN(ILG,IC)   !< Minimum stomatal resistance of vegetation
  !! category \f$[s m^{-1}] (r_{s,min})\f$
  real, intent(in) :: QA50(ILG,IC)   !< Reference value of incoming shortwave
  !! radiation for vegetation category (used in
  !! stomatal resistance calculation) \f$[W m^{-2}] (K \downarrow_{1/2})\f$
  real, intent(in) :: VPDA(ILG,IC)   !< Vapour pressure deficit coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) \f$[ ] (c_{v1})\f$
  real, intent(in) :: VPDB(ILG,IC)   !< Vapour pressure deficit coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) \f$[ ] (c_{v2})\f$
  real, intent(in) :: PSIGA(ILG,IC)   !< Soil moisture suction coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) \f$[ ] (c_{\Psi 1})\f$
  real, intent(in) :: PSIGB(ILG,IC)   !< Soil moisture suction coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) \f$[ ] (c_{\Psi 2})\f$
  real, intent(in) :: FCAN(ILG,IC)   !< Fractional coverage of vegetation category
  !! over bare ground \f$[ ] (X_i)\f$
  real, intent(in) :: FCANS(ILG,IC)   !< Fractional coverage of vegetation category
  !! over snow [ ]
  real, intent(in) :: PAI(ILG,IC)   !< Plant area index of vegetation category
  !! over bare ground \f$[ ] (\Lambda_p)\f$
  real, intent(in) :: PAIS(ILG,IC)   !< Plant area index of vegetation category
  !! over snow \f$[ ] (\Lambda_p)\f$
  real, intent(in) :: ACVDAT(ILG,IC)   !< Optional user-specified value of canopy
  !! visible albedo to override CLASS-calculated
  !! value [ ]
  real, intent(in) :: ACIDAT(ILG,IC)   !< Optional user-specified value of canopy
  !! near-infrared albedo to override CLASS-
  !! calculated value [ ]
  real, intent(in) :: AIL(ILG,IC)   !< Leaf area index of vegetation category over
  !! bare ground \f$[ ] (\Lambda)\f$
  !
  !     * 1-D INPUT ARRAYS.
  !
  real, intent(in) :: FC(ILG)  !< Fractional coverage of canopy over bare
  !! ground [ ]
  real, intent(in) :: FCS(ILG)  !< Fractional coverage of canopy over snow [ ]
  real, intent(in) :: FSNOW(ILG)  !< Diagnosed fractional snow coverage [ ]
  real, intent(in) :: FSNOWC(ILG)  !< Fractional coverage of canopy by frozen water
  !! over snow-free subarea [ ]
  real, intent(in) :: FSNOCS(ILG)  !< Fractional coverage of canopy by frozen water
  !! over snow-covered subarea [ ]
  real, intent(in) :: PSIGND(ILG)  !< Minimum liquid moisture suction in soil layers \f$[m] (\Psi_s)\f$
  real, intent(in) :: FCLOUD(ILG)  !< Fractional cloud cover \f$[ ] (X_i) \f$
  real, intent(in) :: COSZS(ILG)  !< Cosine of solar zenith angle [ ]
  real, intent(in) :: QSWINV(ILG)  !< Visible radiation incident on horizontal
  !! surface \f$[W m^{-2}] (K)\f$
  real, intent(in) :: VPD(ILG)  !< Vapour pressure deficit of air \f$[mb] {\Delta e }\f$
  real, intent(in) :: TA(ILG)  !< Air temperature at reference height \f$[K] (T_a) \f$
  real, intent(in) :: ALVSGC(ILG)  !< Visible/near-IR albedo of bare ground under vegetation [ ]
  real, intent(in) :: ALIRGC(ILG)  !< Visible/near-IR albedo of bare ground under vegetation [ ]
  real, intent(in) :: ALVSSC(ILG)  !< Visible/near-IR albedo of snow under vegetation [ ]
  real, intent(in) :: ALIRSC(ILG)  !< Visible/near-IR albedo of snow under vegetation [ ]
  !
  !     * WORK ARRAYS.
  !
  real, intent(inout) :: CXTEFF(ILG,IC), RCACC(ILG,IC), &
                         RCV(ILG,IC), RCG(ILG,IC), &
                         RCT(ILG), GC(ILG), &
                         TRVS(ILG), TRIR(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: SVF, ALVSCX, ALIRCX, ALVSN, ALIRN, ALVSS, ALIRS, &
          TRTOT, EXPMAX1, EXPMAX2, EXPMAX3, TMP, TRCLRV, &
          TRCLDV, TRCLRT, TRCLDT
  !
  !----------------------------------------------------------------------
  !>
  !! The transmissivity \f$\tau_c\f$ of a vegetation canopy to shortwave
  !! radiation is obtained by applying a form of Beer’s law of
  !! radiation transfer in non-scattering media:
  !! \f$\tau_c = exp[-\kappa \Lambda_p]\f$
  !! where \f$\kappa\f$ is the canopy extinction coefficient and \f$\Lambda_p\f$ is the
  !! plant area index. The extinction coefficient is calculated after
  !! Goudriaan (1988) \cite Gold1958-ng as
  !! \f$\kappa = \epsilon O/cos(Z)\f$
  !! where epsilon is a correction factor less than or equal to 1,
  !! accounting for forward-scattering of radiation and non-random
  !! leaf distributions (i.e. clumping), O represents the mean
  !! projected leaf area fraction perpendicular to the incoming
  !! radiation, and Z is the zenith angle of the incoming radiation.
  !! For crops, grass and needleleaf trees, the distribution of leaf
  !! angles is assumed to be spherical, and thus O = 0.5. For
  !! broadleaf trees, the preferred leaf orientation tends to be
  !! horizontal; thus O = cos(Z).
  !!
  !! In the clear-sky case, incoming shortwave radiation is dominated
  !! by the direct beam; thus, the transmissivity for clear skies
  !! \f$\tau_{c, 0}\f$ is evaluated by simply setting Z to \f$Z_s\f$, the solar zenith
  !! angle. An extensive search of the literature for values of
  !! epsilon appropriate to the four vegetation categories yielded the
  !! following results for the extinction coefficient:
  !!
  !! \f$\kappa = 0.3/cos(Z_s)  (needleleaf trees)\f$
  !! \f$\kappa = 0.4          (broadleaf trees, full canopy)\f$
  !! \f$\kappa = 0.4/cos(Z_s)  (broadleaf trees, leafless)\f$
  !! \f$\kappa = 0.4/cos(Z_s)  (crops and grass)\f$
  !!
  !! In the visible range of the shortwave spectrum, scattering is less important because of high leaf
  !! absorptivities. The following results were obtained for visible radiation:
  !!
  !! \f$\kappa = 0.4/cosZ_s     (needleleaf trees)\f$
  !! \f$\kappa = 0.7            (broadleaf trees, full canopy)\f$
  !! \f$\kappa = 0.4/cosZ_s     (broadleaf trees, leafless)\f$
  !! \f$\kappa = 0.5/cosZ_s     (crops and grass)\f$
  !!
  !! In the case of overcast skies, the hemispherical distribution of
  !! incoming shortwave radiation is modelled using the generally
  !! accepted “standard overcast sky” distribution (e.g. Steven and
  !! Unsworth, 1980), where the shortwave radiation D(Z) emanating
  !! from a sky zenith angle Z is approximated as
  !!
  !! \f$D(Z) = D(0) [(1 + 1.23cosZ)/2.23].\f$
  !!
  !! Integration of the cloudy-sky transmissivity \f$tau_{c, cloudy}\f$ over
  !! the sky hemisphere is performed using a simple weighting function
  !! proposed by Goudriaan (1988) \cite Gold1958-ng :
  !!
  !! \f$\tau_{c, cloudy} = 0.3 \tau_c(Z=15^o) + 0.5 \tau_c(Z=45^o) + 0.2 \tau_c(Z=75^o)\f$
  !!
  !! The albedo \f$\alpha_c\f$ of a vegetation canopy is, like the
  !! transmissivity, dependent in principle on the zenith angle Z of
  !! the incoming radiation and the mean projected leaf area fraction
  !! O. However, in practice the observed diurnal and seasonal
  !! variation of vegetation albedo tends to be slight for closed
  !! canopies or for radiation zenith angles larger than a few
  !! degrees. Since the former is generally true in the tropics and
  !! the latter in the extratropics, the diurnal variation of visible
  !! and near-infrared canopy albedo is neglected in CLASS.
  !!
  !! Corrections are applied to the average albedos to account for the
  !! effects of intercepted snow on the canopy and incomplete canopy
  !! closure. The presence of snow on the canopy will increase the
  !! albedo; this increase will be larger in the visible range of the
  !! spectrum, since the visible albedo of vegetation is typically
  !! very small. Making use of data presented in Bartlett et al. (2015) \cite Bartlett2015-gw,
  !! the visible albedo of snow-covered vegetation is set to
  !! 0.27 and near-infrared albedo to 0.38.
  !!
  !! If the canopy closure is incomplete, a fraction chi of the ground
  !! or snow under the canopy will be visible through gaps in it. This
  !! fraction is equal to the sky view factor of the underlying
  !! surface, which is calculated using an equation analogous to that
  !! for the canopy transmissivity, as an exponential function of the
  !! plant area index \f$\Lambda_p\f$:
  !!
  !! \f$\chi = exp[-c\Lambda_p]\f$
  !!
  !! where c is a constant depending on the vegetation category. The
  !! overall albedo \f$\alpha_T\f$ is calculated as a weighted average of the
  !! canopy albedo and the underlying ground or snow albedo (the
  !! latter weighted by the canopy transmissivity, to account for the
  !! decreased downwelling shortwave radiation at the surface):
  !!
  !! \f$\alpha_T = (1 – \chi) \alpha_c + \chi \tau_c \alpha_0\f$
  !!
  !! where \f$alpha_0\f$ is the albedo of the surface under the canopy.
  !!
  !
  !     * ASSIGN CONSTANT EXPONENTIATION TERMS: EXPMAX1=EXP(-0.4/0.9659),
  !     * EXPMAX2=EXP(-0.4/0.7071), EXPMAX3=EXP(-0.4/0.2588)
  !
  EXPMAX1 = 0.6609  
  EXPMAX2 = 0.5680  
  EXPMAX3 = 0.2132  
  !!
  !> At the beginning of the subroutine, values are assigned to
  !! exponentiation terms and a series of work arrays is
  !! initialized to zero. Then the transmissivity and albedo of each
  !! vegetation category are calculated in turn, first over bare soil
  !! and then over a snow pack. After each of the latter sets of
  !! calculations, consistency checks are carried out to ensure that
  !! the calculated albedos and transmissivities are not less than 0
  !! or greater than 1, and that the transmissivity is not greater
  !! than 90% of the non-reflected radiation.
  !!
  !
  !     * INITIALIZE WORK ARRAYS.
  !
  do I = IL1,IL2 ! loop 50
    RCT(I) = 0.0
    GC(I) = 0.0
    RC(I) = 0.0
  end do ! loop 50
  do J = 1,IC ! loop 100
    do I = IL1,IL2
      CXTEFF(I,J) = 0.0
      RCACC(I,J) = 0.0
      RCG(I,J) = 0.0
      RCV(I,J) = 0.0
    end do
  end do ! loop 100
  !!
  !> The clear-sky and cloudy-sky transmissivities in the visible
  !! range are calculated first, and then the clear-sky and cloudy-sky
  !! transmissivities in the total shortwave spectrum, using the
  !! extinction coefficient equations above. (For broadleaf trees, the
  !! minimum of the full-canopy and leafless calculated values is
  !! used.) The overall transmissivity in the visible range at the
  !! current time step is obtained as an average of the clear and
  !! cloudy-sky values, weighted according to the fractional cloud
  !! cover. The effective overall extinction coefficient CXTEFF in the
  !! visible range is computed for use later on in the stomatal
  !! resistance calculation. The overall transmissivity for the entire
  !! shortwave spectrum at the current time step is similarly obtained
  !! as an average of the clear and cloudy-sky values, weighted
  !! according to the fractional cloud cover. The near-infrared
  !! transmissivity is calculated from the total and visible
  !! transmissivities, making use of the assumption that visible and
  !! near-infrared radiation each typically comprise 50% of the
  !! incoming shortwave radiation. Lastly, the aggregated visible and
  !! near-infrared transmissivities for the bulk canopy are
  !! incremented using the current values weighted by the fractional
  !! coverage of the vegetation category. At the end, the sum is
  !! normalized by the total fractional vegetation coverage in the
  !! subarea.
  !!
  !! For the albedos, first the sky view factor and the near-infrared
  !! albedo of snow-covered vegetation are calculated. Next, a branch
  !! directs further processing depending on the value of the flag
  !! IALC. If IALC = 0, the CLASS-calculated visible and near-infrared
  !! vegetation albedos are used. These albedos are corrected for the
  !! presence of intercepted snow, using an average of the snow-
  !! covered and background snow-free albedos weighted according to
  !! FSNOWC over snow-free, and FSNOCS over snow-covered subareas
  !! (representing the ratio of the amount of snow present on the
  !! canopy relative to the interception capacity). The correction for
  !! incomplete canopy closure is applied to the resulting visible and
  !! near- infrared albedos as described above. If IALC = 1, user-
  !! specified albedos for the vegetation canopy are used instead of
  !! the CLASS values. These are assumed to incorporate the effects of
  !! incomplete canopy closure, but not of the presence of snow. Thus
  !! the CLASS values for the albedo of snow-covered vegetation and
  !! the albedo of snow under the canopy are still used in the
  !! averaging. Finally, the aggregated visible and near-infrared
  !! albedos for the bulk canopy are incremented using the current
  !! values weighted by the fractional coverage of the vegetation
  !! category. At the end, the sum is normalized by the total
  !! fractional vegetation coverage in the subarea.
  !!
  !
  !     * ALBEDO AND TRANSMISSIVITY CALCULATIONS FOR CANOPY OVER
  !     * BARE SOIL.

  do J = 1,IC
    select case (classpfts(J))
    case ('NdlTr') !-------------------------------------------------------
      do I = IL1,IL2 ! loop 150
        if (COSZS(I) > 0. .and. FCAN(I,J) > 0.) then
          TRCLRV = EXP( - 0.4 * PAI(I,J) / COSZS(I))
          !              TMP=MAX(-50.0, -0.4*PAI(I,J)/COSZS(I))  ! JM EDIT
          !              TRCLRV=EXP(TMP)
          TRCLDV = 0.30 * EXP( - 0.4 * PAI(I,J) / 0.9659) + 0.50 * EXP( - 0.4 * &            
                   PAI(I,J) / 0.7071) + 0.20 * EXP( - 0.4 * PAI(I,J) / 0.2588)
          TRCLRT = EXP( - 0.3 * PAI(I,J) / COSZS(I))
          !              TMP=MAX(-50.0,(-0.3*PAI(I,J)/COSZS(I)))    ! JM EDIT
          !              TRCLRT = EXP(TMP)
          TRCLDT = 0.30 * EXP( - 0.3 * PAI(I,J) / 0.9659) + 0.50 * EXP( - 0.3 * &
                   PAI(I,J) / 0.7071) + 0.20 * EXP( - 0.3 * PAI(I,J) / 0.2588)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          if (TRVS(I) > 0.0001) then
            CXTEFF(I,J) = - LOG(TRVS(I)) / MAX(PAI(I,J),1.0E-5)
          else
            CXTEFF(I,J) = CXTLRG
          end if
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCN(I) = TRVSCN(I) + FCAN(I,J) * TRVS(I)
          TRIRCN(I) = TRIRCN(I) + FCAN(I,J) * TRIR(I)

          SVF = EXP(CANEXT(J) * PAI(I,J))
          if (IALC == 0) then
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ALVSC(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ALIRC(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSGC(I)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRGC(I)
          else
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ACVDAT(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ACIDAT(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * ACVDAT(I,J)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * ACIDAT(I,J)
          end if
          ALVSCN(I) = ALVSCN(I) + FCAN(I,J) * ALVSN
          ALIRCN(I) = ALIRCN(I) + FCAN(I,J) * ALIRN
        end if
      end do ! loop 150
      !
    case ('BdlTr','BdlSh') !-------------------------------------------------------
      do I = IL1,IL2 ! loop 250
        if (COSZS(I) > 0. .and. FCAN(I,J) > 0.) then
          TRCLRV = MIN(EXP( - 0.7 * PAI(I,J)),EXP( - 0.4 / COSZS(I)))
          TRCLDV = 0.30 * MIN(EXP( - 0.7 * PAI(I,J)),EXPMAX1) &            
                   + 0.50 * MIN(EXP( - 0.7 * PAI(I,J)),EXPMAX2) &
                   + 0.20 * MIN(EXP( - 0.7 * PAI(I,J)),EXPMAX3)
          TRCLRT = MIN(EXP( - 0.4 * PAI(I,J)),EXP( - 0.4 / COSZS(I)))
          TRCLDT = 0.30 * MIN(EXP( - 0.4 * PAI(I,J)),EXPMAX1) + &
                   0.50 * MIN(EXP( - 0.4 * PAI(I,J)),EXPMAX2) + &
                   0.20 * MIN(EXP( - 0.4 * PAI(I,J)),EXPMAX3)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          if (TRVS(I) > 0.0001) then
            CXTEFF(I,J) = - LOG(TRVS(I)) / MAX(PAI(I,J),1.0E-5)
          else
            CXTEFF(I,J) = CXTLRG
          end if
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCN(I) = TRVSCN(I) + FCAN(I,J) * TRVS(I)
          TRIRCN(I) = TRIRCN(I) + FCAN(I,J) * TRIR(I)

          SVF = EXP(CANEXT(J) * PAI(I,J))
          if (IALC == 0) then
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ALVSC(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ALIRC(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSGC(I)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRGC(I)
          else
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ACVDAT(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ACIDAT(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * ACVDAT(I,J)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * ACIDAT(I,J)
          end if
          ALVSCN(I) = ALVSCN(I) + FCAN(I,J) * ALVSN
          ALIRCN(I) = ALIRCN(I) + FCAN(I,J) * ALIRN
        end if
      end do ! loop 250
      !
    case ('Crops', 'Grass') ! CROPS AND GRASS
      do I = IL1,IL2 ! loop 350
        if (COSZS(I) > 0. .and. FCAN(I,J) > 0.) then
          TRCLRV = EXP( - 0.5 * PAI(I,J) / COSZS(I))
          !              TMP=MAX(-50.0, -0.5*PAI(I,J)/COSZS(I))  ! JM EDIT
          !              TRCLRV=EXP(TMP)
          TRCLDV = 0.30 * EXP( - 0.5 * PAI(I,J) / 0.9659) + 0.50 * EXP( - 0.5 * &
                   PAI(I,J) / 0.7071) + 0.20 * EXP( - 0.5 * PAI(I,J) / 0.2588)
          TRCLRT = EXP( - 0.4 * PAI(I,J) / COSZS(I))
          !              TMP=MAX(-50.0,(-0.4*PAI(I,J)/COSZS(I)))    ! JM EDIT
          !              TRCLRT = EXP(TMP)
          TRCLDT = 0.30 * EXP( - 0.4 * PAI(I,J) / 0.9659) + 0.50 * EXP( - 0.4 * &
                   PAI(I,J) / 0.7071) + 0.20 * EXP( - 0.4 * PAI(I,J) / 0.2588)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          if (TRVS(I) > 0.0001) then
            CXTEFF(I,J) = - LOG(TRVS(I)) / MAX(PAI(I,J),1.0E-5)
          else
            CXTEFF(I,J) = CXTLRG
          end if
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCN(I) = TRVSCN(I) + FCAN(I,J) * TRVS(I)
          TRIRCN(I) = TRIRCN(I) + FCAN(I,J) * TRIR(I)

          SVF = EXP(CANEXT(J) * PAI(I,J))
          if (IALC == 0) then
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ALVSC(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ALIRC(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSGC(I)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRGC(I)
          else ! user-specified values read-in.
            ALVSCX = FSNOWC(I) * ALVSWC + (1.0 - FSNOWC(I)) * ACVDAT(I,J)
            ALIRCX = FSNOWC(I) * ALIRWC + (1.0 - FSNOWC(I)) * ACIDAT(I,J)
            ALVSN = (1.0 - SVF) * ALVSCX + SVF * ACVDAT(I,J)
            ALIRN = (1.0 - SVF) * ALIRCX + SVF * ACIDAT(I,J)
          end if
          ALVSCN(I) = ALVSCN(I) + FCAN(I,J) * ALVSN
          ALIRCN(I) = ALIRCN(I) + FCAN(I,J) * ALIRN
        end if
      end do ! loop 350

    case default
      print * ,'Unknown CLASS PFT in canopyAlbedoTransmiss ',classpfts(J)
      call errorHandler('canopyAlbedoTransmiss', - 1)
    end select
  end do
  !
  !     * TOTAL ALBEDOS.
  !
  IPTBAD = 0
  do I = IL1,IL2 ! loop 450
    if (FC(I) > 0. .and. COSZS(I) > 0.) then
      ALVSCN(I) = ALVSCN(I) / FC(I)
      ALIRCN(I) = ALIRCN(I) / FC(I)
    end if
    if (ALVSCN(I) > 1. .or. ALVSCN(I) < 0.) IPTBAD = I
    if (ALIRCN(I) > 1. .or. ALIRCN(I) < 0.) IPTBAD = I
  end do ! loop 450
  !
  if (IPTBAD /= 0) then
    write(6,6100) IPTBAD,JL,ALVSCN(IPTBAD),ALIRCN(IPTBAD)
6100 format('0AT (I,J) = (',I3,',',I3,'),ALVSCN,ALIRCN = ',2F10.5)
    call errorHandler('canopyAlbedoTransmiss', - 2)
  end if
  !
  !     * TOTAL TRANSMISSIVITIES.
  !
  IPTBAD = 0
  do I = IL1,IL2 ! loop 475
    if (FC(I) > 0. .and. COSZS(I) > 0.) then
      TRVSCN(I) = TRVSCN(I) / FC(I)
      TRIRCN(I) = TRIRCN(I) / FC(I)
      TRVSCN(I) = MIN( TRVSCN(I),0.90 * (1.0 - ALVSCN(I)) )  
      TRIRCN(I) = MIN( TRIRCN(I),0.90 * (1.0 - ALIRCN(I)) )
    end if
    if (TRVSCN(I) > 1. .or. TRVSCN(I) < 0.) IPTBAD = I
    if (TRIRCN(I) > 1. .or. TRIRCN(I) < 0.) IPTBAD = I
  end do ! loop 475
  !
  if (IPTBAD /= 0) then
    write(6,6300) IPTBAD,JL,TRVSCN(IPTBAD),TRIRCN(IPTBAD)
6300 format('0AT (I,J) = (',I3,',',I3,'),TRVSCN,TRIRCN = ',2F10.5)
    call errorHandler('canopyAlbedoTransmiss', - 3)
  end if
  !----------------------------------------------------------------------
  !     * ALBEDO AND TRANSMISSIVITY CALCULATIONS FOR CANOPY OVER SNOW.
  do J = 1,IC
    select case (classpfts(J))
    case ('NdlTr') !-------------------------------------------------------
      do I = IL1,IL2 ! loop 500
        if (COSZS(I) > 0. .and. FCANS(I,J) > 0.) then
          TRCLRV = EXP( - 0.4 * PAIS(I,J) / COSZS(I))   
          !              TMP=MAX(-50.0, -0.4*PAIS(I,J)/COSZS(I))  ! JM EDIT
          !              TRCLRV=EXP(TMP)
          TRCLDV = 0.30 * EXP( - 0.4 * PAIS(I,J) / 0.9659) + 0.50 * EXP( - 0.4 * &
                   PAIS(I,J) / 0.7071) + 0.20 * EXP( - 0.4 * PAIS(I,J) / 0.2588)   
          TRCLRT = EXP( - 0.3 * PAIS(I,J) / COSZS(I))
          !              TMP=MAX(-50.0,(-0.3*PAIS(I,J)/COSZS(I)))    ! JM EDIT
          !              TRCLRT = EXP(TMP)
          TRCLDT = 0.30 * EXP( - 0.3 * PAIS(I,J) / 0.9659) + 0.50 * EXP( - 0.3 * &
                   PAIS(I,J) / 0.7071) + 0.20 * EXP( - 0.3 * PAIS(I,J) / 0.2588)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCS(I) = TRVSCS(I) + FCANS(I,J) * TRVS(I)
          TRIRCS(I) = TRIRCS(I) + FCANS(I,J) * TRIR(I)
          !          END IF
          !  500 CONTINUE
          !
          !      do I=IL1,IL2
          !          IF (COSZS(I)>0. .AND. FCANS(I,J)>0.)             THEN
          if (IALC == 0) then
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ALVSC(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ALIRC(I,J)
          else
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ACVDAT(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ACIDAT(I,J)
          end if
          SVF = EXP(CANEXT(J) * PAIS(I,J))
          ALVSS = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSSC(I)
          ALIRS = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRSC(I)
          ALVSCS(I) = ALVSCS(I) + FCANS(I,J) * ALVSS
          ALIRCS(I) = ALIRCS(I) + FCANS(I,J) * ALIRS
        end if
        !  550 CONTINUE
      end do ! loop 500
      !
    case ('BdlTr','BdlSh') !-------------------------------------------------------
      do I = IL1,IL2 ! loop 600
        if (COSZS(I) > 0. .and. FCANS(I,J) > 0.) then
          TRCLRV = MIN(EXP( - 0.7 * PAIS(I,J)),EXP( - 0.4 / COSZS(I)))
          TRCLDV = 0.30 * MIN(EXP( - 0.7 * PAIS(I,J)),EXPMAX1) &
                   + 0.50 * MIN(EXP( - 0.7 * PAIS(I,J)),EXPMAX2) &
                   + 0.20 * MIN(EXP( - 0.7 * PAIS(I,J)),EXPMAX3)
          TRCLRT = MIN(EXP( - 0.4 * PAIS(I,J)),EXP( - 0.4 / COSZS(I)))
          TRCLDT = 0.30 * MIN(EXP( - 0.4 * PAIS(I,J)),EXPMAX1) + &
                   0.50 * MIN(EXP( - 0.4 * PAIS(I,J)),EXPMAX2) + &
                   0.20 * MIN(EXP( - 0.4 * PAIS(I,J)),EXPMAX3)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCS(I) = TRVSCS(I) + FCANS(I,J) * TRVS(I)
          TRIRCS(I) = TRIRCS(I) + FCANS(I,J) * TRIR(I)
          !          END IF
          !  600 CONTINUE
          !
          !      do I=IL1,IL2
          !          IF (COSZS(I)>0. .AND. FCANS(I,J)>0.)             THEN
          if (IALC == 0) then
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ALVSC(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ALIRC(I,J)
          else
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ACVDAT(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ACIDAT(I,J)
          end if
          SVF = EXP(CANEXT(J) * PAIS(I,J))
          ALVSS = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSSC(I)
          ALIRS = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRSC(I)
          ALVSCS(I) = ALVSCS(I) + FCANS(I,J) * ALVSS
          ALIRCS(I) = ALIRCS(I) + FCANS(I,J) * ALIRS
        end if
        !  650 CONTINUE
      end do ! loop 600
      !
    case ('Crops', 'Grass')  ! CROPS AND GRASS.
      do I = IL1,IL2 ! loop 700
        if (COSZS(I) > 0. .and. FCANS(I,J) > 0.) then
          TRCLRV = EXP( - 0.5 * PAIS(I,J) / COSZS(I))
          !              TMP=MAX(-50.0, -0.5*PAIS(I,J)/COSZS(I))  ! JM EDIT
          !              TRCLRV=EXP(TMP)
          TRCLDV = 0.30 * EXP( - 0.5 * PAIS(I,J) / 0.9659) + 0.50 * EXP( - 0.5 * &
                   PAIS(I,J) / 0.7071) + 0.20 * EXP( - 0.5 * PAIS(I,J) / 0.2588)
          TRCLRT = EXP( - 0.4 * PAIS(I,J) / COSZS(I))
          !              TMP=MAX(-50.0,(-0.4*PAIS(I,J)/COSZS(I)))    ! JM EDIT
          !              TRCLRT = EXP(TMP)
          TRCLDT = 0.30 * EXP( - 0.4 * PAIS(I,J) / 0.9659) + 0.50 * EXP( - 0.4 * &
                   PAIS(I,J) / 0.7071) + 0.20 * EXP( - 0.4 * PAIS(I,J) / 0.2588)
          TRVS(I) = FCLOUD(I) * TRCLDV + (1.0 - FCLOUD(I)) * TRCLRV
          TRTOT = FCLOUD(I) * TRCLDT + (1.0 - FCLOUD(I)) * TRCLRT
          TRIR(I) = 2. * TRTOT - TRVS(I)
          TRVSCS(I) = TRVSCS(I) + FCANS(I,J) * TRVS(I)
          TRIRCS(I) = TRIRCS(I) + FCANS(I,J) * TRIR(I)
          !          END IF
          !  700 CONTINUE

          !      do I=IL1,IL2
          !          IF (COSZS(I)>0. .AND. FCANS(I,J)>0.)             THEN
          if (IALC == 0) then
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ALVSC(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ALIRC(I,J)
          else
            ALVSCX = FSNOCS(I) * ALVSWC + (1.0 - FSNOCS(I)) * ACVDAT(I,J)
            ALIRCX = FSNOCS(I) * ALIRWC + (1.0 - FSNOCS(I)) * ACIDAT(I,J)
          end if
          SVF = EXP(CANEXT(J) * PAIS(I,J))
          ALVSS = (1.0 - SVF) * ALVSCX + SVF * TRVS(I) * ALVSSC(I)
          ALIRS = (1.0 - SVF) * ALIRCX + SVF * TRIR(I) * ALIRSC(I)
          ALVSCS(I) = ALVSCS(I) + FCANS(I,J) * ALVSS
          ALIRCS(I) = ALIRCS(I) + FCANS(I,J) * ALIRS
        end if
        !  750 CONTINUE
      end do ! loop 700
    case default
      print * ,'Unknown CLASS PFT in canopyAlbedoTransmiss ',classpfts(J)
      call errorHandler('canopyAlbedoTransmiss', - 4)
    end select

  end do
  !
  !     * TOTAL ALBEDOS AND CONSISTENCY CHECKS.
  !
  IPTBAD = 0
  do I = IL1,IL2 ! loop 775
    if (FCS(I) > 0. .and. COSZS(I) > 0.) then
      ALVSCS(I) = ALVSCS(I) / FCS(I)
      ALIRCS(I) = ALIRCS(I) / FCS(I)
    end if
    if (ALVSCS(I) > 1. .or. ALVSCS(I) < 0.) IPTBAD = I
    if (ALIRCS(I) > 1. .or. ALIRCS(I) < 0.) IPTBAD = I
  end do ! loop 775
  !
  !     * TOTAL TRANSMISSIVITIES AND CONSISTENCY CHECKS.
  !
  IPTBAD = 0
  JPTBAD = 0
  do I = IL1,IL2 ! loop 800
    if (FCS(I) > 0. .and. COSZS(I) > 0.) then
      TRVSCS(I) = TRVSCS(I) / FCS(I)
      TRIRCS(I) = TRIRCS(I) / FCS(I)
      TRVSCS(I) = MIN( TRVSCS(I),0.90 * (1.0 - ALVSCS(I)) )
      TRIRCS(I) = MIN( TRIRCS(I),0.90 * (1.0 - ALIRCS(I)) )
    end if
    if (TRVSCS(I) > 1. .or. TRVSCS(I) < 0.) IPTBAD = I
    if (TRIRCS(I) > 1. .or. TRIRCS(I) < 0.) IPTBAD = I
    if ((1. - ALVSCN(I) - TRVSCN(I)) < 0.) then
      JPTBAD = 1000 + I
      JPTBDI = I
    end if
    if ((1. - ALVSCS(I) - TRVSCS(I)) < 0.) then
      JPTBAD = 2000 + I
      JPTBDI = I
    end if
    if ((1. - ALIRCN(I) - TRIRCN(I)) < 0.) then
      JPTBAD = 3000 + I
      JPTBDI = I
    end if
    if ((1. - ALIRCS(I) - TRIRCS(I)) < 0.) then
      JPTBAD = 4000 + I
      JPTBDI = I
    end if
  end do ! loop 800
  !
  if (IPTBAD /= 0) then
    write(6,6400) IPTBAD,JL,TRVSCS(IPTBAD),TRIRCS(IPTBAD)
6400 format('0AT (I,J) = (',I3,',',I3,'),TRVSCS,TRIRCS = ',2F10.5)
    call errorHandler('canopyAlbedoTransmiss', - 5)
  end if
  !
  if (IPTBAD /= 0) then
    write(6,6200) IPTBAD,JL,ALVSCS(IPTBAD),ALIRCS(IPTBAD)
6200 format('0AT (I,J) = (',I3,',',I3,'),ALVSCS,ALIRCS = ',2F10.5)
    call errorHandler('canopyAlbedoTransmiss', - 6)
  end if
  !
  if (JPTBAD /= 0) then
    write(6,6500) JPTBDI,JL,JPTBAD
6500 format('0AT (I,J) = (',I3,',',I3,'),JPTBAD =  ',I5)
    call errorHandler('canopyAlbedoTransmiss', - 7)
  end if
  !-----------------------------------------------------------------------
  !>
  !! In the final section, the stomatal resistance \f$r_c\f$ of the
  !! vegetation canopy is determined. Based on the analysis of Schulze
  !! et al. (1995), the unstressed stomatal resistance \f$r_{c, u}\f$ for a
  !! given vegetation category can be calculated as a function of the
  !! incoming visible shortwave radiation \f$K\downarrow\f$:
  !!
  !! \f$r_{c, u} = r_{s, min} \kappa_e / ln[{K\downarrow + K\downarrow_{1/2} / \kappa_e} / { K\downarrow exp(-\kappa_e \Lambda) + K\downarrow_{1/2} / \kappa_e }]\f$
  !!
  !! where \f$r_{s, min}\f$ is the minimum stomatal resistance for the
  !! vegetation category, \f$\kappa_e\f$ is the extinction coeffient for
  !! visible radiation (CXTEFF above), \f$\Lambda\f$ is the leaf area index, and
  !! \f$K\downarrow_{1/2}\f$ is the value of \f$K\downarrow\f$ at which \f$r_{c, u} = 2 r_{s, min}.\f$
  !!
  !! Suboptimum environmental conditions for transpiration may lead to
  !! stresses on the plant, causing the stomatal resistance to be
  !! greater than its unstressed value. The effects of these stresses
  !! are modelled by defining functions of the air temperature \f$T_a\f$, the
  !! air vapour pressure deficit \f$\Delta e\f$, and the soil moisture suction
  !! \f$\Psi_s\f$. These functions are used to derive \f$r_{c, i}\f$ of each
  !! vegetation category on the basis of \f$r_{c, u, i}\f$:
  !!
  !! \f$r_{c, i} = f(T_a) f( \Delta_e) f( \Psi_s) * r_{c, u, i}\f$
  !!
  !! The air temperature function \f$f(T_a)\f$ has a value of 1 for
  !! temperatures between \f$5^o C\f$ and \f$40^o C\f$, and an
  !! arbitrary large value of 250 for temperatures less than \f$-5^o C\f$
  !! and greater than \f$50^o C\f$. Between these points it
  !! varies in a linear fashion. For the vapour pressure deficit
  !! function \f$f( \Delta e)\f$, two alternate forms are provided, after Oren et
  !! al. (1999) and Wu et al. (2000) respectively:
  !!
  !! \f$f( \Delta e) = [( \Delta e/10.0)^{cv2}]/c_{v1}\f$ and
  !! \f$f(\Delta e) = 1/[exp(-c_{v1} \Delta e/10.0)]\f$
  !!
  !! where \f$c_{v1}\f$ and \f$c_{v2}\f$ are parameters depending on the vegetation
  !! category. If \f$c_{v2}\f$ is greater than zero, the first form is used
  !! if not, the second form is used. The soil moisture suction
  !! function \f$f(\Psi_s)\f$ is expressed, following Choudhury and Idso
  !! (1985) \cite Choudhury1985-mm and Fisher et al. (1981) \cite Fisher1981-xf, as:
  !!
  !! \f$f(\Psi_s) = 1 + (\Psi_s / c_{\Psi 1})^{c \Psi 2}\f$
  !!
  !! where \f$c_{\Psi 1}\f$ and \f$c_{\Psi 2}\f$ are parameters depending on the vegetation
  !! category. Finally, the aggregated stomatal resistance for the
  !! canopy over the bare ground subarea is obtained as a weighted
  !! average over the vegetation categories. (It is assumed that
  !! transpiration is suppressed when snow is present under the
  !! canopy, so rc for this subarea is set to a large number). Since,
  !! following the electrical analogy, resistances act in parallel,
  !! the aggregated resistance for the subarea of canopy over bare
  !! ground is obtained as an average of inverses:
  !!
  !! \f$X/r_c = \Sigma (X_i / r_{c, i})\f$
  !!
  !! The calculations described above pertaining to stomatal
  !! resistances are performed in loops 850, 900 and 950. In the 850
  !! loop, \f$f(T_a)\f$ is evaluated. In the 900 loop, for each vegetation
  !! category in turn, \f$f(\Delta e)\f$, \f$f(\Psi_s)\f$ and \f$r_{c, i}\f$ are determined. \f$r_{c, i}\f$
  !! is assigned upper and lower limits of 5000 and 10 \f$s m^{-1}\f$
  !! respectively, and the accumulated stomatal resistance for the
  !! canopy is incremented by \f$X_i / r_{c, i}\f$. In loop 950, if the incoming
  !! visible radiation is small, the stomatal resistances are likewise
  !! set to a large number. Otherwise, the stomatal resistance of
  !! vegetation over snow is set to a large number, and the
  !! normalized, aggregated stomatal resistance of vegetation over
  !! soil is obtained by inverting the accumulated value.
  !!
  !
  !     * BULK STOMATAL RESISTANCES FOR CANOPY OVERLYING SNOW AND CANOPY
  !     * OVERLYING BARE SOIL.
  !
  do I = IL1,IL2 ! loop 850
    if ((FCS(I) + FC(I)) > 0.0) then
      if (TA(I) <= 268.15) then  
        RCT(I) = 250.   
      else if (TA(I) < 278.15) then  
        RCT(I) = 1. / (1. - (278.15 - TA(I)) * .1)
      else if (TA(I) > 313.15) then  
        if (TA(I) >= 323.15) then
          RCT(I) = 250.   
        else
          RCT(I) = 1. / (1. - (TA(I) - 313.15) * 0.1)  
        end if
      else
        RCT(I) = 1.
      end if
    end if
  end do ! loop 850
  !
  do J = 1,IC ! loop 900
    do I = IL1,IL2
      if (FCAN(I,J) > 0.) then
        if (VPD(I) > 0. .and. VPDA(I,J) > 0.0) then
          if (ABS(VPDB(I,J)) > 1.0E-5) then
            RCV(I,J) = MAX(1.,((VPD(I) / 10.) ** VPDB(I,J)) / &
                       VPDA(I,J))
          else
            RCV(I,J) = 1. / EXP( - VPDA(I,J) * VPD(I) / 10.)
          end if
        else
          RCV(I,J) = 1.0
        end if
        if (PSIGA(I,J) > 0.0) then
          RCG(I,J) = 1. + (PSIGND(I) / PSIGA(I,J)) ** PSIGB(I,J)
        else
          RCG(I,J) = 1.0
        end if
        if (QSWINV(I) > 0.01 .and. COSZS(I) > 0. .and. &
            CXTEFF(I,J) > 1.0E-5 .and. RCG(I,J) < 1.0E5) then
          RCACC(I,J) = MIN(CXTEFF(I,J) * RSMIN(I,J) / LOG((QSWINV(I) + &
                       QA50(I,J) / CXTEFF(I,J)) / (QSWINV(I) * EXP( - CXTEFF(I,J) * &
                       PAI(I,J)) + QA50(I,J) / CXTEFF(I,J))) * RCV(I,J) * RCG(I,J) * &
                       RCT(I),5000.)
          RCACC(I,J) = MAX(RCACC(I,J),10.0)
        else
          RCACC(I,J) = 5000.
        end if
        RC(I) = RC(I) + FCAN(I,J) / RCACC(I,J)
      end if
    end do
  end do ! loop 900
  !
  do I = IL1,IL2 ! loop 950
    if ((FCS(I) + FC(I)) > 0.) then
      if (QSWINV(I) < 2.0) then
        RCS(I) = 5000.0
        RC(I) = 5000.0
      else
        RCS(I) = 5000.0
        if (RC(I) > 0.) then
          RC(I) = FC(I) / RC(I)
        else
          RC(I) = 5000.0
        end if
      end if
    else
      RC(I) = 0.0
      RCS(I) = 0.0
    end if
  end do ! loop 950
  !
  return
end subroutine canopyAlbedoTransmiss
