!> \file
!> Calculates visible and near-IR ground albedos.
!! @author D. Verseghy, M. Lazare
!
subroutine groundAlbedo (ALVSG, ALIRG, ALVSGC, ALIRGC, & ! Formerly GRALB
                         ALGWV, ALGWN, ALGDV, ALGDN, &
                         THLIQ, FSNOW, ALVSU, ALIRU, FCMXU, &
                         AGVDAT, AGIDAT, FG, ISAND, &
                         ILG, IG, IL1, IL2, JL, IALG)
  !
  !     * DEC 15/16 - D.VERSEGHY. ASSIGN ROCK ALBEDO USING SOIL COLOUR
  !     *                         INDEX INSTEAD OF VIA LOOKUP TABLE.
  !     * JAN 16/15 - D.VERSEGHY. CORRECT ACCOUNTING FOR URBAN ALBEDO.
  !     * FEB 09/15 - D.VERSEGHY. New version for gcm18 and class 3.6:
  !     *                         - Wet and dry albedoes for EACH of
  !     *                           visible and near-ir are passed in
  !     *                           instead of ALGWET and ALGDRY. These
  !     *                           are used to calculate ALISG and ALIRG.
  !     * NOV 16/13 - M.LAZARE.   FINAL VERSION FOR GCM17:
  !     *                         - REMOVE UNNECESSARY LOWER BOUND
  !     *                           OF 1.E-5 ON "FURB".
  !     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND
  !     *                         CANOPY-COVERED GROUND.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * SEP 04/03 - D.VERSEGHY. RATIONALIZE CALCULATION OF URBAN
  !     *                         ALBEDO.
  !     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF USER-
  !     *                         SPECIFIED VALUES TO GROUND ALBEDO.
  !     *                         PASS IN ICE AND ORGANIC ALBEDOS
  !     *                         VIA NEW COMMON BLOCK "CLASS8".
  !     * FEB 07/02 - D.VERSEGHY. REVISIONS TO BARE SOIL ALBEDO
  !     *                         CALCULATIONS; REMOVAL OF SOIL
  !     *                         MOISTURE EXTRAPOLATION TO SURFACE.
  !     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         CALCULATE SOIL ALBEDO FROM PERCENT
  !     *                         SAND CONTENT RATHER THAN FROM COLOUR
  !     *                         INDEX.
  !     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         FIX BUG TO CALCULATE GROUND ALBEDO
  !     *                         UNDER CANOPIES AS WELL AS OVER BARE
  !     *                         SOIL.
  !     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.
  !     *                         "CALL ABORT" CHANGED TO "CALL errorHandler",
  !     *                         TO ENABLE RUNNING RUN ON PC'S.
  !     * FEB 12/93 - D.VERSEGHY/M.LAZARE. INCREASE DRY SOIL ALBEDO TO
  !     *                                  0.45 FROM 0.35.
  !     * MAR 13/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CLASS - VERSION 2.0.
  !     *                         CODE FOR MODEL VERSION GCM7U (WITH
  !     *                         CANOPY).
  !     * APR 11/89 - D.VERSEGHY. CALCULATE VISIBLE AND NEAR-IR SOIL
  !     *                         ALBEDOS BASED ON TEXTURE AND SURFACE
  !     *                         WETNESS. (SET TO ICE ALBEDOS OVER
  !     *                         CONTINENTAL ICE SHEETS.)
  !
  use classicParams, only : ALVSI, ALIRI, ALVSO, ALIRO

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IG, IL1, IL2, JL, IALG
  integer  :: IPTBAD, I
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: ALVSG (ILG)  !< Visible albedo of bare ground \f$[ ] (\alpha_{g, VIS})\f$
  real, intent(out) :: ALIRG (ILG)  !< Near-IR albedo of bare ground \f$[ ] (\alpha_{g, NIR})\f$
  real, intent(out) :: ALVSGC (ILG) !< Visible albedo of ground under vegetation canopy [ ]
  real, intent(out) :: ALIRGC (ILG) !< Near-IR albedo of ground under vegetation canopy [ ]
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: ALGWV (ILG)    !< Visible albedo of wet soil for modelled area \f$[ ] (alpha_{g, NIR, wet})\f$
  real, intent(in) :: ALGWN (ILG)    !< Near-infrared albedo of wet soil for modelled area  \f$[ ] (alpha_{g, NIR, wet})\f$
  real, intent(in) :: ALGDV (ILG)    !< Visible albedo of dry soil for modelled area \f$[ ] (alpha_{g, VIS, dry})\f$
  real, intent(in) :: ALGDN (ILG)    !< Near-infrared albedo of dry soil for modelled area  \f$[ ] (alpha_{g, NIR, dry})\f$
  real, intent(in) :: THLIQ (ILG,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(in) :: ALVSU (ILG)    !< Visible albedo of urban part of modelled area \f$[ ] (alpha_{u, VIS})\f$
  real, intent(in) :: ALIRU (ILG)    !< Near-IR albedo of urban part of modelled area \f$[ ] (alpha_{u, NIR})\f$
  real, intent(in) :: FCMXU (ILG)    !< Fractional coverage of urban part of modelled area \f$[ ] (X_u)\f$
  real, intent(in) :: AGVDAT(ILG)    !< Assigned value of visible albedo of ground – optional [ ]
  real, intent(in) :: AGIDAT(ILG)    !< Assigned value of near-IR albedo of ground – optional [ ]
  real, intent(in) :: FG    (ILG)    !< Fractional coverage of bare soil on modelled area \f$[ ] (X_g)\f$
  real, intent(in) :: FSNOW (ILG)    !< Fractional coverage of snow on modelled area  [  ]
  !
  integer, intent(in) :: ISAND (ILG, IG)!< Soil type flag based on sand content, assigned in subroutine soilProperties
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: FURB, ALBSOL
  !
  !---------------------------------------------------------------------
  !>
  !! If the ISAND flag for the surface soil layer is greater than or
  !! equal to zero (indicating mineral soil), first the urban area not
  !! covered by snow is evaluated.  Next the visible and near-IR open ground
  !! albedos, \f$\alpha_{g, VIS}\f$ and \f$\alpha_{g, NIR}\f$, for each wavelength range, are calculated on the basis of the
  !! wet and dry ground albedos \f$\alpha_{g, wet}\f$ and \f$\alpha_{g, dry}\f$ which were assigned
  !! for the modelled area in soilProperties. Idso et al. (1975) found a
  !! correlation between the soil liquid moisture content in the top
  !! 10 cm of soil (represented in CLASS by that of the first soil
  !! layer, \f$\theta_{1, 1}\f$) and the total surface albedo \f$\alpha_{g, T}\f$: for
  !! water contents less than 0.22 \f$m^3 m^{-3}\f$, \f$\alpha_{g, T}\f$ took on the value of
  !! \f$\alpha_{g, dry}\f$; for liquid water contents greater than 0.26 \f$m^3 m^{-3}\f$,
  !! \f$\alpha_{g, T}\f$ took on the value of \f$\alpha_{g, wet}\f$. For values of \f$\theta_{1, 1}\f$ between
  !! these two limits, a linear relationship is assumed:
  !!
  !! \f$[\alpha_{g, T} - \alpha_{g, dry} ] / [\theta_{l, 1} - 0.22] = [\alpha_{g, wet} - \alpha_{g, dry} ]/[0.26 - 0.22]\f$
  !!
  !! Thus, in groundAlbedo for each of the two wavelength ranges \f$\alpha_{g}\f$ is calculated as follows:
  !!
  !! \f$\alpha_{g, T} = \alpha_{g, dry}\f$            \f$\theta_{l, 1} \leq 0.22 \f$ \n
  !! \f$\alpha_{g, T} = \theta_{l, 1} [\alpha_{g, wet} - \alpha_{g, dry} ]/0.04 - 5.50 [\alpha_{g, wet} - \alpha_{g, dry} ] + \alpha_{g, dry} \f$
  !! \f$0.22 < \theta_{l, 1} < 0.26 \f$ \n
  !! \f$\alpha_{g, T} = \alpha{g, wet} \f$             \f$  0.26 \leq \theta_{l, 1}\f$
  !!
  !! Afterwards, a correction is applied to \f$\alpha_{g, VIS}\f$ and \f$\alpha_{g, NIR}\f$ in order to
  !! account for the possible presence of urban surfaces over the
  !! modelled area. Visible and near-IR albedos are assigned for local
  !! urban areas, \f$\alpha_{g, VIS}\f$ and \f$\alpha_{g, NIR}\f$, as part of the background data
  !! (see the section on “Data Requirements”). A weighted average over the bare soil area \f$X_g\f$ is
  !! calculated from the fractional snow-free urban area \f$X_u\f$ as:
  !!
  !! \f$\alpha_{g, VIS} = [X_u \alpha_{u, VIS} + (X_g-X_u) \alpha_{g, VIS}] / X_g \f$\n
  !! \f$\alpha_{g, NIR} = [X_u \alpha_{u, NIR} + (1.0-X_u) \alpha_{g, NIR}] / X_g \f$
  !!
  !! If the soil on the modelled area is not mineral, i.e. if the
  !! ISAND flag is less than zero, \f$\alpha_{g, VIS}\f$ and \f$\alpha_{g, NIR}\f$ are determined as
  !! follows:
  !!
  !! If ISAND = -2, indicating organic soil, \f$\alpha_{g, VIS}\f$ and \f$\alpha_{g, NIR}\f$ are
  !! assigned values of 0.05 and 0.30 respectively from the lookup
  !! tables in the block data subroutine soilPropertiesD, corresponding to
  !! average measured values reported in Comer et al. (2000) \cite Comer2000-mz.
  !!
  !! If ISAND = -3, indicating rock at the surface, \f$\alpha_{g, VIS}\f$ and
  !! \f$\alpha_{g, NIR}\f$ are assigned the dry ground values from soilProperties.
  !!
  !! If ISAND = -4, indicating continental ice sheet or glacier, \f$\alpha_{g, VIS}\f$
  !! and \f$\alpha_{g, NIR}\f$ are assigned values of 0.95 and 073 from soilPropertiesD,
  !! reflecting values reported for Antarctica (e.g. Sellers, 1974).
  !!
  !! The above calculations are all performed if the flag IALG is set
  !! to zero. If IALG is set to one, indicating that assigned ground
  !! albedos are to be used instead of calculated values, \f$\alpha_{g, VIS}\f$ and
  !! \f$\alpha_{g, NIR}\f$ are set to the assigned values AGVDAT and AGIDAT
  !! respectively.
  !!
  !! Lastly, the ground values of visible and near-IR albedo under the
  !! vegetation canopy are currently set equal to the open values
  !! (this approach is under review).
  !!
  IPTBAD = 0
  do I = IL1,IL2
    if (IALG == 0) then
      if (ISAND(I,1) >= 0) then
        FURB = FCMXU(I) * (1.0 - FSNOW(I))

        if (THLIQ(I,1) >= 0.26) then
          ALVSG(I) = ALGWV(I)
          ALIRG(I) = ALGWN(I)
        else if (THLIQ(I,1) <= 0.22) then
          ALVSG(I) = ALGDV(I)
          ALIRG(I) = ALGDN(I)
        else
          ALVSG(I) = THLIQ(I,1) * (ALGWV(I) - ALGDV(I)) / 0.04 + &
                     ALGDV(I) - 5.50 * (ALGWV(I) - ALGDV(I))
          ALIRG(I) = THLIQ(I,1) * (ALGWN(I) - ALGDN(I)) / 0.04 + &
                     ALGDN(I) - 5.50 * (ALGWN(I) - ALGDN(I))
        end if
        !
        if (FG(I) > 0.001) then
          ALVSG(I) = ((FG(I) - FURB) * ALVSG(I) + FURB * ALVSU(I)) / FG(I)
          ALIRG(I) = ((FG(I) - FURB) * ALIRG(I) + FURB * ALIRU(I)) / FG(I)
        end if
        if (ALVSG(I) > 1.0 .or. ALVSG(I) < 0.0) IPTBAD = I
        if (ALIRG(I) > 1.0 .or. ALIRG(I) < 0.0) IPTBAD = I

      else if (ISAND(I,1) == - 4) then
        ALVSG(I) = ALVSI
        ALIRG(I) = ALIRI
      else if (ISAND(I,1) == - 3) then
        ALVSG(I) = ALGDV(I)
        ALIRG(I) = ALGDN(I)
      else if (ISAND(I,1) == - 2) then
        ALVSG(I) = ALVSO
        ALIRG(I) = ALIRO
      end if
    else if (IALG == 1) then
      ALVSG(I) = AGVDAT(I)
      ALIRG(I) = AGIDAT(I)
    end if
    ALVSGC(I) = ALVSG(I)
    ALIRGC(I) = ALIRG(I)
  end do ! loop 100
  !
  if (IPTBAD /= 0) then
    write(6,6100) IPTBAD,JL,ALVSG(IPTBAD),ALIRG(IPTBAD)
6100 format('0AT (I,J) = (',I3,',',I3,'),ALVSG,ALIRG = ',2F10.5)
    call errorHandler('groundAlbedo', - 1)
  end if

  return
end subroutine groundAlbedo
