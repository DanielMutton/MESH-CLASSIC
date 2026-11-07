
!> \file
!! Assigns thermal and hydraulic properties to soil layers based on sand/clay content, or soil type.
!! Also calculates permeable thickness of soil layers, and wet and dry surface albedo for mineral soils.
!! @author D. Verseghy, M. Lazare, V. Fortin, Y. Wu, J. Melton
!!
subroutine soilPropertiesGAT (THPOR, THLRET, THLMIN, BI, PSISAT, GRKSAT, & ! Formerly CLASSB
                              THLRAT, HCPS, TCS, THFC, THLW, PSIWLT, &
                              DELZW, ZBOTW, ALGWV, ALGWN, ALGDV, ALGDN, &
                              SAND, CLAY, ORGM, SOCI, DELZ, ZBOT, SDEPTH, &
                              ISAND, IGDR, NML, IL1, IL2, IG, ipeatland)
  !
  !     * APR 4/17  - J. Melton   TCFINE was here in place of TCCLAY,
  !                               Change to TCCLAY for consistency with rest of model.
  !     * DEC 16/16 - D.VERSEGHY. REMOVE OPTION FOR USING OLD SOIL
  !     *                         WET AND DRY ALBEDOS DETERMINED BY
  !     *                         SOIL TEXTURE ("IGRALB" SWITCH)
  !     * SEP 28/16 - J. Melton.  Finish bringing in Yuanqiao Wu's peatland work.
  !     * JUN 24/15 - J. MELTON.  PASS IN IGRALB SO THAT WE CAN SKIP
  !                               USING SOCI IF IGRALB IS 0.
  !     * FEB 09/15 - D.VERSEGHY. New version for gcm18 and class 3.6:
  !     *                         - PSIWLT is now a constant value
  !     *                           of 150, the same as used with CTEM.
  !     *                         - A new field THLW (wilting point)
  !     *                           is defined here, the same as used
  !     *                           with CTEM. It is passed out and
  !     *                           used in energyBudgetDriver.
  !     *                         - Wet and dry albedoes for EACH of
  !     *                           visible and near-ir are defined
  !     *                           using a lookup-table based on the
  !     *                           new soil colour index field SOCI
  !     *                           which is now passed in. These
  !     *                           are known as {ALGWV, ALGWN, ALGDV, ALGDN}
  !     *                           and replace {ALGWET, ALGDRY}.
  !     * JAN 15/15 - D.VERSEGHY. CHANGE PSIWLT FOR MINERAL SOILS
  !     *                         TO A CONSTANT VALUE OF 150 M.
  !     *                         AND ADD NEW VARIABLE THLW.
  !     * AUG 25/14 - M.LAZARE.   PASS IN NEW WET AND DRY SOIL
  !     *                         BRIGHTNESS FIELDS FROM CLM.
  !     * NOV 16/13 - M.LAZARE.   FINAL VERSION FOR GCM17:
  !     *                         - REVERT BACK TO CLASS2.7
  !     *                           SPECIFICATION FOR "ALGWET".
  !     * NOV 11/11 - M.LAZARE.   - IMPLEMENT CTEM CHOICE OF
  !     *                           ALGDRY DETERMINED BY ADDED
  !     *                           PASSED SWITCH "ICTEMMOD".
  !     * OCT 18/11 - M.LAZARE.   - REMOVE UNUSED "IORG".
  !     *                         - CHANGE "THSAND", "THORG"
  !     *                           AND "THFINE" FROM ARRAYS
  !     *                           (INTERNAL ONLY) TO SCALAR.
  !     *                         - IGDR NOW PASSED OUT TO BE
  !     *                           USED IN waterFlowInfiltrate/waterFlowNonInfiltrate/waterBaseflow.
  !     *                         - PASS IN IL1 AND IL2 TO
  !     *                           DEFINE LOOPS.
  !     * OCT 08/11 - M.LAZARE.   ALGDRY CHANGED BACK TO FORMULA
  !     *                         USED IN GCM15I (.0056->.0046).
  !     * SEP 27/11 - D.VERSEGHY. CONSTRAIN DELZW TO BE >= 5 CMS
  !     *                         TO AVOID PROBLEMATIC UNDERSHOOTS.
  !     * AUG 25/11 - D.VERSEGHY. USE THFC FORMULATION FOR BOTTOM
  !     *                         LAYER AT BOTTOM OF SOIL PERMEABLE
  !     *                         DEPTH.
  !     * DEC 23/09 - V.FORTIN.   REVISE CALCULATION OF THFC FOR
  !     *                         BOTTOM LAYER IN MINERAL SOILS
  !     *                         ACCORDING TO SOULIS ET AL. (2009).
  !     * JAN 06/09 - D.VERSEGHY. REVERSE ORDER OF 200 AND 300 LOOPS.
  !     * DEC 11/07 - D.VERSEGHY. CHANGE CALCULATION OF TCS FROM
  !     *                         GEOMETRIC MEAN TO LINEAR MEAN.
  !     * FEB 07/07 - D.VERSEGHY. SET THFC TO THLRET FOR ORGANIC SOILS
  !     *                         STREAMLINE SOME CALCULATIONS.
  !     * SEP 15/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 300 LOOP.
  !     * APR 06/05 - D.VERSEGHY. MOVE CALCULATION OF GRKTLD
  !     *                         INTO waterFlowInfiltrate; REVISED CALCULATION
  !     *                         OF ALGDRY (WITH M.LAZARE).
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * SEP 04/03 - D.VERSEGHY. PERMEABLE THICKNESS OF SOIL
  !     *                         LAYERS CONSTRAINED TO >= 1 MM
  !     *                         PROTECT SENSITIVE CALCULATIONS
  !     *                         AGAINST ROUNDOFF ERROR.
  !     * JUN 28/02 - D.VERSEGHY. ASSIGN SOIL HYDROLOGICAL AND
  !     *                         THERMAL PROPERTIES BASED ON
  !     *                         SAND, CLAY AND ORGANIC MATTER
  !     *                         CONTENT.
  !
  use classicParams, only : thpmoss, thrmoss, thmmoss, bmoss, psismoss, &
                             grksmoss, hcpmoss, THPORG, THRORG, THMORG, BORG, &
                             PSISORG, GRKSORG, TCICE, TCSAND, TCCLAY, TCOM, &
                             RHOSOL, RHOOM, HCPICE, HCPOM, HCPSND, HCPCLY, &
                             ALWV, ALWN, ALDV, ALDN


  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: NML, IL1, IL2, IG
  integer :: I, J, k
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: THPOR (NML,IG) !< Pore volume \f$[m^3 m^{-3} ] ( \theta_p)\f$
  real, intent(out) :: THLRET(NML,IG) !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ] (\theta_{ret} )\f$
  real, intent(out) :: THLMIN(NML,IG) !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ] (\theta_{min} )\f$
  real, intent(out) :: BI    (NML,IG) !< Clapp and Hornberger empirical parameter [ ] (b)
  real, intent(out) :: PSISAT(NML,IG) !< Soil moisture suction at saturation [m] \f$(\Psi_{sat} )\f$
  real, intent(out) :: GRKSAT(NML,IG) !< Hydraulic conductivity of soil at saturation \f$[m s^{-1} ] (K_{sat} )\f$
  real, intent(out) :: THLRAT(NML,IG) !< Fractional saturation of soil at half the saturated hydraulic conductivity [ ] \f$(f_{inf} )\f$
  real, intent(out) :: HCPS  (NML,IG) !< Volumetric heat capacity of soil matter \f$[J m^{-3} K^{-1} ] (C_g)\f$
  real, intent(out) :: TCS   (NML,IG) !< Thermal conductivity of soil \f$[W m^{-1} K^{-1} ] (\tau_g)\f$
  real, intent(out) :: THFC  (NML,IG) !< Field capacity \f$[m^3 m^{-3} ] (\theta_{fc} )\f$
  real, intent(out) :: THLW  (NML,IG) !< Soil water content at wilting point, \f$[m^3 m^{-3} ] (\theta_{wilt})\f$
  real, intent(out) :: PSIWLT(NML,IG) !< Soil moisture suction at wilting point [m] \f$(\Psi_{wilt} )\f$
  real, intent(out) :: DELZW (NML,IG) !< Thickness of permeable part of soil layer [m]
  real, intent(out) :: ZBOTW (NML,IG) !< Depth of bottom of permeable part of soil layer [m]
  real, intent(out) :: ALGWV (NML)    !< Visible albedo of wet soil for modelled area  [  ]
  real, intent(out) :: ALGWN (NML)    !< Near-IR albedo of wet soil for modelled area  [  ]
  real, intent(out) :: ALGDV (NML)    !< Visible albedo of dry soil for modelled area  [  ]
  real, intent(out) :: ALGDN (NML)    !< Near-IR albedo of dry soil for modelled area  [  ]
  !
  integer, intent(inout) :: ISAND (NML,IG) !< Sand content flag
  integer, intent(inout) :: IGDR  (NML) !< Index of soil layer in which bedrock is encountered
  !
  !     * INPUT ARRAYS.
  !
  integer, intent(in) :: ipeatland(nml)  !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
  real, intent(in) :: SAND  (NML,IG) !< Percent sand content of soil layer [percent] \f$(X_{sand} )\f$
  real, intent(in) :: CLAY  (NML,IG) !< Percent clay content of soil layer [percent] \f$(X_{clay} )\f$
  real, intent(in) :: ORGM  (NML,IG) !< Percent organic matter content of soil layer [percent]
  real, intent(in) :: DELZ  (IG)       !< Thickness of soil layer [m]
  real, intent(in) :: ZBOT  (IG)       !< Depth of bottom of soil layer [m]
  real, intent(in) :: SDEPTH(NML)    !< Permeable depth of soil column (depth to bedrock) [m] \f$(z_b)\f$
  real, intent(in) :: SOCI  (NML)    !< Soil colour index
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: VSAND, VORG, VFINE, VTOT, AEXP, ABC, THSAND, THFINE, THORG
  !
  !---------------------------------------------------------------------
  !
  !>
  !! In the first section of code, two integer :: flags are evaluated: IGDR, the index of the soil layer in which the
  !! bottom of the soil permeable depth, \f$z_b\f$ , occurs; and ISAND, the value of the SAND variable for each soil
  !! layer converted to an integer. ISAND is used throughout the CLASS code as a flag to determine which
  !! branches of code to execute.
  !!

!MESH  do M = 1,IM ! loop 50
    do I = IL1,IL2
      IGDR(I) = 1
    end do
!MESH  end do ! loop 50
  !
  do J = 1,IG ! loop 100
!MESH    do M = 1,IM
      do I = IL1,IL2
        ISAND (I,J) = NINT(SAND(I,J))
        if (ISAND(I,J) > - 3) IGDR(I) = J
      end do
!MESH    end do
  end do ! loop 100
  !
  !>
  !! In loop 200, calculations are done to determine the permeable thickness DELZW of each soil layer. If
  !! the land cover is an ice sheet (indicated by an ISAND value of -4 in the top layer), the value of DELZW
  !! in each layer is set to DELZ, the standard thickness of the corresponding thermal layer, and the soil flag is
  !! set to -4. If the layer consists of rock, indicated by an ISAND value of -3, DELZW is set to 0. For soil
  !! layers where \f$z_b\f$ occurs below the bottom of the layer, DELZW is set to DELZ; for layers where \f$z_b\f$ occurs
  !! near the top of the layer, DELZW is set to 0 and ISAND is set to -3. For the soil layer containing \f$z_b\f$ ,
  !! DELZW is set to the distance from the top of the layer to \f$z_b\f$ , and is further constrained to be \f$\geq\f$ 5 cm, to
  !! avoid overshoots in the water movement calculations. For all layers, the distance to the bottom of its
  !! respective permeable depth, ZBOTW, is set to the depth of the top of the layer plus DELZW. At the end of the loop,
  !! if the soil is a mineral one, the wet and dry visible and near-IR soil albedo values are assigned using the lookup
  !! tables ALWV, ALWN, ALDV and ALDN found in the DATA section of the subroutine, on the basis of a global gridded
  !! soil “colour” (i.e. reflectivity) index SOCI, documented in Lawrence and Chase (2007) \cite Lawrence2007-bc and
  !! Oleson et al. (2010) \cite Oleson2010-c88.
  !!

!MESH  do M = 1,IM ! loop 200
    do I = IL1,IL2
      do J = 1,IG ! loop 150
        if (ISAND(I,1) == - 4) then
          DELZW(I,J) = DELZ(J)
          ISAND(I,J) = - 4
        else if (ISAND(I,J) == - 3) then
          DELZW(I,J) = 0.0
        else if (SDEPTH(I) >= ZBOT(J)) then
          DELZW(I,J) = DELZ(J)
        else if (SDEPTH(I) < (ZBOT(J) - DELZ(J) + 0.025)) then
          DELZW(I,J) = 0.0
          ISAND(I,J) = - 3
        else
          DELZW(I,J) = MAX(0.05,(SDEPTH(I) - (ZBOT(J) - DELZ(J))))
        end if
        ZBOTW(I,J) = MAX(0.0,ZBOT(J) - DELZ(J)) + DELZW(I,J)
      end do ! loop 150
      if (SAND(I,1) >= - 3.5) then
        ALGWV(I) = ALWV(NINT(SOCI(I)))
        ALGWN(I) = ALWN(NINT(SOCI(I)))
        ALGDV(I) = ALDV(NINT(SOCI(I)))
        ALGDN(I) = ALDN(NINT(SOCI(I)))
      else
        ALGWV(I) = 0.0
        ALGWN(I) = 0.0
        ALGDV(I) = 0.0
        ALGDN(I) = 0.0
      end if
    end do
!MESH  end do ! loop 200
  !
  !>
  !! In loop 300, various thermal and hydraulic soil properties are assigned to each of the soil layers,
  !! depending on soil type. Values of ISAND greater than or equal to zero indicate mineral soil. The pore volume \f$\theta_p\f$ ,
  !! the saturated hydraulic conductivity \f$K_{sat}\f$ , and the soil moisture suction at saturation \f$\Psi\f$ sat are calculated from
  !! the percentage sand content \f$X_{sand}\f$ , and the hydraulic parameter b is calculated from the percentage clay
  !! content \f$X_{clay}\f$ , based on empirical relationships given in Cosby et al. (1984) \cite Cosby1984-jc
  !!
  !! \f$\theta_p = (-0.126 X_{sand} +48.9)/100.0\f$
  !!
  !! \f$b = 0.159 X_{clay} + 2.91\f$
  !!
  !! \f$\Psi_{sat} = 0.01 exp(-0.0302 X_{sand} + 4.33)\f$
  !!
  !! \f$K_{sat} = 7.0556 x 10 -6 exp(0.0352 X_{sand} - 2.035)\f$
  !!
  !! The fractional saturation of the soil at half the saturated hydraulic conductivity, \f$f_{inf}\f$ , is calculated by
  !! inverting the Clapp and Hornberger (1978) \cite Clapp1978-898 expression relating hydraulic conductivity \f$K\f$ to liquid water
  !! content of the soil \f$\theta_l\f$ :
  !! \f$K = K_{sat} (\theta_l / \theta_p) (2b + 3)\f$
  !!
  !! Thus,
  !! \f$f_{inf} = 0.51/(2b+3)\f$
  !!
  !! The residual soil liquid water content remaining after evaporation or freezing, \f$\theta_{min}\f$ , and the liquid water
  !! retention capacity, \f$\theta_{ret}\f$ , are both set for mineral soils to a textbook value of 0.04.
  !! The volumetric sand, silt + clay and organic matter components of the soil matrix are derived by
  !! converting the percent values to volume fractions. The overall volumetric heat capacity of the soil
  !! material, \f$C_g\f$ , is then calculated as a weighted average:
  !! \f$C_g = \Sigma (C_{sand} \theta_{sand} + C_{fine} \theta_{fine} + C_{org} \theta_{org} )/(1 - \theta_p)\f$
  !! where the subscript "fine" refers to the silt and clay particles taken together. The thermal conductivity of
  !! the soil material, \f$\tau_g\f$ , is likewise calculated as a weighted average over the thermal conductivities of the
  !! components:
  !! \f$\tau_g = \Sigma (\tau_{sand} \theta_{sand} + \tau_{fine} \theta_{fine} + \tau_{org} \theta_{org} )/(1 - \theta_p)\f$
  !!
  !! The field capacity \f$\theta_{fc}\f$ , that is, the liquid water content of the soil at which gravitational drainage effectively
  !! ceases, is calculated by setting the expression for K above to a value of 0.1 mm \f$d^{-1}\f$ , and solving for the
  !! liquid water content:
  !! \f$\theta_{fc} = \theta_p (1.157 x 10^{-9} /K_{sat} )^{1/(2b + 3)}\f$
  !!
  !! The only exception is the field capacity of the lowest permeable layer in mineral soils (layer IGDR), which
  !! is determined using an expression developed by Soulis et al. (2010), which takes into account the
  !! permeable depth of the whole overlying soil:
  !! \f$\theta_{fc} = \theta_p /(b-1) \bullet (\Psi_{sat} b/ z_b)^{1/b} \bullet [(3b+2)^{(b-1)/b} - (2b+2)^{(b-1)/b} ]\f$
  !!
  !! The soil moisture suction \f$\Psi_{wilt}\f$ at the wilting point (the liquid water content at which plant roots can no
  !! longer draw water from the soil) is assigned a textbook value of 150.0 m. The liquid water content at the wilting point,
  !! \f$ \theta_{wilt} \f$, is determined from \f$\Psi_{wilt}\f$ by inverting the calculated from the saturated soil moisture
  !! suction using the Clapp and Hornberger (1978) \cite Clapp1978-898 expression for soil moisture suction as a function of
  !! liquid water content: \f$ \theta_{wilt} = \theta_p (\theta_{wilt} / \Psi_{sat} )^{-1/b} \f$
  !!
  !! Organic soils are flagged with an ISAND value of -2. For these soils, the variables \f$\theta_p\f$ , b, \f$K_{sat}\f$ , \f$\Psi_{sat}\f$ , \f$\theta_{min}\f$ , and
  !! \f$\theta_{ret}\f$ are assigned values based on the peat texture (fibric, hemic or sapric). These values are obtained from
  !! the arrays in common block CLASS5 (see the soilProperties documentation above). The current default is
  !! to assume the first soil layer as fibric, the second as hemic, and any lower layers as sapric. The volumetric
  !! heat capacity and thermal conductivity are set to textbook values for organic matter; the field capacity is
  !! set equal to the retention capacity; \f$f_{inf}\f$ is obtained as above; and \f$\Psi_{wilt}\f$ is assume to apply at a liquid water
  !! content of \f$\theta_{min}\f$ .
  !!
  !! In the cases of rock soils and ice sheets (respectively flagged with ISAND values of -3 and -4), all of the
  !! above variables are set to zero except for the volumetric heat capacity and the thermal conductivity, which
  !! are assigned values representative of rock or ice.
  !!

  do J = 1,IG
!MESH    do M = 1,IM
      do I = IL1,IL2
        if (ISAND(I,J) == - 4) then
          THPOR (I,J) = 0.0
          THLRET(I,J) = 0.0
          THLMIN(I,J) = 0.0
          BI    (I,J) = 0.0
          PSISAT(I,J) = 0.0
          GRKSAT(I,J) = 0.0
          THLRAT(I,J) = 0.0
          TCS(I,J) = TCICE
          HCPS(I,J) = HCPICE
          THFC(I,J) = 0.0
          THLW(I,J) = 0.0
          PSIWLT(I,J) = 0.0
        else if (ISAND(I,J) == - 3) then
          THPOR (I,J) = 0.0
          THLRET(I,J) = 0.0
          THLMIN(I,J) = 0.0
          BI    (I,J) = 0.0
          PSISAT(I,J) = 0.0
          GRKSAT(I,J) = 0.0
          THLRAT(I,J) = 0.0
          HCPS(I,J) = HCPSND
          TCS(I,J) = TCSAND
          THFC(I,J) = 0.0
          THLW(I,J) = 0.0
          PSIWLT(I,J) = 0.0
        else if (ISAND(I,J) == - 2) then

          if (j == 1) then
            k = 1
          else if (j == 2 .or. j == 3) then
            k = 2
          else
            k = 3
          end if
          THPOR (I,J) = THPORG(k)
          THLRET(I,J) = THRORG(k)
          THLMIN(I,J) = THMORG(k)
          BI    (I,J) = BORG(k)
          PSISAT(I,J) = PSISORG(k)
          GRKSAT(I,J) = GRKSORG(k)
          THLRAT(I,J) = 0.5 ** (1.0 / (2.0 * BI(I,J) + 3.0))
          HCPS(I,J) = HCPOM
          TCS(I,J) = TCOM

          ! If tile is a peatland, then set up as such.
          ! this setup assumes that the soil has 10 cm
          ! layers in the top 1 meter.
          ! FLAG: We should consider here that we don't
          ! necessarily want the entire soil column as
          ! peat so should also include the peat depth ! JM Oct 2016.
          if (ipeatland(i) > 0) then ! Peatland flag, 1 = bog, 2 = fen
            if (j == 1) then ! First layer is moss
              thpor(i,j)  = thpmoss
              thlret(i,j) = thrmoss
              thlmin(i,j) = thmmoss
              bi(i,j)     = bmoss
              psisat(i,j) = psismoss
              grksat(i,j) = grksmoss
              hcps(i,j) = hcpmoss
              tcs(i,j) = tcom
            else if (j == 2) then ! Next treated as soil and is fibric peat
              thpor(i,j)  = thporg(1)
              thlret(i,j) = throrg(1)
              thlmin(i,j) = thmorg(1)
              bi(i,j)     = borg(1)
              psisat(i,j) = psisorg(1)
              grksat(i,j) = grksorg(1)
            else if (j >= 3 .and. j <= 5) then ! These layers are considered hemic peat
              thpor(i,j)  = thporg(2)
              thlret(i,j) = throrg(2)
              thlmin(i,j) = thmorg(2)
              bi(i,j)     = borg(2)
              psisat(i,j) = psisorg(2)
              grksat(i,j) = grksorg(2)
            else ! Remainder are sapric peat
              thpor(i,j)  = thporg(3)
              thlret(i,j) = throrg(3)
              thlmin(i,j) = thmorg(3)
              bi(i,j)     = borg(3)
              psisat(i,j) = psisorg(3)
              grksat(i,j) = grksorg(3)
            end if
            thlrat(i,j) = 0.5 ** (1.0 / (2.0 * bi(i,j) + 3.0))
          end if
          THFC(I,J) = THLRET(I,J)
          THLW(I,J) = THLMIN(I,J)
          PSIWLT(I,J) = PSISAT(I,J) * (THLMIN(I,J) / &
                          THPOR(I,J)) ** ( - BI(I,J))

        else if (SAND(I,J) >= 0.) then
          THPOR (I,J) = ( - 0.126 * SAND(I,J) + 48.9) / 100.0
          THLRET(I,J) = 0.04
          THLMIN(I,J) = 0.04
          BI    (I,J) = 0.159 * CLAY(I,J) + 2.91
          PSISAT(I,J) = 0.01 * EXP( - 0.0302 * SAND(I,J) + 4.33)
          GRKSAT(I,J) = 7.0556E-6 * (EXP(0.0352 * SAND(I,J) - 2.035))
          THLRAT(I,J) = 0.5 ** (1.0 / (2.0 * BI(I,J) + 3.0))
          VSAND = SAND(I,J) / (RHOSOL * 100.0)
          VORG = ORGM(I,J) / (RHOOM * 100.0)
          VFINE = (100.0 - SAND(I,J) - ORGM(I,J)) / (RHOSOL * 100.0)
          VTOT = VSAND + VFINE + VORG
          THSAND = (1.0 - THPOR(I,J)) * VSAND / VTOT
          THORG = (1.0 - THPOR(I,J)) * VORG / VTOT
          THFINE = 1.0 - THPOR(I,J) - THSAND - THORG
          HCPS(I,J) = (HCPSND * THSAND + HCPCLY * THFINE + &
                        HCPOM * THORG) / (1.0 - THPOR(I,J))
          TCS(I,J) = (TCSAND * THSAND + TCOM * THORG + &
                       TCCLAY * THFINE) / (1.0 - THPOR(I,J))
          if (J /= IGDR(I)) then
            THFC(I,J) = THPOR(I,J) * (1.157E-9 / GRKSAT(I,J)) ** &
                          (1.0 / (2.0 * BI(I,J) + 3.0))
          else
            AEXP = (BI(I,J) - 1.0) / BI(I,J)
            ABC = (3.0 * BI(I,J) + 2.0) ** AEXP - &
                  (2.0 * BI(I,J) + 2.0) ** AEXP
            THFC(I,J) = (ABC * THPOR(I,J) / (BI(I,J) - 1.0)) * &
                          (PSISAT(I,J) * BI(I,J) / SDEPTH(I)) ** &
                          (1.0 / BI(I,J))
          end if
          PSIWLT(I,J) = 150.0
          THLW(I,J) = THPOR(I,J) * (PSIWLT(I,J) / PSISAT(I,J)) ** &
                        ( - 1.0 / BI(I,J))
        end if
      end do
!MESH    end do
  end do ! loop 300
  !
  return
end subroutine soilPropertiesGAT
