!> \file
!! Diagnoses snowpack visible and near-IR albedos given the
!! all-wave albedo at the current time step. Calculates snowpack
!! transmissivity for shortwave radiation.
!
subroutine snowAlbedoTransmiss (ALVSSN, ALIRSN, ALVSSC, ALIRSC, ALBSNO, & ! Formerly SNOWALBA
                                TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO, &
                                REFSN, BCSN, SNO, CSZ, ZSNOW, FSNOW, ASVDAT, ASIDAT, &
                                ALVSG, ALIRG, &
                                ILG, IG, IL1, IL2, JL, IALS, NBS, ISNOALB)
  !
  !     * JAN 27/16 - D.VERSEGHY. REFINE CALCULATIONS OF ALVSSN AND ALIRSN.
  !     * NOV 16/13 - J.COLE.     Final version for gcm17:
  !     *                         - Fixes to get the proper BC mixing ratio in
  !     *                           snow, which required passing in and using
  !     *                           the snow density RHON.
  !     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,
  !     *                         WHICH IS BASED ON 4-BAND SOLAR.
  !     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF
  !     *                         ALVSSN AND ALIRSN.
  !     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND
  !     *                         CANOPY-COVERED SNOW.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF
  !     *                         USER-SPECIFIED VALUES TO SNOW
  !     *                         ALBEDO.
  !     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         SPECIFY LOCATION OF ICE SHEETS
  !     *                         BY SOIL TEXTURE ARRAY RATHER
  !     *                         THAN BY SOIL COLOUR INDEX.
  !     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.
  !     *                         CALL ABORT CHANGED TO CALL errorHandler TO
  !     *                         ENABLE RUNNING ON PC'S.
  !     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -
  !     *                         DIVIDE PREVIOUS SUBROUTINE
  !     *                         "SNOALB" INTO "snowAlbedoTransmiss" AND
  !     *                         "snowAging" AND VECTORIZE.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO
  !     *                         VISIBLE AND NEAR-IR PORTIONS
  !     *                         CALCULATE TRANSMISSIVITY TO
  !     *                         SHORTWAVE RADIATION.
  !
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IG, IL1, IL2, JL, IALS
  integer, intent(in) :: NBS       !< Number of modelled shortwave radiation wavelength bands
  integer, intent(in) :: ISNOALB   !< Switch to model snow albedo in two or more wavelength bands
  integer             :: IPTBAD, I, IB
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out)   :: ALVSSC(ILG)  !< Visible albedo of snow on ground under vegetation canopy [ ]
  real, intent(out)   :: ALIRSC(ILG)  !< Near-IR albedo of snow on ground under vegetation canopy [ ]
  real, intent(out)   :: ALSNO (ILG,NBS) !< Albedo of snow in each modelled wavelength band  [  ]
  real, intent(out)   :: TRSNOWG(ILG,NBS)!< Transmissivity of snow in bare areas to shortwave radiation \f$[ ] (\tau_{s,g})\f$
  !
  real, intent(inout) :: TRSNOWC(ILG) !< Transmissivity of snow under vegetation to shortwave radiation
  !< \f$[ ] (\tau_{s, c})\f$
  real, intent(inout) :: ALVSSN(ILG)  !< Visible albedo of snow pack on bare ground \f$[ ]
  !< (alpha_{s, VIS})\f$
  real, intent(inout) :: ALIRSN(ILG)  !< Near-IR albedo of snow pack on bare ground \f$[ ]
  !< (alpha_{s, NIR})\f$


  !
  !     * INPUT ARRAYS.
  !
  real, intent(in)    :: ALVSG  (ILG)    !< Near-IR albedo f bare ground  [  ]
  real, intent(in)    :: ALIRG  (ILG)    !< Visible albedo of bare ground  [  ]
  real, intent(in)    :: FSDB(ILG,NBS)!< Direct solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
  real, intent(in)    :: FSFB(ILG,NBS)!< Diffuse solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
  real, intent(in)    :: ZSNOW (ILG)  !< Depth of snow \f$[m] (z_s)\f$
  real, intent(in)    :: FSNOW (ILG)  !< Fractional coverage of snow on grid cell [ ]
  real, intent(in)    :: ASVDAT(ILG)  !< Assigned value of visible albedo of snow pack –
  !< optional [ ]
  real, intent(in)    :: ASIDAT(ILG)  !< Assigned value of near-IR albedo of snow pack –
  !< optional [ ]
  real, intent(in)    :: REFSN (ILG)  !< Snow grain size  [m]
  real, intent(in)    :: BCSN  (ILG)  !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
  real, intent(in)    :: CSZ   (ILG)  !< Cosine of solar zenith angle  [  ]
  real, intent(in)    :: SNO   (ILG)  !< Mass of snow pack \f$[kg m^{-2}]\f$
  real, intent(in)    :: RHOSNO(ILG)  !< Density of snow pack  \f$[kg m^{-3}]\f$
  real, intent(inout) :: ALBSNO(ILG)  !< All-wave albedo of snow pack \f$[ ] (\alpha_{s, T})\f$

  !
  !     * LOCAL ARRAYS
  !
  real    :: SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF (ILG,NBS), &
             TRDIR(ILG,NBS), TRDIF (ILG,NBS), REFSNO(ILG), BCSNO(ILG)
  integer :: C_FLAG(ILG)
  !
  !     * CONSTANTS.
  !
  real    :: WDIRCT, WDIFF
  integer :: SUM_C_FLAG
  !------------------------------------------------------------------

  IPTBAD = 0
  do I = IL1,IL2 ! loop 100
    if (ALBSNO(I) < 0.50 .and. ALBSNO(I) > 0.499) ALBSNO(I) = 0.50
    if (FSNOW(I) > 0.0 .and. IALS == 0) then
      if (ALBSNO(I) > 0.70) then
        ALVSSN(I) = 0.7857 * ALBSNO(I) + 0.2900
        ALIRSN(I) = 1.2142 * ALBSNO(I) - 0.2900
      else
        ALVSSN(I) = 0.9706 * ALBSNO(I) + 0.1347
        ALIRSN(I) = 1.0294 * ALBSNO(I) - 0.1347
      end if
      if (ALVSSN(I) > 0.999 .or. ALVSSN(I) < 0.001) IPTBAD = I
      if (ALIRSN(I) > 0.999 .or. ALIRSN(I) < 0.001) IPTBAD = I
    else if (FSNOW(I) > 0.0 .and. IALS == 1) then
      ALVSSN(I) = ASVDAT(I)
      ALIRSN(I) = ASIDAT(I)
    end if
    ALVSSC(I) = ALVSSN(I)
    ALIRSC(I) = ALIRSN(I)
    TRSNOWC(I) = EXP( - 25.0 * ZSNOW(I))
  end do ! loop 100
  !
  if (IPTBAD /= 0) then
    write(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)
6100 format('0AT (I,J) = (',I3,',',I3,'),ALVSSN,ALIRSN = ',2F10.5)
    call errorHandler('snowAlbedoTransmiss', - 1)
  end if
  !
  if (ISNOALB == 0) then
    do I = IL1,IL2
      ALSNO(I,1) = ALVSSN(I)
      ALSNO(I,2) = ALIRSN(I)
      ALSNO(I,3) = ALIRSN(I)
      ALSNO(I,4) = ALIRSN(I)

      TRSNOWG(I,1:NBS) = TRSNOWC(I)
    end do ! I
  else if (ISNOALB == 1) then
    do IB = 1,NBS
      do I = IL1,IL2
        if (IB == 1) then
          SALBG(I,IB) = ALVSG(I)
          ALSNO(I,IB) = ALVSSN(I)
        else
          SALBG(I,IB) = ALIRG(I)
          ALSNO(I,IB) = ALIRSN(I)
        end if
      end do ! I
    end do ! IB
    SUM_C_FLAG = 0
    do I = IL1,IL2
      if (ZSNOW(I) > 0.0) then
        C_FLAG(I) = 1
      else
        C_FLAG(I) = 0
      end if
      SUM_C_FLAG = SUM_C_FLAG + C_FLAG(I)
    end do ! I

    if (IALS == 0) then
      if (SUM_C_FLAG > 0) then
        !>
        !! Convert the units of the snow grain size and BC mixing ratio
        !! Snow grain size from meters to microns and BC from \f$kg BC/m^3\f$ to ng BC/kg SNOW
        !!
        do I = IL1,IL2
          if (C_FLAG(I) == 1) then
            REFSNO(I) = REFSN(I) * 1.0E6
            BCSNO(I)  = (BCSN(I) / RHOSNO(I)) * 1.0E12
          end if
        end do ! I

        call SNOW_ALBVAL(ALDIF, & ! OUTPUT
                         ALDIR, &
                         CSZ, & ! INPUT
                         SALBG, &
                         BCSNO, &
                         REFSNO, &
                         SNO, &
                         C_FLAG, &
                         IL1, &
                         IL2, &
                         ILG, &
                         NBS)

        call SNOW_TRANVAL(TRDIF, & ! OUTPUT
                          TRDIR, &
                          CSZ, &  ! INPUT
                          SALBG, &
                          BCSNO, &
                          REFSNO, &
                          SNO, &
                          C_FLAG, &
                          IL1, &
                          IL2, &
                          ILG, &
                          NBS)

        do IB = 1,NBS
          do I = IL1,IL2
            if (C_FLAG(I) == 1) then
              WDIRCT = FSDB(I,IB) &
                       / (FSDB(I,IB) + FSFB(I,IB) + 1.E-10)
              WDIFF  = 1.0 - WDIRCT
              ALSNO(I,IB) = ALDIF (I,IB) * WDIFF &
                            + ALDIR(I,IB) * WDIRCT
              TRSNOWG(I,IB) = TRDIF (I,IB) * WDIFF &
                              + TRDIR(I,IB) * WDIRCT
            end if ! C_FLAG
          end do ! I
        end do ! IB
      else ! SUM_C_FLAG == 0
        do I = IL1,IL2
          ALSNO(I,1)     = ALVSSN(I)
          ALSNO(I,2:NBS) = ALIRSN(I)
          TRSNOWG(I,1:NBS) = TRSNOWC(I)
        end do ! I
      end if ! SUM_C_FLAG
    else if (IALS == 1) then
      do I = IL1,IL2
        ALSNO(I,1) = ASVDAT(I)
        ALSNO(I,2) = ASIDAT(I)
        ALSNO(I,3) = ASIDAT(I)
        ALSNO(I,4) = ASIDAT(I)

        TRSNOWG(I,1:NBS) = TRSNOWC(I)
      end do ! I
    end if ! IALS
  end if ! ISNOALB
  return
end subroutine snowAlbedoTransmiss
!> \file
!! @author D. Verseghy, J. Cole, M. Lazare
!!
!! In subroutine snowAging, called at the end of waterBudgetDriver, the change of
!! total snow albedo over the current time step is calculated using
!! an empirical exponential decay function, which has different
!! coefficients depending on whether the snow is dry or melting. In
!! this subroutine, if the ISNOALB switch is set to 0, the visible and
!! near-IR components of the snow albedo are diagnosed from the total albedo.
!! According to the literature (Aguado, 1985 \cite Aguado1985-fv ; Robinson and Kukla, 1984; Dirmhirn and
!! Eaton, 1975 \cite Dirmhirn1975-vx), the following represent typical snow albedos for
!! fresh snow, old dry snow and melting snow:
!!
!! \f[
!! \begin{array} { | l | c | c | c | }
!! \hline
!!              & \text{Total albedo}  & \text{Visible albedo} & \text{Near-IR albedo} \\ \hline
!! \text{Fresh snow}   &      0.84     &      0.95      &      0.73      \\ \hline
!! \text{Old dry snow} &      0.70     &      0.84      &      0.56      \\ \hline
!! \text{Melting snow} &      0.50     &      0.62      &      0.38      \\ \hline
!! \end{array}
!! \f]
!!
!! The same decay function is assumed to apply to all three albedo
!! ranges, so the relative location of the visible and near-IR
!! albedos, \f$\alpha_{s, VIS}\f$ and \f$\alpha_{s, NIR}\f$, on the decay curve will be analogous
!! to that of the total albedo, \f$\alpha_{s, T}\f$. Thus, for dry snow:
!!
!! \f$[\alpha_{s, VIS} - 0.84]/[0.95-0.84] = [\alpha_{s, T} - 0.70]/[0.84-0.70]\f$
!! \f$[\alpha_{s, NIR} - 0.56]/[0.73-0.56] = [\alpha_{s, T} - 0.70]/[0.84-0.70]\f$
!!
!! or, simplifying:
!!
!! \f$\alpha_{s, VIS} = 0.7857 \alpha_{s, T} + 0.2900\f$
!! \f$\alpha_{s, NIR} = 1.2142 \alpha_{s, T} - 0.2900\f$
!!
!! For melting snow:
!!
!! [\f$\alpha_{s, VIS} - 0.62]/[0.95-0.62] = [\alpha_{s, T} - 0.50]/[0.84-0.50]\f$
!! [\f$\alpha_{s, NIR} - 0.38]/[0.73-0.38] = [\alpha_{s, T} - 0.50]/[0.84-0.50]\f$
!!
!! or, simplifying:
!!
!! \f$\alpha_{s, VIS} = 0.9706 \alpha_{s, T} + 0.1347\f$
!! \f$\alpha_{s, NIR} = 1.0294 \alpha_{s, T} - 0.1347\f$
!!
!! The above calculations are performed if the flag IALS is set to
!! zero. If IALS is set to one, indicating that assigned snow
!! albedos are to be used instead of calculated values, \f$\alpha_{s, VIS}\f$ and
!! \f$\alpha_{s, NIR}\f$ are set to the assigned values ASVDAT and ASIDAT
!! respectively. The sub-canopy values of visible and near-IR albedo
!! are currently set equal to the open snowpack values (this is
!! expected to change if a canopy litterfall parametrization is
!! developed).
!!
!! The transmissivity of snow under vegetation \f$\tau_{s, c}\f$ is
!! then calculated from the snow depth
!! ZSNOW using Beer’s law, with an empirical extinction coefficient
!! of \f$25.0 m^{-1}\f$ derived from the literature (Grenfell and Maykut,
!! 1977 \cite Grenfell1977-pi ; Thomas, 1963):
!!
!! \f$\tau_{s, c} = exp[-25.0 z_s]\f$
!!
!! If the ISNOALB switch is set to zero, the value of ALSNO in the first
!! wavelength band is set to the previously calculated value of \f$\alpha_{s, VIS}\f$
!! and the values for the remaining bands are set to the previously calculated value
!! of \f$\alpha_{s, NIR}\f$; the value of \f$\tau_{s, g}\f$  is set to \f$\tau_{s, c}\f$.
!! If the ISNOALB switch is set to 1, a new parameterization for the snow albedo and
!! transmissivity in four shortwave radiation bands (one visible and three near-IR)
!! is used, according to Cole et al. (2017).  This parameterization incorporates
!! the effects of snow grain size and black carbon content, and makes use of lookup
!! tables contained in the CCCma subroutines SNOW_ALBVAL and SNOW_TRANVAL.
