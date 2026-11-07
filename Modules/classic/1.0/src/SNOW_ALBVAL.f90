!> \file
!> Computes the direct and diffuse snow albedo using lookup table and current snow conditions.
!! @author J. Cole
!!
subroutine SNOW_ALBVAL (albdif, albdir, & ! OUTPUT
                        smu, salb, bc_conc, snow_reff, & ! INPUT
                        swe, c_ind, il1, il2, ilg, nbnd)
  !
  !     * JAN 24, 2013 - J. COLE
  !                    - COMPUTES THE DIRECT AND DIFFUSE SNOW ALBEDO
  !                      USING LOOKUP TABLE AND CURRENT SNOW CONDITIONS.
  !
  implicit none
  !
  ! INPUT
  !
  real, intent(in)    :: smu(ILG) !< COSINE OF THE SOLAR ZENITH ANGLE [UNITLESS]
  real, intent(in)    :: bc_conc(ILG) !< CONCENTRATION OF BLACK CARBON IN THE SNOW PACK [NG (BC)/KG (SNOW)]
  real, intent(in)    :: snow_reff(ILG) !< EFFECTIVE RADIUS OF THE SNOW GRAIN [MICRONS]
  real, intent(in)    :: swe(ILG) !< SNOW WATER EQUIVALENT (SNOWPACK DENSITY*SNOW PACK DEPTH) [KG/M^2]
  real, intent(in)    :: salb(ILG,NBND) !< ALBEDO OF THE UNDERLYING SURFACE [UNITLESS]
  integer, intent(in) :: c_ind(ILG) !< INDICATOR THAT A CALCULATION SHOULD BE PERFORMED FOR THIS POINT 1-YES, 0-NO
  integer, intent(in) :: il1 !< STARTING POINT FOR ALBEDO CALCULATIONS
  integer, intent(in) :: il2 !< ENDING POINT FOR ALBEDO CALCULATIONS
  integer, intent(in) :: ilg !< NUMBER OF POINTS FOR WHICH TO COMPUTE ALBEDOS
  integer, intent(in) :: nbnd !< NUMBER OF WAVELENGTH INTERVALS FOR WHICH TO COMPUTE THE ALBEDOS
  !
  ! OUTPUT
  !
  real, intent(out) :: albdif (ilg,nbnd) !< DIFFUSE SNOW ALBEDO (AKA WHITE SKY ALBEDO)
  real, intent(out) :: albdir(ilg,nbnd) !< DIRECT BEAM SNOW ALBEDO (AKA BLACK SKY ALBEDO)
  !
  ! LOCAL
  !
  real, DIMENSION(ILG,2)  :: wsmu, wbc, wreff, wswe
  real                    :: wsalb(2), wtt, snow_reff_l
  integer, DIMENSION(ILG) :: ismu, ibc, ireff, iswe
  integer                 :: ib, i, isalb, iismu, iisalb, iibc, iireff, iiswe, mvidx
  !
  ! CONSTANTS
  !
  integer,PARAMETER :: &
                       nsmu     = 10, &
                       nsalb    = 11, &
                       nbc      = 20, &
                       nreff    = 10, &
                       nswe     = 11, &
                       nbnd_lut = 4

  real,PARAMETER :: & ! STATE VALUES FOR LUT
                    LSALB(NSALB)      = (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/), &
                    LSMU(NSMU)        = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/), &
                    LSNOW_REFF(NREFF) = (/50.0,75.0,100.0,150.0,200.0,275.0,375.0,500.0,700.0,1000.0/), &
                    LSWE(NSWE)        = (/0.1,0.25,0.65,1.7,4.4,12.0,30.0,75.0,200.0,500.0,5000.0/), &
                    LBC_CONC(NBC)     = (/0.0,1.0, &
                    5.0,10.0, &
                    50.0,100.0, &
                    500.0,1000.0, &
                    5000.0,10000.0, &
                    50000.0,100000.0, &
                    250000.0,500000.0,750000.0,1000000.0, &
                    2500000.0,5000000.0,7500000.0,10000000.0/)

  real,DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) :: &
                                                        albdif_lut, &
                                                        albdir_lut

  integer :: &
             snow_alb_lut_init

  COMMON /SNOWALBLUT/ albdif_lut,albdir_lut,snow_alb_lut_init

  ! ABORT IF THE LUT HAS NOT BEEN READ IN
  if (snow_alb_lut_init /= 1) then
    write(6, * ) 'SNOW ALBEDO LUT HAS NOT BEEN INITIALIZED', &
               snow_alb_lut_init
    call errorHandler('SNOW_ALBVAL', - 1)
  end if

  ! ABORT IF THE NUMBER OF BANDS IN THE LUT DOES NOT MATCH SIZE PASSED IN
  if (nbnd_lut /= nbnd) then
    write(6, * ) 'MISMATCH IN NUMBER OF WAVELENGTH INTERVALS'
    call errorHandler('SNOW_ALBVAL', - 2)
  end if

  ! COMPUTE THE ALBEDOS USING LINEAR INTERPOLATION

  ! COMPUTE THE INTERPOLATION WEIGHTS AND POINTS ONCE AND REUSE FOR
  ! ALBEDO INTERPOLATION FOR EACH BAND.
  ! HAVE A CHECK TO SET THE WEIGHTS DEPENDING IF THE INPUT IS
  ! OUTSIDE OR INSIDE THE LOOKUP TABLE RANGE

  do i = il1,il2
    if (c_ind(i) == 1) then
      snow_reff_l = snow_reff(i)
      ismu(i)     = mvidx(LSMU,nsmu,smu(i))
      ibc(i)      = mvidx(LBC_CONC,nbc,bc_conc(i))
      ireff(i)    = mvidx(LSNOW_REFF,nreff,snow_reff(i))
      iswe(i)     = mvidx(LSWE,nswe,swe(i))

      if (smu(i) <= LSMU(1)) then
        wsmu(i,2) = 0.0
        wsmu(i,1) = 1.0 - wsmu(i,2)
      else if (smu(i) > LSMU(NSMU)) then
        wsmu(i,2) = 1.0
        wsmu(i,1) = 1.0 - wsmu(i,2)
      else
        wsmu(i,2) = (smu(i) - LSMU(ismu(i))) &
                    / (LSMU(ismu(i) + 1) - LSMU(ismu(i)))
        wsmu(i,1) = 1.0 - wsmu(i,2)
      end if

      if (bc_conc(i) <= LBC_CONC(1)) then
        wbc(i,2) = 0.0
        wbc(i,1) = 1.0 - wbc(i,2)
      else if (bc_conc(i) > LBC_CONC(NBC)) then
        wbc(i,2) = 1.0
        wbc(i,1) = 1.0 - wbc(i,2)
      else
        wbc(i,2) = (bc_conc(i) - LBC_CONC(ibc(i))) &
                   / (LBC_CONC(ibc(i) + 1) - LBC_CONC(ibc(i)))
        wbc(i,1) = 1.0 - wbc(i,2)
      end if

      if (snow_reff_l <= LSNOW_REFF(1)) then
        wreff(i,2) = 0.0
        wreff(i,1) = 1.0 - wreff(i,2)
      else if (snow_reff_l > LSNOW_REFF(NREFF)) then
        wreff(i,2) = 1.0
        wreff(i,1) = 1.0 - wreff(i,2)
      else
        wreff(i,2) = (snow_reff_l - LSNOW_REFF(ireff(i))) &
                     / (LSNOW_REFF(ireff(i) + 1) &
                     - LSNOW_REFF(ireff(i)))
        wreff(i,1) = 1.0 - wreff(i,2)
      end if

      if (swe(i) <= LSWE(1)) then
        wswe(i,2) = 0.0
        wswe(i,1) = 1.0 - wswe(i,2)
      else if (swe(i) > LSWE(NSWE)) then
        wswe(i,2) = 1.0
        wswe(i,1) = 1.0 - wswe(i,2)
      else
        wswe(i,2) = (swe(i) - LSWE(iswe(i))) &
                    / (LSWE(iswe(i) + 1) - LSWE(iswe(i)))
        wswe(i,1) = 1.0 - wswe(i,2)
      end if
    end if
  end do ! i

  do ib = 1,nbnd
    do i = il1,il2
      if (c_ind(i) == 1) then

        isalb = mvidx(LSALB,nsalb,salb(i,ib))

        if (salb(i,ib) <= LSALB(1)) then
          wsalb(2) = 0.0
          wsalb(1) = 1.0 - wsalb(2)
        else if (salb(i,ib) > LSALB(NSALB)) then
          wsalb(2) = 1.0
          wsalb(1) = 1.0 - wsalb(2)
        else
          wsalb(2) = (salb(i,ib) - LSALB(isalb)) &
                     / (LSALB(isalb + 1) - LSALB(isalb))
          wsalb(1) = 1.0 - wsalb(2)
        end if

        albdir(i,ib) = 0.0
        albdif (i,ib) = 0.0

        do iisalb = isalb,isalb + 1
          do iismu = ismu(i),ismu(i) + 1
            do iireff = ireff(i),ireff(i) + 1
              do iiswe = iswe(i),iswe(i) + 1
                do iibc = ibc(i),ibc(i) + 1

                  wtt = wsmu(i,iismu - ismu(i) + 1) &
                        * wreff(i,iireff - ireff(i) + 1) &
                        * wswe(i,iiswe - iswe(i) + 1) &
                        * wbc(i,iibc - ibc(i) + 1) &
                        * wsalb(iisalb - isalb + 1)

                  albdif (i,ib) = albdif (i,ib) + wtt &
                                  * albdif_lut(iibc,iiswe,iireff,iismu,iisalb,ib)
                  albdir(i,ib) = albdir(i,ib) + wtt &
                                 * albdir_lut(iibc,iiswe,iireff,iismu,iisalb,ib)

                end do ! iibc
              end do  ! iiswe
            end do     ! iireff
          end do        ! iismu
        end do           ! iisalb

        if (albdif (i,ib) > 1.0 .or. albdif (i,ib) < 0.0) then
          write(6, * ) 'Bad albdif ',i,ib,smu(i),bc_conc(i), &
                      snow_reff(i),swe(i),salb(i,ib),albdif (i,ib)
          write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
          call errorHandler('SNOW_ALBVAL', - 3)
        end if
        if (albdir(i,ib) > 1.0 .or. albdir(i,ib) < 0.0) then
          write(6, * ) 'Bad albdir ',i,ib,smu(i),bc_conc(i), &
                        snow_reff(i),swe(i),salb(i,ib),albdir(i,ib)
          write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
          call errorHandler('SNOW_ALBVAL', - 3)
        end if
      else
        albdif (i,ib) = - 999.0
        albdir(i,ib) = - 999.0
      end if
    end do ! i
  end do ! ib

  return
end subroutine SNOW_ALBVAL
!> \file
!> This subroutine computes the direct and diffuse snow albedo using a
!! lookup table and information about the current snow pack state.
!! Albedos are computed for each solar radiation wavelength intervals
!! so a total of 8 albedos will be returned.  These albedos can then be
!! used to compute the total snow albedo based on the by weighting
!! the results by the direct beam fraction of the incident solar radiation.
