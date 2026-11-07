!> \file
!> Transfers information between the 'gathered' and 'scattered' form of the CTEM data arrays.
!! @author R. Li, J. Melton, E. Chan
!!
module ctemGatherScatter

  implicit none

  public :: ctemg1
  public :: ctems2
  public :: ctemg2

contains

  !> \ingroup ctemgatherscatter_ctemg1
  !! @{
  !> Performs initial 'gather' operation on CTEM variables for consistency
  !! with physics variables gather operations.
  !!
  !! ctemg1 converts variables from the 'row' format (nlat, nmos, ...)
  !! to the 'gat' format (ilg, ...) which is what the model calculations
  !! are performed on. The ctemg1 subroutine is used to transform the
  !! read in state variables (which come in with the 'row' format from the
  !! various input files).
  !! @author R. Li, J. Melton
  subroutine ctemg1 (gleafmasgat, bleafmasgat, stemmassgat, &  ! Out
                     rootmassgat, fcancmxgat, zbtwgat, & ! Out
                     dlzwgat, sdepgat, ailcggat, & ! Out
                     ailcbgat, ailcgat, zolncgat, rmatcgat, & ! Out
                     rmatctemgat, slaigat, bmasveggat, cmasvegcgat, & ! Out
                     veghghtgat, rootdpthgat, alvsctmgat, alirctmgat, & ! Out
                     paicgat, slaicgat, faregat, & ! Out
                     ipeatlandgat, maxAnnualActLyrGAT, & ! Out
                     tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, & ! Out
                     tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, & ! Out
                     tracerMossCMassgat, tracerMossLitrMassgat, & ! Out                     
                     twarmmgat, tcoldmgat, gdd5gat, & ! Out
                     ariditygat, srplsmongat, defctmongat, anndefctgat, & ! Out
                     annsrplsgat, annpcpgat, dry_season_lengthgat, & ! Out
                     litrmsmossgat, Cmossmasgat, dmossgat, & ! Out                     
                     pandaysgat, lfstatusgat, slopefracgat, pstemmassgat, & ! Out
                     pgleafmassgat,litrmassgat, soilcmasgat, grwtheffgat, & ! Out
                     ilmos, jlmos, iwmos, jwmos, nml, &! In
                     gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, &! In
                     fcancmxrow, zbtwrot, dlzwrot, sdeprot, &! In
                     ailcgrow, ailcbrow, ailcrow, zolncrow, &! In
                     rmatcrow, rmatctemrow, slairow, bmasvegrow, &! In
                     cmasvegcrow, veghghtrow, rootdpthrow, alvsctmrow, &! In
                     alirctmrow, paicrow, slaicrow, FAREROT, &! In
                     ipeatlandrow, maxAnnualActLyrROT, &! In
                     tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &! In
                     tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &! In
                     tracerMossCMassrot, tracerMossLitrMassrot, & ! In
                     twarmmrow, tcoldmrow, gdd5row, & ! In
                     aridityrow, srplsmonrow, defctmonrow, anndefctrow, & ! In
                     annsrplsrow, annpcprow, dry_season_lengthrow, & ! In
                     litrmsmossrow, Cmossmasrow, dmossrow, & ! In
                     pandaysrow, lfstatusrow, slopefracrow, pstemmassrow, & ! In
                     pgleafmassrow, litrmassrow, soilcmasrow, grwtheffrow)! In

    !
    !     2 May 2019  - Convert to f90 and put in this module.
    !     J . Melton
    !     22  Jul 2013  - Add in module for parameters
    !     J. Melton
    !
    !     July 52009   - gather operation on ctem variables for consistency
    !                   with class' tiled version
    !     Rong Li
    !
    use classicParams, only : nlat, nmos, ilg, ignd, ican, icp1, icc, iccp2

    implicit none

    real, intent(out) :: gleafmasgat(ilg,icc)
    real, intent(out) :: bleafmasgat(ilg,icc)
    real, intent(out) :: stemmassgat(ilg,icc)
    real, intent(out) ::   rootmassgat(ilg,icc)
    real, intent(out) :: fcancmxgat(ilg,icc)
    real, intent(out) :: zbtwgat(ilg,ignd)
    real, intent(out) :: dlzwgat(ilg,ignd)
    real, intent(out) :: sdepgat(ilg)
    real, intent(out) :: ailcggat(ilg,icc)
    real, intent(out) :: ailcbgat(ilg,icc)
    real, intent(out) :: ailcgat(ilg,ican)
    real, intent(out) :: zolncgat(ilg,ican)
    real, intent(out) :: rmatcgat(ilg,ican,ignd)
    real, intent(out) :: faregat(ilg)
    real, intent(out) :: rmatctemgat(ilg,icc,ignd)
    real, intent(out) :: slaigat(ilg,icc)
    real, intent(out) :: bmasveggat(ilg,icc)
    real, intent(out) :: cmasvegcgat(ilg,ican)
    real, intent(out) :: veghghtgat(ilg,icc)
    real, intent(out) :: rootdpthgat(ilg,icc)
    real, intent(out) :: alvsctmgat(ilg,ican)
    real, intent(out) :: alirctmgat(ilg,ican)
    real, intent(out) :: paicgat(ilg,ican)
    real, intent(out) :: slaicgat(ilg,ican)
    real, intent(out) :: maxAnnualActLyrGAT(ilg)
    integer, intent(out) :: ipeatlandgat(ilg)
    real, intent(out) :: tracermossCMassgat(ilg)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(out) :: tracermossLitrMassgat(ilg)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(out) :: tracergLeafMassgat(ilg,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerbLeafMassgat(ilg,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerstemMassgat(ilg,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerrootMassgat(ilg,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerlitrMassgat(ilg,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: tracersoilCMassgat(ilg,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: twarmmgat(ilg)
    real, intent(out) :: tcoldmgat(ilg)
    real, intent(out) :: gdd5gat(ilg)
    real, intent(out) :: ariditygat(ilg)
    real, intent(out) :: srplsmongat(ilg)
    real, intent(out) :: defctmongat(ilg)
    real, intent(out) :: anndefctgat(ilg) 
    real, intent(out) :: annsrplsgat(ilg)
    real, intent(out) :: annpcpgat(ilg)
    real, intent(out) :: dry_season_lengthgat(ilg)
    real, intent(out) :: litrmsmossgat(ilg)
    real, intent(out) :: Cmossmasgat(ilg)
    real, intent(out) :: dmossgat(ilg)
    integer, intent(out) :: pandaysgat(ilg,icc)
    integer, intent(out) :: lfstatusgat(ilg,icc)
    real, intent(out) :: slopefracgat(ilg,8)
    real, intent(out) :: pstemmassgat(ilg,icc)
    real, intent(out) :: pgleafmassgat(ilg,icc)
    ! COMBAK PERLAY
    real, intent(out)  :: soilcmasgat(ilg,iccp2)
    real, intent(out)  :: litrmassgat(ilg,iccp2)
    !real, intent(out)  :: soilcmasgat(ilg,iccp2,ignd)
    !real, intent(out)  :: litrmassgat(ilg,iccp2,ignd)
    ! COMBAK PERLAY
    real, intent(out) :: grwtheffgat(ilg,icc)

    integer, intent(in) :: ilmos (ilg)
    integer, intent(in) :: jlmos  (ilg)
    integer, intent(in) :: iwmos  (ilg)
    integer, intent(in) :: jwmos (ilg)
    integer, intent(in) :: nml
    real, intent(in) :: gleafmasrow(nlat,nmos,icc)
    real, intent(in) :: bleafmasrow(nlat,nmos,icc)
    real, intent(in) :: stemmassrow(nlat,nmos,icc)
    real, intent(in) :: rootmassrow(nlat,nmos,icc)
    real, intent(in) :: fcancmxrow(nlat,nmos,icc)
    real, intent(in) :: zbtwrot(nlat,nmos,ignd)
    real, intent(in) :: dlzwrot(nlat,nmos,ignd)
    real, intent(in) :: sdeprot(nlat,nmos)
    real, intent(in) :: ailcgrow(nlat,nmos,icc)
    real, intent(in) :: ailcbrow(nlat,nmos,icc)
    real, intent(in) :: ailcrow(nlat,nmos,ican)
    real, intent(in) :: zolncrow(nlat,nmos,ican)
    real, intent(in) :: rmatcrow(nlat,nmos,ican,ignd)
    real, intent(in) :: FAREROT(nlat,nmos)
    real, intent(in) :: rmatctemrow(nlat,nmos,icc,ignd)
    real, intent(in) :: slairow(nlat,nmos,icc)
    real, intent(in) :: bmasvegrow(nlat,nmos,icc)
    real, intent(in) :: cmasvegcrow(nlat,nmos,ican)
    real, intent(in) :: veghghtrow(nlat,nmos,icc)
    real, intent(in) :: rootdpthrow(nlat,nmos,icc)
    real, intent(in) :: alvsctmrow(nlat,nmos,ican)
    real, intent(in) :: alirctmrow(nlat,nmos,ican)
    real, intent(in) :: paicrow(nlat,nmos,ican)
    real, intent(in) :: slaicrow(nlat,nmos,ican)
    real, intent(in) :: maxAnnualActLyrROT(nlat,nmos)
    integer, intent(in) ::  ipeatlandrow(nlat,nmos)
    real, intent(in) :: tracermossCMassrot(nlat,nmos)     !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(in) :: tracermossLitrMassrot(nlat,nmos)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: tracergLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerbLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerstemMassrot(nlat,nmos,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerrootMassrot(nlat,nmos,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerlitrMassrot(nlat,nmos,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracersoilCMassrot(nlat,nmos,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in)  :: twarmmrow(nlat,nmos)
    real, intent(in)  :: tcoldmrow(nlat,nmos)
    real, intent(in)  :: gdd5row(nlat,nmos)
    real, intent(in)  :: aridityrow(nlat,nmos)
    real, intent(in)  :: srplsmonrow(nlat,nmos)
    real, intent(in)  :: defctmonrow(nlat,nmos)
    real, intent(in)  :: anndefctrow(nlat,nmos)
    real, intent(in)  :: annsrplsrow(nlat,nmos)
    real, intent(in)  :: annpcprow(nlat,nmos)                         
    real, intent(in)  :: dry_season_lengthrow(nlat,nmos)
    real, intent(in)  :: litrmsmossrow(nlat,nmos)
    real, intent(in)  :: Cmossmasrow(nlat,nmos)
    real, intent(in)  :: dmossrow(nlat,nmos)
    integer, intent(in) :: pandaysrow(nlat,nmos,icc)
    integer, intent(in) :: lfstatusrow(nlat,nmos,icc)
    real, intent(in)  :: slopefracrow(nlat,nmos,8)
    real, intent(in)  :: pstemmassrow(nlat,nmos,icc)
    real, intent(in)  :: pgleafmassrow(nlat,nmos,icc)
    ! COMBAK PERLAY
    real, intent(in)  :: soilcmasrow(nlat,nmos,iccp2)
    real, intent(in)  :: litrmassrow(nlat,nmos,iccp2)
    !real, intent(in)  :: soilcmasrow(nlat,nmos,iccp2,ignd)
    !real, intent(in)  :: litrmassrow(nlat,nmos,iccp2,ignd)
    ! COMBAK PERLAY
    real, intent(in) :: grwtheffrow(nlat,nmos,icc)

    ! Local
    integer ::  k, l, m, n

    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      sdepgat(k) = sdeprot(ilmos(k),jlmos(k))
      ipeatlandgat(k) = ipeatlandrow(ilmos(k),jlmos(k))
      faregat(k) = FAREROT(ilmos(k),jlmos(k))
      maxAnnualActLyrGAT(k) = maxAnnualActLyrROT(ilmos(k),jlmos(k))
      tracerMossCMassgat(k) = tracerMossCMassrot(ilmos(k),jlmos(k))
      tracermossLitrMassgat(k) = tracermossLitrMassrot(ilmos(k),jlmos(k))
      twarmmgat(k)    = twarmmrow(ilmos(k),jlmos(k))
      tcoldmgat(k)    = tcoldmrow(ilmos(k),jlmos(k))
      gdd5gat(k)      = gdd5row(ilmos(k),jlmos(k))
      ariditygat(k)   = aridityrow(ilmos(k),jlmos(k))
      srplsmongat(k)  = srplsmonrow(ilmos(k),jlmos(k))
      defctmongat(k)  = defctmonrow(ilmos(k),jlmos(k))
      anndefctgat(k)  = anndefctrow(ilmos(k),jlmos(k))
      annsrplsgat(k)  = annsrplsrow(ilmos(k),jlmos(k))
      annpcpgat(k)    = annpcprow(ilmos(k),jlmos(k))
      dry_season_lengthgat(k) = &
                                dry_season_lengthrow(ilmos(k),jlmos(k))
      litrmsmossgat(k) =  litrmsmossrow(ilmos(k),jlmos(k))
      Cmossmasgat(k) = Cmossmasrow(ilmos(k),jlmos(k))
      dmossgat(k) = dmossrow(ilmos(k),jlmos(k))
      ipeatlandgat(k) = ipeatlandrow(ilmos(k),jlmos(k))
      do n = 1,8
        slopefracgat(k,n) = slopefracrow(ilmos(k),jlmos(k),n)
      end do
    end do ! loop 100

    do l = 1,icc ! loop 101
      do k = 1,nml
        gleafmasgat(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
        bleafmasgat(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
        stemmassgat(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
        rootmassgat(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
        fcancmxgat(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
        ailcggat(k,l)    = ailcgrow(ilmos(k),jlmos(k),l)
        ailcbgat(k,l)    = ailcbrow(ilmos(k),jlmos(k),l)
        slaigat(k,l)     = slairow(ilmos(k),jlmos(k),l)
        bmasveggat(k,l)  = bmasvegrow(ilmos(k),jlmos(k),l)
        veghghtgat(k,l)  = veghghtrow(ilmos(k),jlmos(k),l)
        rootdpthgat(k,l) = rootdpthrow(ilmos(k),jlmos(k),l)
        tracergLeafMassgat(k,l) = tracergLeafMassrot(ilmos(k),jlmos(k),l)
        tracerbLeafMassgat(k,l) = tracerbLeafMassrot(ilmos(k),jlmos(k),l)
        tracerStemMassgat(k,l) = tracerStemMassrot(ilmos(k),jlmos(k),l)
        tracerRootMassgat(k,l) = tracerRootMassrot(ilmos(k),jlmos(k),l)
        pandaysgat(k,l)  = pandaysrow(ilmos(k),jlmos(k),l)
        lfstatusgat(k,l) = lfstatusrow(ilmos(k),jlmos(k),l)
        pstemmassgat(k,l) = pstemmassrow(ilmos(k),jlmos(k),l)
        pgleafmassgat(k,l) = pgleafmassrow(ilmos(k),jlmos(k),l)
        grwtheffgat(k,l) = grwtheffrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 101

    do l = 1,ican ! loop 201
      do k = 1,nml
        ailcgat(k,l)     = ailcrow(ilmos(k),jlmos(k),l)
        zolncgat(k,l)    = zolncrow(ilmos(k),jlmos(k),l)
        cmasvegcgat(k,l) = cmasvegcrow(ilmos(k),jlmos(k),l)
        alvsctmgat(k,l)  = alvsctmrow(ilmos(k),jlmos(k),l)
        alirctmgat(k,l)  = alirctmrow(ilmos(k),jlmos(k),l)
        paicgat(k,l)     = paicrow(ilmos(k),jlmos(k),l)
        slaicgat(k,l)    = slaicrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 201

    do l = 1,ignd ! loop 250
      do k = 1,nml
        zbtwgat(k,l) = zbtwrot(ilmos(k),jlmos(k),l)
        dlzwgat(k,l) = dlzwrot(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 250

    do l = 1,icc ! loop 280
      do m = 1,ignd
        do k = 1,nml
          rmatctemgat(k,l,m) = rmatctemrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 280

    do l = 1,ican ! loop 290
      do m = 1,ignd
        do k = 1,nml
          rmatcgat(k,l,m) = rmatcrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 290

    do l = 1,iccp2
      do k = 1,nml
        ! COMBAK PERLAY
        do m = 1,ignd
          tracerSoilCMassgat(k,l,m) = tracerSoilCMassrot(ilmos(k),jlmos(k),l,m)
          tracerLitrMassgat(k,l,m) = tracerLitrMassrot(ilmos(k),jlmos(k),l,m)
        !litrmassgat(k,l,m) = litrmassrow(ilmos(k),jlmos(k),l,m)
        !soilcmasgat(k,l,m) = soilcmasrow(ilmos(k),jlmos(k),l,m)
        end do
        litrmassgat(k,l) = litrmassrow(ilmos(k),jlmos(k),l)
        soilcmasgat(k,l) = soilcmasrow(ilmos(k),jlmos(k),l)
        ! COMBAK PERLAY
      end do
    end do
    return

  end subroutine ctemg1
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup ctemgatherscatter_ctems2
  !! @{
  !> Performs subsequent scatter operation on biogeochemical variables
  !!
  !! ctems2 converts variables from the 'gat' format to the
  !! 'row' format, which is suitable for writing to output/restart
  !! files. If a variable is not written to either of those files,
  !! there is no need to scatter the variable as it will be in the
  !! correct format for model calclations ('gat').
  !!
  !! @author R. Li, J. Melton, E. Chan
  subroutine ctems2 (fcancmxrow, rmatcrow, zolncrow, paicrow, &
                     ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                     ailcgsrow, fcancsrow, fcancrow, rmatctemrow, &
                     co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                     co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                     cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                     rmlcgvegrow, canresrow, sdeprow, ch4concrow, &
                     sandrow, clayrow, orgmrow, &
                     anvegrow, rmlvegrow, tbaraccrow_m, prbfrhucgrd, &
                     extnprobgrd, pfcancmxrow, nfcancmxrow, &
                     stemmassrow, rootmassrow, litrmassrow, gleafmasrow, &
                     bleafmasrow, soilcmasrow, ailcbrow, flhrlossrow, &
                     pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                     lyrotmasrow, tymaxlairow, vgbiomasrow, gavgltmsrow, &
                     stmhrlosrow, bmasvegrow, colddaysrow, rothrlosrow, &
                     alvsctmrow, alirctmrow, gavglairow, npprow, &
                     neprow, hetroresrow, autoresrow, soilresprow, &
                     rmrow, rgrow, nbprow, litresrow, &
                     socresrow, gpprow, dstcemlsrow, litrfallrow, &
                     humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                     litrfallvegrow, humiftrsvegrow, &
                     rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                     tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                     afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                     burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                     lucsocinrow, nppvegrow, dstcemls3row, &
                     farerow, gavgscmsrow, &
                     rmlvegaccrow, rmsvegrow, rmrvegrow, rgvegrow, &
                     vgbiomas_vegrow, gppvegrow, nepvegrow, &
                     fcanrow, pftexistrow, &
                     emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                     emit_h2row, emit_noxrow, emit_n2orow, emit_pm25row, &
                     emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                     btermrow, ltermrow, mtermrow, &
                     nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                     soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                     ch4WetSpecrow, wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                     twarmmrow, tcoldmrow, gdd5row, &
                     aridityrow, srplsmonrow, defctmonrow, anndefctrow, &
                     annsrplsrow, annpcprow, dry_season_lengthrow, &
                     anmossrow, rmlmossrow, gppmossrow, armossrow, &
                     nppmossrow, peatdeprow, litrmsmossrow, Cmossmasrow, &
                     dmossrow, ipeatlandrow, pddrow, wetfrac_presrow, & ! thlqaccrow_m, thicaccrow_m, &
                     tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &
                     tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &
                     tracerMossCMassrot, tracerMossLitrMassrot, &
                     ilmos, jlmos, iwmos, jwmos, &
                     nml, fcancmxgat, rmatcgat, zolncgat, paicgat, &
                     ailcgat, ailcggat, cmasvegcgat, slaicgat, &
                     ailcgsgat, fcancsgat, fcancgat, rmatctemgat, &
                     co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                     co2i2csgat, xdiffusgat, slaigat, cfluxcggat, &
                     cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                     rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                     sandgat, claygat, orgmgat, &
                     anveggat, rmlveggat, tbaraccgat_m, prbfrhucgat, &
                     extnprobgat, pfcancmxgat, nfcancmxgat, &
                     stemmassgat, rootmassgat, litrmassgat, gleafmasgat, &
                     bleafmasgat, soilcmasgat, ailcbgat, flhrlossgat, &
                     pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                     lyrotmasgat, tymaxlaigat, vgbiomasgat, gavgltmsgat, &
                     stmhrlosgat, bmasveggat, colddaysgat, rothrlosgat, &
                     alvsctmgat, alirctmgat, gavglaigat, nppgat, &
                     nepgat, hetroresgat, autoresgat, soilrespgat, &
                     rmgat, rggat, nbpgat, litresgat, &
                     socresgat, gppgat, dstcemlsgat, litrfallgat, &
                     humiftrsgat, veghghtgat, rootdpthgat, rmlgat, &
                     litrfallveggat, humiftrsveggat, &
                     rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                     tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                     afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                     burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                     lucsocingat, nppveggat, dstcemls3gat, &
                     faregat, gavgscmsgat, &
                     rmlvegaccgat, rmsveggat, rmrveggat, rgveggat, &
                     vgbiomas_veggat, gppveggat, nepveggat, &
                     fcangat, pftexistgat, &
                     emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                     emit_h2gat, emit_noxgat, emit_n2ogat, emit_pm25gat, &
                     emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                     btermgat, ltermgat, mtermgat, &
                     nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                     soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                     ch4WetSpecgat, wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                     twarmmgat, tcoldmgat, gdd5gat, &
                     ariditygat, srplsmongat, defctmongat, anndefctgat, &
                     annsrplsgat, annpcpgat, dry_season_lengthgat, &
                     anmossgat, rmlmossgat, gppmossgat, armossgat, &
                     nppmossgat, peatdepgat, litrmsmossgat, Cmossmasgat, &
                     dmossgat, ipeatlandgat, pddgat, wetfrac_presgat, &
                     tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, &
                     tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, &
                     tracerMossCMassgat, tracerMossLitrMassgat)!, thlqaccgat_m, thicaccgat_m)

    !     Dec 232016     Remove thlqaccXXX_m/thicaccXXX_m
    !     Ed Chan
    !
    !     July 122013    Bring in the ctem params use statement
    !     J. Melton

    !      August 4, 2009 scatter operation on CTEM variables.
    !      Rong Li
    !
    use classicParams, only : nlat, nmos, ilg, ignd, ican, icp1, &
                                icc, iccp2, iccp1

    implicit none

    !
    !
    real, intent(out) :: fcancmxrow(nlat,nmos,icc), rmatcrow(nlat,nmos,ican,ignd), &
                         zolncrow(nlat,nmos,ican), paicrow(nlat,nmos,ican), &
                         ailcrow(nlat,nmos,ican), ailcgrow(nlat,nmos,icc), &
                         cmasvegcrow(nlat,nmos,ican), slaicrow(nlat,nmos,ican), &
                         ailcgsrow(nlat,nmos,icc), fcancsrow(nlat,nmos,icc), &
                         fcancrow(nlat,nmos,icc), rmatctemrow(nlat,nmos,icc,ignd), &
                         co2concrow(nlat,nmos), co2i1cgrow(nlat,nmos,icc), &
                         co2i1csrow(nlat,nmos,icc), co2i2cgrow(nlat,nmos,icc), &
                         co2i2csrow(nlat,nmos,icc), xdiffus(nlat), &
                         slairow(nlat,nmos,icc), cfluxcgrow(nlat,nmos), &
                         cfluxcsrow(nlat,nmos), ancsvegrow(nlat,nmos,icc), &
                         ancgvegrow(nlat,nmos,icc), rmlcsvegrow(nlat,nmos,icc), &
                         rmlcgvegrow(nlat,nmos,icc), canresrow(nlat,nmos), &
                         sdeprow(nlat,nmos), fcanrow(nlat,nmos,icp1), &
                         ch4concrow(nlat,nmos), wetfrac_presrow(nlat,nmos)
    !
    real, intent(out) :: sandrow(nlat,nmos,ignd), clayrow(nlat,nmos,ignd), &
                         orgmrow(nlat,nmos,ignd)
    !
    real, intent(out) :: anvegrow(nlat,nmos,icc), rmlvegrow(nlat,nmos,icc)
    !
    real, intent(out) :: prbfrhucgrd(nlat), extnprobgrd(nlat), &
                         tbaraccrow_m(nlat,nmos,ignd), &
                         pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc), &
                         stemmassrow(nlat,nmos,icc), rootmassrow(nlat,nmos,icc), &
                         pstemmassrow(nlat,nmos,icc), pgleafmassrow(nlat,nmos,icc), &
                         gleafmasrow(nlat,nmos,icc), bleafmasrow(nlat,nmos,icc), &
                         ailcbrow(nlat,nmos,icc), flhrlossrow(nlat,nmos,icc), &
                         ! COMBAK PERLAY
                         litrmassrow(nlat,nmos,iccp2), &
                         soilcmasrow(nlat,nmos,iccp2)
    ! &      litrmassrow(nlat,nmos,iccp2,ignd), &
    ! &      soilcmasrow(nlat,nmos,iccp2,ignd), &
    ! COMBAK PERLAY

    !
    integer, intent(out) :: pandaysrow(nlat,nmos,icc), lfstatusrow(nlat,nmos,icc), &
                            colddaysrow(nlat,nmos,2)

    logical, intent(out) :: pftexistrow(nlat,nmos,icc)
    !
    real, intent(out) :: grwtheffrow(nlat,nmos,icc), lystmmasrow(nlat,nmos,icc), &
                         lyrotmasrow(nlat,nmos,icc), tymaxlairow(nlat,nmos,icc), &
                         vgbiomasrow(nlat,nmos), gavgltmsrow(nlat,nmos), &
                         stmhrlosrow(nlat,nmos,icc), bmasvegrow(nlat,nmos,icc), &
                         rothrlosrow(nlat,nmos,icc), &
                         alvsctmrow(nlat,nmos,ican), alirctmrow(nlat,nmos,ican), &
                         gavglairow(nlat,nmos)
    !
    real, intent(out) :: npprow(nlat,nmos), neprow(nlat,nmos), &
                         hetroresrow(nlat,nmos), autoresrow(nlat,nmos), &
                         soilresprow(nlat,nmos), rmrow(nlat,nmos), rgrow(nlat,nmos), &
                         nbprow(nlat,nmos), litresrow(nlat,nmos), &
                         socresrow(nlat,nmos), gpprow(nlat,nmos), &
                         dstcemlsrow(nlat,nmos), litrfallrow(nlat,nmos), &
                         humiftrsrow(nlat,nmos), veghghtrow(nlat,nmos,icc), &
                         litrfallvegrow(nlat,nmos,icc), &
                         ! COMBAK PERLAY
                         humiftrsvegrow(nlat,nmos,iccp2), &
                         ! &      humiftrsvegrow(nlat,nmos,iccp2,ignd), &
                         ! COMBAK PERLAY
                         rootdpthrow(nlat,nmos,icc), rmlrow(nlat,nmos), &
                         rmsrow(nlat,nmos), rmrrow(nlat,nmos), &
                         tltrleafrow(nlat,nmos,icc), tltrstemrow(nlat,nmos,icc), &
                         tltrrootrow(nlat,nmos,icc), leaflitrrow(nlat,nmos,icc), &
                         roottemprow(nlat,nmos,icc), afrleafrow(nlat,nmos,icc), &
                         afrstemrow(nlat,nmos,icc), afrrootrow(nlat,nmos,icc), &
                         wtstatusrow(nlat,nmos,icc), ltstatusrow(nlat,nmos,icc), &
                         burnfracrow(nlat,nmos), smfuncvegrow(nlat,nmos,icc), &
                         lucemcomrow(nlat,nmos), lucltrinrow(nlat,nmos), &
                         lucsocinrow(nlat,nmos), nppvegrow(nlat,nmos,icc), &
                         dstcemls3row(nlat,nmos)
    !
    !     fire variables
    !
    real, intent(out) :: emit_co2row(nlat,nmos,icc), emit_corow(nlat,nmos,icc), &
                         emit_ch4row(nlat,nmos,icc), emit_nmhcrow(nlat,nmos,icc), &
                         emit_h2row(nlat,nmos,icc), emit_noxrow(nlat,nmos,icc), &
                         emit_n2orow(nlat,nmos,icc), emit_pm25row(nlat,nmos,icc), &
                         emit_tpmrow(nlat,nmos,icc), emit_tcrow(nlat,nmos,icc), &
                         emit_ocrow(nlat,nmos,icc), emit_bcrow(nlat,nmos,icc), &
                         burnvegfrow(nlat,nmos,icc), btermrow(nlat,nmos,icc), &
                         ltermrow(nlat,nmos), mtermrow(nlat,nmos,icc)

    real, intent(out) :: farerow(nlat,nmos), gavgscmsrow(nlat,nmos)
    !
    real, intent(out) :: rmlvegaccrow(nlat,nmos,icc), rmsvegrow(nlat,nmos,icc), &
                         rmrvegrow(nlat,nmos,icc), rgvegrow(nlat,nmos,icc)
    !
    real, intent(out) :: vgbiomas_vegrow(nlat,nmos,icc)
    !
    real, intent(out) :: gppvegrow(nlat,nmos,icc), nepvegrow(nlat,nmos,iccp1), &
                         nbpvegrow(nlat,nmos,iccp1), hetroresvegrow(nlat,nmos,iccp1), &
                         autoresvegrow(nlat,nmos,icc), &
                         ! COMBAK PERLAY
                         litresvegrow(nlat,nmos,iccp2), &
                         soilcresvegrow(nlat,nmos,iccp2)
    ! &      litresvegrow(nlat,nmos,iccp2,ignd), &
    ! &      soilcresvegrow(nlat,nmos,iccp2,ignd)
    ! COMBAK PERLAY

    real, intent(in) :: fcancmxgat(ilg,icc), rmatcgat(ilg,ican,ignd), &
                        zolncgat(ilg,ican), paicgat(ilg,ican), &
                        ailcgat(ilg,ican), ailcggat(ilg,icc), &
                        cmasvegcgat(ilg,ican), slaicgat(ilg,ican), &
                        ailcgsgat(ilg,icc), fcancsgat(ilg,icc), &
                        fcancgat(ilg,icc), rmatctemgat(ilg,icc,ignd), &
                        co2concgat(ilg), co2i1cggat(ilg,icc), &
                        co2i1csgat(ilg,icc), co2i2cggat(ilg,icc), &
                        co2i2csgat(ilg,icc), xdiffusgat(ilg), &
                        slaigat(ilg,icc), cfluxcggat(ilg), &
                        cfluxcsgat(ilg), ancsveggat(ilg,icc), &
                        ancgveggat(ilg,icc), rmlcsveggat(ilg,icc), &
                        rmlcgveggat(ilg,icc), canresgat(ilg), &
                        sdepgat(ilg), fcangat(ilg,icp1), &
                        ch4concgat(ilg), wetfrac_presgat(ilg)
    !
    real, intent(in) :: sandgat(ilg,ignd), claygat(ilg,ignd), &
                        orgmgat(ilg,ignd)
    !
    real, intent(in) :: anveggat(ilg,icc), rmlveggat(ilg,icc)
    !
    real, intent(in) :: prbfrhucgat(ilg), extnprobgat(ilg), &
                        tbaraccgat_m(ilg,ignd), &
                        pfcancmxgat(ilg,icc), nfcancmxgat(ilg,icc), &
                        stemmassgat(ilg,icc), rootmassgat(ilg,icc), &
                        pstemmassgat(ilg,icc), pgleafmassgat(ilg,icc), &
                        gleafmasgat(ilg,icc), bleafmasgat(ilg,icc), &
                        ailcbgat(ilg,icc), flhrlossgat(ilg,icc), &
                        ! COMBAK PERLAY
                        litrmassgat(ilg,iccp2), &
                        soilcmasgat(ilg,iccp2)
    ! &      litrmassgat(ilg,iccp2,ignd), &
    ! &      soilcmasgat(ilg,iccp2,ignd), &
    ! COMBAK PERLAY

    integer, intent(in) :: pandaysgat(ilg,icc), lfstatusgat(ilg,icc), &
                           colddaysgat(ilg,2)

    logical, intent(in) :: pftexistgat(ilg,icc)
    !
    real, intent(in) :: grwtheffgat(ilg,icc), lystmmasgat(ilg,icc), &
                        lyrotmasgat(ilg,icc), tymaxlaigat(ilg,icc), &
                        vgbiomasgat(ilg), gavgltmsgat(ilg), &
                        stmhrlosgat(ilg,icc), bmasveggat(ilg,icc), &
                        rothrlosgat(ilg,icc), alvsctmgat(ilg,ican), &
                        alirctmgat(ilg,ican), gavglaigat(ilg)
    !
    real, intent(in) :: nppgat(ilg), nepgat(ilg), &
                        hetroresgat(ilg), autoresgat(ilg), &
                        soilrespgat(ilg), rmgat(ilg), rggat(ilg), &
                        nbpgat(ilg), litresgat(ilg), &
                        socresgat(ilg), gppgat(ilg), &
                        dstcemlsgat(ilg), litrfallgat(ilg), &
                        humiftrsgat(ilg), veghghtgat(ilg,icc), &
                        litrfallveggat(ilg,icc), &
                        ! COMBAK PERLAY
                        humiftrsveggat(ilg,iccp2), &
                        ! humiftrsveggat(ilg,iccp2,ignd), &
                        ! COMBAK PERLAY
                        rootdpthgat(ilg,icc), rmlgat(ilg), &
                        rmsgat(ilg), rmrgat(ilg), &
                        tltrleafgat(ilg,icc), tltrstemgat(ilg,icc), &
                        tltrrootgat(ilg,icc), leaflitrgat(ilg,icc), &
                        roottempgat(ilg,icc), afrleafgat(ilg,icc), &
                        afrstemgat(ilg,icc), afrrootgat(ilg,icc), &
                        wtstatusgat(ilg,icc), ltstatusgat(ilg,icc), &
                        burnfracgat(ilg), smfuncveggat(ilg,icc), &
                        lucemcomgat(ilg), lucltringat(ilg), &
                        lucsocingat(ilg), nppveggat(ilg,icc), &
                        dstcemls3gat(ilg)
    !
    !      fire variables
    real, intent(in) :: emit_co2gat(ilg,icc), emit_cogat(ilg,icc), &
                        emit_ch4gat(ilg,icc), emit_nmhcgat(ilg,icc), &
                        emit_h2gat(ilg,icc), emit_noxgat(ilg,icc), &
                        emit_n2ogat(ilg,icc), emit_pm25gat(ilg,icc), &
                        emit_tpmgat(ilg,icc), emit_tcgat(ilg,icc), &
                        emit_ocgat(ilg,icc), emit_bcgat(ilg,icc), &
                        burnvegfgat(ilg,icc), btermgat(ilg,icc), &
                        ltermgat(ilg), mtermgat(ilg,icc)

    !
    real, intent(in) :: faregat(ilg)
    real, intent(in) :: gavgscmsgat(ilg)
    !
    real, intent(in) :: rmlvegaccgat(ilg,icc), rmsveggat(ilg,icc), &
                        rmrveggat(ilg,icc), rgveggat(ilg,icc)
    !
    real, intent(in) :: vgbiomas_veggat(ilg,icc)
    !
    real, intent(in) :: gppveggat(ilg,icc), nepveggat(ilg,iccp1), &
                        nbpveggat(ilg,iccp1), hetroresveggat(ilg,iccp1), &
                        autoresveggat(ilg,icc), &
                        ! COMBAK PERLAY
                        litresveggat(ilg,iccp2), &
                        soilcresveggat(ilg,iccp2)
    !       litresveggat(ilg,iccp2,ignd), &
    ! &      soilcresveggat(ilg,iccp2,ignd)
    ! COMBAK PERLAY

    !   Methane related variables
    real, intent(out)  :: ch4WetSpecrow(nlat,nmos), wetfdynrow(nlat,nmos), &
                          ch4WetDynrow(nlat,nmos), ch4soillsrow(nlat,nmos)
    real, intent(in)   :: ch4WetSpecgat(ilg), wetfdyngat(ilg), &
                          ch4WetDyngat(ilg), ch4soillsgat(ilg)

    real, intent(out) :: twarmmrow(nlat,nmos), tcoldmrow(nlat,nmos), &
                         gdd5row(nlat,nmos), aridityrow(nlat,nmos),  &
                         srplsmonrow(nlat,nmos), defctmonrow(nlat,nmos),  &
                         anndefctrow(nlat,nmos), annsrplsrow(nlat,nmos),  &
                         annpcprow(nlat,nmos), dry_season_lengthrow(nlat,nmos)
                         
    real, intent(in) :: twarmmgat(ilg), tcoldmgat(ilg), gdd5gat(ilg), &
                        ariditygat(ilg), srplsmongat(ilg), defctmongat(ilg), &
                        anndefctgat(ilg), annsrplsgat(ilg), annpcpgat(ilg), &
                        dry_season_lengthgat(ilg)

    !   Peatland variables
    real, intent(out) :: anmossrow(nlat,nmos), rmlmossrow(nlat,nmos),  &
                         gppmossrow(nlat,nmos), armossrow(nlat,nmos),  &
                         nppmossrow(nlat,nmos), peatdeprow(nlat,nmos),  &
                         litrmsmossrow(nlat,nmos), Cmossmasrow(nlat,nmos),  &
                         dmossrow(nlat,nmos), pddrow(nlat,nmos)
                         
    real, intent(in) :: anmossgat(ilg), rmlmossgat(ilg), gppmossgat(ilg), &
                        armossgat(ilg), nppmossgat(ilg), peatdepgat(ilg), &
                        litrmsmossgat(ilg), Cmossmasgat(ilg), dmossgat(ilg), &
                        pddgat(ilg)

    integer, intent(out) :: ipeatlandrow(nlat,nmos)
    integer, intent(in) :: ipeatlandgat(ilg)

    ! allocated with nlat,nmos,...:
    real, intent(out) :: tracermossCMassrot(nlat,nmos)     !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(out) :: tracermossLitrMassrot(nlat,nmos)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(out) :: tracergLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerbLeafMassrot(nlat,nmos,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerstemMassrot(nlat,nmos,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerrootMassrot(nlat,nmos,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: tracerlitrMassrot(nlat,nmos,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(out) :: tracersoilCMassrot(nlat,nmos,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    real, intent(in) :: tracermossCMassgat(ilg)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(in) :: tracermossLitrMassgat(ilg)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: tracergLeafMassgat(ilg,icc)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerbLeafMassgat(ilg,icc)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerstemMassgat(ilg,icc)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerrootMassgat(ilg,icc)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerlitrMassgat(ilg,iccp2,ignd)       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracersoilCMassgat(ilg,iccp2,ignd)      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    
    integer, intent(in) :: nml
    integer, intent(in)  :: ilmos(ilg), jlmos(ilg), iwmos(ilg), jwmos(ilg) !     gather-scatter index arrays.
    !
    integer ::  k,l,m

    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      sdeprow(ilmos(k),jlmos(k))        = sdepgat(k)
      co2concrow(ilmos(k),jlmos(k))     = co2concgat(k)
      ch4concrow(ilmos(k),jlmos(k))     = ch4concgat(k)
      cfluxcgrow(ilmos(k),jlmos(k))     = cfluxcggat(k)
      cfluxcsrow(ilmos(k),jlmos(k))     = cfluxcsgat(k)
      canresrow(ilmos(k),jlmos(k))      = canresgat(k)
      xdiffus(ilmos(k))                 = xdiffusgat(k)
      prbfrhucgrd(ilmos(k))             = prbfrhucgat(k)
      extnprobgrd(ilmos(k))             = extnprobgat(k)
      vgbiomasrow(ilmos(k),jlmos(k))    = vgbiomasgat(k)
      gavgltmsrow(ilmos(k),jlmos(k))    = gavgltmsgat(k)
      gavglairow(ilmos(k),jlmos(k))     = gavglaigat(k)
      npprow(ilmos(k),jlmos(k))         = nppgat(k)
      neprow(ilmos(k),jlmos(k))         = nepgat(k)
      hetroresrow(ilmos(k),jlmos(k))    = hetroresgat(k)
      autoresrow(ilmos(k),jlmos(k))     = autoresgat(k)
      soilresprow(ilmos(k),jlmos(k))    = soilrespgat(k)
      rmrow(ilmos(k),jlmos(k))          = rmgat(k)
      rgrow(ilmos(k),jlmos(k))          = rggat(k)
      nbprow(ilmos(k),jlmos(k))         = nbpgat(k)
      litresrow(ilmos(k),jlmos(k))      = litresgat(k)
      socresrow(ilmos(k),jlmos(k))      = socresgat(k)
      gpprow(ilmos(k),jlmos(k))         = gppgat(k)
      dstcemlsrow(ilmos(k),jlmos(k))    = dstcemlsgat(k)
      litrfallrow(ilmos(k),jlmos(k))    = litrfallgat(k)
      humiftrsrow(ilmos(k),jlmos(k))    = humiftrsgat(k)
      rmlrow(ilmos(k),jlmos(k))         = rmlgat(k)
      rmsrow(ilmos(k),jlmos(k))         = rmsgat(k)
      rmrrow(ilmos(k),jlmos(k))         = rmrgat(k)
      burnfracrow(ilmos(k),jlmos(k))    = burnfracgat(k)
      ltermrow(ilmos(k),jlmos(k))       = ltermgat(k)
      lucemcomrow(ilmos(k),jlmos(k))    = lucemcomgat(k)
      lucltrinrow(ilmos(k),jlmos(k))    = lucltringat(k)
      lucsocinrow(ilmos(k),jlmos(k))    = lucsocingat(k)
      dstcemls3row(ilmos(k),jlmos(k))   = dstcemls3gat(k)
      farerow(ilmos(k),jlmos(k))        = faregat(k)
      gavgscmsrow(ilmos(k),jlmos(k))    = gavgscmsgat(k)
      !          wetfracrow(ilmos(k))              = wetfracgat(k)
      !          slopefrac_row(ilmos(k))            = slopefrac_gat(k)
      ch4WetSpecrow(ilmos(k),jlmos(k))     = ch4WetSpecgat(k)
      wetfdynrow(ilmos(k),jlmos(k))     = wetfdyngat(k)
      ch4WetDynrow(ilmos(k),jlmos(k))     = ch4WetDyngat(k)
      ch4soillsrow(ilmos(k),jlmos(k))   = ch4soillsgat(k)
      wetfrac_presrow(ilmos(k),jlmos(k)) = wetfrac_presgat(k)

      twarmmrow(ilmos(k),jlmos(k)) = twarmmgat(k)
      tcoldmrow(ilmos(k),jlmos(k)) = tcoldmgat(k)
      gdd5row(ilmos(k),jlmos(k)) = gdd5gat(k)
      aridityrow(ilmos(k),jlmos(k)) = ariditygat(k)
      srplsmonrow(ilmos(k),jlmos(k)) = srplsmongat(k)
      defctmonrow(ilmos(k),jlmos(k)) = defctmongat(k)
      anndefctrow(ilmos(k),jlmos(k)) = anndefctgat(k)
      annsrplsrow(ilmos(k),jlmos(k)) = annsrplsgat(k)
      annpcprow(ilmos(k),jlmos(k)) = annpcpgat(k)
      dry_season_lengthrow(ilmos(k),jlmos(k)) = dry_season_lengthgat(k)

      tracerMossCMassrot(ilmos(k),jlmos(k)) = tracerMossCMassgat(k)
      tracermossLitrMassrot(ilmos(k),jlmos(k)) = tracermossLitrMassgat(k)
      !
    end do ! loop 100
    !
    do l = 1,icc ! loop 101
      do k = 1,nml
        smfuncvegrow(ilmos(k),jlmos(k),l) = smfuncveggat(k,l)
        mtermrow(ilmos(k),jlmos(k),l)     = mtermgat(k,l)
        btermrow(ilmos(k),jlmos(k),l)     = btermgat(k,l)
        ailcgrow(ilmos(k),jlmos(k),l)     = ailcggat(k,l)
        ailcgsrow(ilmos(k),jlmos(k),l)    = ailcgsgat(k,l)
        co2i1cgrow(ilmos(k),jlmos(k),l)   = co2i1cggat(k,l)
        co2i1csrow(ilmos(k),jlmos(k),l)   = co2i1csgat(k,l)
        co2i2cgrow(ilmos(k),jlmos(k),l)   = co2i2cggat(k,l)
        co2i2csrow(ilmos(k),jlmos(k),l)   = co2i2csgat(k,l)
        slairow(ilmos(k),jlmos(k),l)      = slaigat(k,l)
        anvegrow(ilmos(k),jlmos(k),l)     = anveggat(k,l)
        rmlvegrow(ilmos(k),jlmos(k),l)    = rmlveggat(k,l)
        pfcancmxrow(ilmos(k),jlmos(k),l)  = pfcancmxgat(k,l)
        nfcancmxrow(ilmos(k),jlmos(k),l)  = nfcancmxgat(k,l)
        fcancmxrow(ilmos(k),jlmos(k),l)   = fcancmxgat(k,l)
        stemmassrow(ilmos(k),jlmos(k),l)  = stemmassgat(k,l)
        rootmassrow(ilmos(k),jlmos(k),l)  = rootmassgat(k,l)
        pstemmassrow(ilmos(k),jlmos(k),l)  = pstemmassgat(k,l)
        pgleafmassrow(ilmos(k),jlmos(k),l)  = pgleafmassgat(k,l)
        gleafmasrow(ilmos(k),jlmos(k),l)  = gleafmasgat(k,l)
        pgleafmassrow(ilmos(k),jlmos(k),l)  = pgleafmassgat(k,l)
        bleafmasrow(ilmos(k),jlmos(k),l)  = bleafmasgat(k,l)
        ailcbrow(ilmos(k),jlmos(k),l)     = ailcbgat(k,l)
        flhrlossrow(ilmos(k),jlmos(k),l)  = flhrlossgat(k,l)
        pandaysrow(ilmos(k),jlmos(k),l)   = pandaysgat(k,l)
        lfstatusrow(ilmos(k),jlmos(k),l)  = lfstatusgat(k,l)
        grwtheffrow(ilmos(k),jlmos(k),l)  = grwtheffgat(k,l)
        lystmmasrow(ilmos(k),jlmos(k),l)  = lystmmasgat(k,l)
        lyrotmasrow(ilmos(k),jlmos(k),l)  = lyrotmasgat(k,l)
        tymaxlairow(ilmos(k),jlmos(k),l)  = tymaxlaigat(k,l)
        stmhrlosrow(ilmos(k),jlmos(k),l)  = stmhrlosgat(k,l)
        bmasvegrow(ilmos(k),jlmos(k),l)   = bmasveggat(k,l)
        rothrlosrow(ilmos(k),jlmos(k),l)  = rothrlosgat(k,l)
        veghghtrow(ilmos(k),jlmos(k),l)   = veghghtgat(k,l)
        rootdpthrow(ilmos(k),jlmos(k),l)  = rootdpthgat(k,l)
        tltrleafrow(ilmos(k),jlmos(k),l)  = tltrleafgat(k,l)
        tltrstemrow(ilmos(k),jlmos(k),l)  = tltrstemgat(k,l)
        tltrrootrow(ilmos(k),jlmos(k),l)  = tltrrootgat(k,l)
        leaflitrrow(ilmos(k),jlmos(k),l)  = leaflitrgat(k,l)
        roottemprow(ilmos(k),jlmos(k),l)  = roottempgat(k,l)
        afrleafrow(ilmos(k),jlmos(k),l)   = afrleafgat(k,l)
        afrstemrow(ilmos(k),jlmos(k),l)   = afrstemgat(k,l)
        afrrootrow(ilmos(k),jlmos(k),l)   = afrrootgat(k,l)
        wtstatusrow(ilmos(k),jlmos(k),l)  = wtstatusgat(k,l)
        ltstatusrow(ilmos(k),jlmos(k),l)  = ltstatusgat(k,l)
        nppvegrow(ilmos(k),jlmos(k),l)    = nppveggat(k,l)
        rmlvegaccrow(ilmos(k),jlmos(k),l) = rmlvegaccgat(k,l)
        rmsvegrow(ilmos(k),jlmos(k),l)    = rmsveggat(k,l)
        rmrvegrow(ilmos(k),jlmos(k),l)    = rmrveggat(k,l)
        rgvegrow(ilmos(k),jlmos(k),l)     = rgveggat(k,l)
        gppvegrow(ilmos(k),jlmos(k),l)    = gppveggat(k,l)
        vgbiomas_vegrow(ilmos(k),jlmos(k),l) = vgbiomas_veggat(k,l)
        autoresvegrow(ilmos(k),jlmos(k),l) = autoresveggat(k,l)
        pftexistrow(ilmos(k),jlmos(k),l) = pftexistgat(k,l)
        ancsvegrow(ilmos(k),jlmos(k),l) = ancsveggat(k,l)
        ancgvegrow(ilmos(k),jlmos(k),l) = ancgveggat(k,l)
        rmlcsvegrow(ilmos(k),jlmos(k),l) = rmlcsveggat(k,l)
        rmlcgvegrow(ilmos(k),jlmos(k),l) = rmlcgveggat(k,l)
        litrfallvegrow(ilmos(k),jlmos(k),l) = litrfallveggat(k,l)

        !         fire variables
        emit_co2row(ilmos(k),jlmos(k),l)    = emit_co2gat(k,l)
        emit_corow(ilmos(k),jlmos(k),l)     = emit_cogat(k,l)
        emit_ch4row(ilmos(k),jlmos(k),l)    = emit_ch4gat(k,l)
        emit_nmhcrow(ilmos(k),jlmos(k),l)   = emit_nmhcgat(k,l)
        emit_h2row(ilmos(k),jlmos(k),l)     = emit_h2gat(k,l)
        emit_noxrow(ilmos(k),jlmos(k),l)    = emit_noxgat(k,l)
        emit_n2orow(ilmos(k),jlmos(k),l)    = emit_n2ogat(k,l)
        emit_pm25row(ilmos(k),jlmos(k),l)   = emit_pm25gat(k,l)
        emit_tpmrow(ilmos(k),jlmos(k),l)    = emit_tpmgat(k,l)
        emit_tcrow(ilmos(k),jlmos(k),l)     = emit_tcgat(k,l)
        emit_ocrow(ilmos(k),jlmos(k),l)     = emit_ocgat(k,l)
        emit_bcrow(ilmos(k),jlmos(k),l)     = emit_bcgat(k,l)
        burnvegfrow(ilmos(k),jlmos(k),l)    = burnvegfgat(k,l)

        tracergLeafMassrot(ilmos(k),jlmos(k),l) = tracergLeafMassgat(k,l)
        tracerbLeafMassrot(ilmos(k),jlmos(k),l) = tracerbLeafMassgat(k,l)
        tracerStemMassrot(ilmos(k),jlmos(k),l) = tracerStemMassgat(k,l)
        tracerRootMassrot(ilmos(k),jlmos(k),l) = tracerRootMassgat(k,l)
      end do
    end do ! loop 101
    !
    do l = 1,iccp1 ! loop 102
      do k = 1,nml
        hetroresvegrow(ilmos(k),jlmos(k),l) = hetroresveggat(k,l)
        nepvegrow(ilmos(k),jlmos(k),l)    = nepveggat(k,l)
        nbpvegrow(ilmos(k),jlmos(k),l)    = nbpveggat(k,l)
      end do
    end do ! loop 102

    do l = 1,iccp2 ! loop 103
      do k = 1,nml
        ! COMBAK PERLAY
        litrmassrow(ilmos(k),jlmos(k),l) = litrmassgat(k,l)
        soilcmasrow(ilmos(k),jlmos(k),l) = soilcmasgat(k,l)
        litresvegrow(ilmos(k),jlmos(k),l) = litresveggat(k,l)
        soilcresvegrow(ilmos(k),jlmos(k),l) = soilcresveggat(k,l)
        humiftrsvegrow(ilmos(k),jlmos(k),l) = humiftrsveggat(k,l)
        do m = 1,ignd
          ! litrmassrow(ilmos(k),jlmos(k),l,m) = litrmassgat(k,l,m)
          ! soilcmasrow(ilmos(k),jlmos(k),l,m) = soilcmasgat(k,l,m)
          ! litresvegrow(ilmos(k),jlmos(k),l,m) = litresveggat(k,l,m)
          ! soilcresvegrow(ilmos(k),jlmos(k),l,m)=soilcresveggat(k,l,m)
          ! humiftrsvegrow(ilmos(k),jlmos(k),l,m) = humiftrsveggat(k,l,m)
          ! COMBAK PERLAY
          tracerSoilCMassrot(ilmos(k),jlmos(k),l,m) = tracerSoilCMassgat(k,l,m)
          tracerLitrMassrot(ilmos(k),jlmos(k),l,m) = tracerLitrMassgat(k,l,m)
        end do
      end do
    end do ! loop 103

    do l = 1,2 ! loop 106
      do k = 1,nml
        colddaysrow(ilmos(k),jlmos(k),l) = colddaysgat(k,l)
      end do
    end do ! loop 106
    !
    do l = 1,ican ! loop 201
      do k = 1,nml
        ailcrow(ilmos(k),jlmos(k),l)    = ailcgat(k,l)
        zolncrow(ilmos(k),jlmos(k),l)   = zolncgat(k,l)
        cmasvegcrow(ilmos(k),jlmos(k),l) = cmasvegcgat(k,l)
        paicrow(ilmos(k),jlmos(k),l)    = paicgat(k,l)
        slaicrow(ilmos(k),jlmos(k),l)   = slaicgat(k,l)
        alvsctmrow(ilmos(k),jlmos(k),l) = alvsctmgat(k,l)
        alirctmrow(ilmos(k),jlmos(k),l) = alirctmgat(k,l)
      end do
    end do ! loop 201
    !
    do l = 1,ignd ! loop 250
      do k = 1,nml
        sandrow(ilmos(k),jlmos(k),l)     = sandgat(k,l)
        clayrow(ilmos(k),jlmos(k),l)     = claygat(k,l)
        orgmrow(ilmos(k),jlmos(k),l)     = orgmgat(k,l)
        tbaraccrow_m(ilmos(k),jlmos(k),l) = tbaraccgat_m(k,l)
        ! thlqaccrow_m(ilmos(k),jlmos(k),l)= thlqaccgat_m(k,l)
        ! thicaccrow_m(ilmos(k),jlmos(k),l)= thicaccgat_m(k,l)
      end do
    end do ! loop 250
    !
    do l = 1,icc ! loop 280
      do m = 1,ignd
        do k = 1,nml
          rmatctemrow(ilmos(k),jlmos(k),l,m) = rmatctemgat(k,l,m)
        end do
      end do
    end do ! loop 280
    !
    do l = 1,ican ! loop 290
      do m = 1,ignd
        do k = 1,nml
          rmatcrow(ilmos(k),jlmos(k),l,m) = rmatcgat(k,l,m)
        end do
      end do
    end do ! loop 290

    !     this class variable is scattered here, but it gathered in classg,
    !     not in ctemg2. jm jan 82013.
    do l = 1,icp1 ! loop 300
      do k = 1,nml
        fcanrow(ilmos(k),jlmos(k),l)     = fcangat(k,l)
      end do
    end do ! loop 300

    !    scatter peatland variables----------------------------------------\
    do k = 1,nml ! loop 400
      anmossrow(ilmos(k),jlmos(k)) = anmossgat(k)
      rmlmossrow(ilmos(k),jlmos(k)) = rmlmossgat(k)
      gppmossrow(ilmos(k),jlmos(k)) = gppmossgat(k)
      armossrow(ilmos(k),jlmos(k)) = armossgat(k)
      nppmossrow(ilmos(k),jlmos(k)) = nppmossgat(k)
      peatdeprow (ilmos(k),jlmos(k))  = peatdepgat(k)
      litrmsmossrow(ilmos(k),jlmos(k)) = litrmsmossgat(k)
      Cmossmasrow(ilmos(k),jlmos(k)) = Cmossmasgat(k)
      dmossrow(ilmos(k),jlmos(k)) = dmossgat(k)
      ipeatlandrow(ilmos(k),jlmos(k)) = ipeatlandgat(k)
      pddrow(ilmos(k),jlmos(k)) = pddgat(k)
    end do ! loop 400

    return

  end subroutine ctems2
  !! @}
  ! ------------------------------------------------------------------------------------
  !> \ingroup ctemgatherscatter_ctems2
  !! @{
  !> Performs subsequent 'gather' operation on CTEM variables for consistency
  !! with physics variables gather operations.
  !!
  !! ctemg2 takes variables in the 'row' format (nlat, nmos, ...)
  !! and converts them to the 'gat' format (ilg, ...). At present
  !! ctemg2 is bloated with many variables that do not require
  !! gathering. This subroutine should ideally be only used for
  !! state variables that are updated from external files as
  !! the run progresses. Since the model calculations operate
  !! on the 'gat' form, any other variables need not be gathered
  !! as they will already be in the correct format from the previous
  !! model timestep.
  !!
  !! @author R. Li, Y. Wu, E. Chan, J. Melton
  !!
  subroutine ctemg2 (fcancmxgat, rmatcgat, zolncgat, &
                     paicgat, &
                     ailcgat, ailcggat, cmasvegcgat, slaicgat, &
                     ailcgsgat, fcancsgat, fcancgat, rmatctemgat, &
                     co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                     co2i2csgat, xdiffusgat, slaigat, cfluxcggat, &
                     cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                     rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                     sandgat, claygat, orgmgat, &
                     anveggat, rmlveggat, tbaraccgat_m, prbfrhucgat, &
                     extnprobgat, pfcancmxgat, nfcancmxgat, &
                     stemmassgat, rootmassgat, litrmassgat, gleafmasgat, &
                     bleafmasgat, soilcmasgat, ailcbgat, flhrlossgat, &
                     pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                     lyrotmasgat, tymaxlaigat, vgbiomasgat, gavgltmsgat, &
                     stmhrlosgat, bmasveggat, colddaysgat, rothrlosgat, &
                     alvsctmgat, alirctmgat, gavglaigat, nppgat, &
                     nepgat, hetroresgat, autoresgat, soilrespgat, &
                     rmgat, rggat, nbpgat, litresgat, &
                     socresgat, gppgat, dstcemlsgat, litrfallgat, &
                     humiftrsgat, veghghtgat, rootdpthgat, rmlgat, &
                     rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                     tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                     afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                     burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                     lucsocingat, dstcemls3gat, popdingat, &
                     faregat, gavgscmsgat, rmlvegaccgat, pftexistgat, &
                     rmsveggat, rmrveggat, rgveggat, vgbiomas_veggat, &
                     gppveggat, nepveggat, &
                     emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                     emit_h2gat, emit_noxgat, emit_n2ogat, emit_pm25gat, &
                     emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                     btermgat, ltermgat, mtermgat, daylgat, dayl_maxgat, &
                     nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                     soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                     ch4WetSpecgat, slopefracgat, &
                     wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                     twarmmgat, tcoldmgat, gdd5gat, &
                     ariditygat, srplsmongat, defctmongat, anndefctgat, &
                     annsrplsgat, annpcpgat, dry_season_lengthgat, &
                     anmossgat, rmlmossgat, gppmossgat, armossgat, nppmossgat, &
                     litrmsmossgat, peatdepgat, Cmossmasgat, dmossgat, &
                     ipeatlandgat, pddgat, tracerCO2gat, &
                     ilmos, jlmos, iwmos, jwmos, &
                     nml, fcancmxrow, rmatcrow, zolncrow, paicrow, &
                     ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                     ailcgsrow, fcancsrow, fcancrow, rmatctemrow, &
                     co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                     co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                     cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                     rmlcgvegrow, canresrow, sdeprow, ch4concrow, &
                     sandrow, clayrow, orgmrow, &
                     anvegrow, rmlvegrow, tbaraccrow_m, prbfrhucrow, &
                     extnprobrow, pfcancmxrow, nfcancmxrow, &
                     stemmassrow, rootmassrow, litrmassrow, gleafmasrow, &
                     bleafmasrow, soilcmasrow, ailcbrow, flhrlossrow, &
                     pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                     lyrotmasrow, tymaxlairow, vgbiomasrow, gavgltmsrow, &
                     stmhrlosrow, bmasvegrow, colddaysrow, rothrlosrow, &
                     alvsctmrow, alirctmrow, gavglairow, npprow, &
                     neprow, hetroresrow, autoresrow, soilresprow, &
                     rmrow, rgrow, nbprow, litresrow, &
                     socresrow, gpprow, dstcemlsrow, litrfallrow, &
                     humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                     rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                     tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                     afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                     burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                     lucsocinrow, dstcemls3row, popdinrow, &
                     farerow, gavgscmsrow, rmlvegaccrow, pftexistrow, &
                     rmsvegrow, rmrvegrow, rgvegrow, vgbiomas_vegrow, &
                     gppvegrow, nepvegrow, &
                     emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                     emit_h2row, emit_noxrow, emit_n2orow, emit_pm25row, &
                     emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                     btermrow, ltermrow, mtermrow, daylrow, dayl_maxrow, &
                     nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                     soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                     ch4WetSpecrow, slopefracrow, &
                     wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                     twarmmrow, tcoldmrow, gdd5row, &
                     aridityrow, srplsmonrow, defctmonrow, anndefctrow, &
                     annsrplsrow, annpcprow, dry_season_lengthrow, &
                     anmossrow, rmlmossrow, gppmossrow, armossrow, nppmossrow, &
                     litrmsmossrow, peatdeprow, Cmossmasrow, dmossrow, &
                     ipeatlandrow, pddrow, tracerCO2rot)
    !               thlqaccrow_m, thicaccrow_m, ipeatlandrow, pddrow)


    !     Dec 232016     Remove thlqaccXXX_m/thicaccXXX_m
    !     Ed Chan
    !     March 192015   Gathering of peatland variables
    !     Yuanqiao Wu
    !
    !     July 122013    Bring in the ctem params use statement
    !     J. Melton

    !     July 282009    Gather operation on CTEM variables.
    !     Rong Li
    !
    use classicParams,      only : nlat, nmos, ilg, ignd, ican, icp1, &
                                   icc, iccp2, iccp1

    implicit none
    !
    !     * integer :: constants.
    !
    integer :: k, l, m, n, j
    integer, intent(in) :: nml
    !
    !     * gather-scatter index arrays.
    !
    integer, intent(in) :: ilmos(ilg), jlmos(ilg), iwmos(ilg), jwmos(ilg)
    !
    !
    real, intent(out) :: fcancmxgat(ilg,icc), rmatcgat(ilg,ican,ignd), &
                         zolncgat(ilg,ican), paicgat(ilg,ican), &
                         ailcgat(ilg,ican), ailcggat(ilg,icc), &
                         cmasvegcgat(ilg,ican), slaicgat(ilg,ican), &
                         ailcgsgat(ilg,icc), fcancsgat(ilg,icc), &
                         fcancgat(ilg,icc), rmatctemgat(ilg,icc,ignd), &
                         co2concgat(ilg), co2i1cggat(ilg,icc), &
                         co2i1csgat(ilg,icc), co2i2cggat(ilg,icc), &
                         co2i2csgat(ilg,icc), xdiffusgat(ilg), &
                         slaigat(ilg,icc), cfluxcggat(ilg), &
                         cfluxcsgat(ilg), ancsveggat(ilg,icc), &
                         ancgveggat(ilg,icc), rmlcsveggat(ilg,icc), &
                         rmlcgveggat(ilg,icc), canresgat(ilg), &
                         sdepgat(ilg), ch4concgat(ilg)
    !
    real, intent(out) :: sandgat(ilg,ignd), claygat(ilg,ignd), &
                         orgmgat(ilg,ignd)
    !
    real, intent(out) :: anveggat(ilg,icc), rmlveggat(ilg,icc)
    !
    real, intent(out) :: prbfrhucgat(ilg), extnprobgat(ilg), &
                         daylgat(ilg), dayl_maxgat(ilg), &
                         tbaraccgat_m(ilg,ignd), &
                         pfcancmxgat(ilg,icc), nfcancmxgat(ilg,icc), &
                         stemmassgat(ilg,icc), rootmassgat(ilg,icc), &
                         pstemmassgat(ilg,icc), pgleafmassgat(ilg,icc), &
                         gleafmasgat(ilg,icc), bleafmasgat(ilg,icc), &
                         ailcbgat(ilg,icc), flhrlossgat(ilg,icc), &
                         ! COMBAK PERLAY
                         soilcmasgat(ilg,iccp2), litrmassgat(ilg,iccp2)
    ! soilcmasgat(ilg,iccp2,ignd), litrmassgat(ilg,iccp2,ignd)
    ! COMBAK PERLAY
    !
    integer, intent(out) :: pandaysgat(ilg,icc), lfstatusgat(ilg,icc), &
                            colddaysgat(ilg,2)

    logical, intent(out) :: pftexistgat(ilg,icc)
    !
    real, intent(out) :: grwtheffgat(ilg,icc), lystmmasgat(ilg,icc), &
                         lyrotmasgat(ilg,icc), tymaxlaigat(ilg,icc), &
                         vgbiomasgat(ilg), gavgltmsgat(ilg), &
                         stmhrlosgat(ilg,icc), bmasveggat(ilg,icc), &
                         rothrlosgat(ilg,icc), &
                         alvsctmgat(ilg,ican), alirctmgat(ilg,ican), &
                         gavglaigat(ilg)
    !
    real, intent(out) :: nppgat(ilg), nepgat(ilg), &
                         hetroresgat(ilg), autoresgat(ilg), &
                         soilrespgat(ilg), rmgat(ilg), rggat(ilg), &
                         nbpgat(ilg), litresgat(ilg), &
                         socresgat(ilg), gppgat(ilg), &
                         dstcemlsgat(ilg), litrfallgat(ilg), &
                         humiftrsgat(ilg), veghghtgat(ilg,icc), &
                         rootdpthgat(ilg,icc), rmlgat(ilg), &
                         rmsgat(ilg), rmrgat(ilg), &
                         tltrleafgat(ilg,icc), tltrstemgat(ilg,icc), &
                         tltrrootgat(ilg,icc), leaflitrgat(ilg,icc), &
                         roottempgat(ilg,icc), afrleafgat(ilg,icc), &
                         afrstemgat(ilg,icc), afrrootgat(ilg,icc), &
                         wtstatusgat(ilg,icc), ltstatusgat(ilg,icc), &
                         burnfracgat(ilg), smfuncveggat(ilg,icc), &
                         lucemcomgat(ilg), lucltringat(ilg), &
                         lucsocingat(ilg), &
                         dstcemls3gat(ilg), &
                         tracerCO2gat(ilg)
    !
    !      fire emission variables
    real, intent(out) :: emit_co2gat(ilg,icc), emit_cogat(ilg,icc), &
                         emit_ch4gat(ilg,icc), emit_nmhcgat(ilg,icc), &
                         emit_h2gat(ilg,icc), emit_noxgat(ilg,icc), &
                         emit_n2ogat(ilg,icc), emit_pm25gat(ilg,icc), &
                         emit_tpmgat(ilg,icc), emit_tcgat(ilg,icc), &
                         emit_ocgat(ilg,icc), emit_bcgat(ilg,icc), &
                         burnvegfgat(ilg,icc), btermgat(ilg,icc), &
                         ltermgat(ilg), mtermgat(ilg,icc), &
                         popdingat(ilg)
    !
    real, intent(out) :: faregat(ilg)
    real, intent(out) :: gavgscmsgat(ilg)
    !
    real, intent(out) :: rmlvegaccgat(ilg,icc), rmsveggat(ilg,icc), &
                         rmrveggat(ilg,icc), rgveggat(ilg,icc)
    !
    real, intent(out) :: vgbiomas_veggat(ilg,icc)
    !
    real, intent(out) :: gppveggat(ilg,icc), nepveggat(ilg,iccp1), &
                         nbpveggat(ilg,iccp1), hetroresveggat(ilg,iccp1), &
                         autoresveggat(ilg,icc), &
                         ! COMBAK PERLAY
                         litresveggat(ilg,iccp2), &
                         soilcresveggat(ilg,iccp2)
    !        litresveggat(ilg,iccp2,ignd), &
    ! &      soilcresveggat(ilg,iccp2,ignd)
    ! COMBAK PERLAY
    !
    real, intent(in) :: fcancmxrow(nlat,nmos,icc), rmatcrow(nlat,nmos,ican,ignd), &
                        zolncrow(nlat,nmos,ican), paicrow(nlat,nmos,ican), &
                        ailcrow(nlat,nmos,ican), ailcgrow(nlat,nmos,icc), &
                        cmasvegcrow(nlat,nmos,ican), slaicrow(nlat,nmos,ican), &
                        ailcgsrow(nlat,nmos,icc), fcancsrow(nlat,nmos,icc), &
                        fcancrow(nlat,nmos,icc), rmatctemrow(nlat,nmos,icc,ignd), &
                        co2concrow(nlat,nmos), co2i1cgrow(nlat,nmos,icc), &
                        co2i1csrow(nlat,nmos,icc), co2i2cgrow(nlat,nmos,icc), &
                        co2i2csrow(nlat,nmos,icc), xdiffus(nlat), &
                        slairow(nlat,nmos,icc), cfluxcgrow(nlat,nmos), &
                        cfluxcsrow(nlat,nmos), ancsvegrow(nlat,nmos,icc), &
                        ancgvegrow(nlat,nmos,icc), rmlcsvegrow(nlat,nmos,icc), &
                        rmlcgvegrow(nlat,nmos,icc), canresrow(nlat,nmos), &
                        sdeprow(nlat,nmos), ch4concrow(nlat,nmos)
    !
    real, intent(in) :: sandrow(nlat,nmos,ignd), clayrow(nlat,nmos,ignd), &
                        orgmrow(nlat,nmos,ignd)
    !
    real, intent(in) :: anvegrow(nlat,nmos,icc), rmlvegrow(nlat,nmos,icc)
    !
    real, intent(in) :: prbfrhucrow(nlat,nmos), extnprobrow(nlat,nmos), &
                        daylrow(nlat), dayl_maxrow(nlat), &
                        tbaraccrow_m(nlat,nmos,ignd), &
                        pfcancmxrow(nlat,nmos,icc), nfcancmxrow(nlat,nmos,icc), &
                        stemmassrow(nlat,nmos,icc), rootmassrow(nlat,nmos,icc), &
                        pstemmassrow(nlat,nmos,icc), pgleafmassrow(nlat,nmos,icc), &
                        gleafmasrow(nlat,nmos,icc), bleafmasrow(nlat,nmos,icc), &
                        ailcbrow(nlat,nmos,icc), flhrlossrow(nlat,nmos,icc), &
                        ! COMBAK PERLAY
                        soilcmasrow(nlat,nmos,iccp2), &
                        litrmassrow(nlat,nmos,iccp2)
    ! &      soilcmasrow(nlat,nmos,iccp2,ignd), &
    ! &      litrmassrow(nlat,nmos,iccp2,ignd)
    ! COMBAK PERLAY
    !
    integer, intent(in) :: pandaysrow(nlat,nmos,icc), lfstatusrow(nlat,nmos,icc), &
                           colddaysrow(nlat,nmos,2)

    logical, intent(in) :: pftexistrow(nlat,nmos,icc)
    !
    real, intent(in) :: grwtheffrow(nlat,nmos,icc), lystmmasrow(nlat,nmos,icc), &
                        lyrotmasrow(nlat,nmos,icc), tymaxlairow(nlat,nmos,icc), &
                        vgbiomasrow(nlat,nmos), gavgltmsrow(nlat,nmos), &
                        stmhrlosrow(nlat,nmos,icc), bmasvegrow(nlat,nmos,icc), &
                        rothrlosrow(nlat,nmos,icc), &
                        alvsctmrow(nlat,nmos,ican), alirctmrow(nlat,nmos,ican), &
                        gavglairow(nlat,nmos)
    !
    real, intent(in) :: npprow(nlat,nmos), neprow(nlat,nmos), &
                        hetroresrow(nlat,nmos), autoresrow(nlat,nmos), &
                        soilresprow(nlat,nmos), rmrow(nlat,nmos), &
                        rgrow(nlat,nmos), &
                        nbprow(nlat,nmos), litresrow(nlat,nmos), &
                        socresrow(nlat,nmos), gpprow(nlat,nmos), &
                        dstcemlsrow(nlat,nmos), litrfallrow(nlat,nmos), &
                        humiftrsrow(nlat,nmos), veghghtrow(nlat,nmos,icc), &
                        rootdpthrow(nlat,nmos,icc), rmlrow(nlat,nmos), &
                        rmsrow(nlat,nmos), rmrrow(nlat,nmos), &
                        tltrleafrow(nlat,nmos,icc), tltrstemrow(nlat,nmos,icc), &
                        tltrrootrow(nlat,nmos,icc), leaflitrrow(nlat,nmos,icc), &
                        roottemprow(nlat,nmos,icc), afrleafrow(nlat,nmos,icc), &
                        afrstemrow(nlat,nmos,icc), afrrootrow(nlat,nmos,icc), &
                        wtstatusrow(nlat,nmos,icc), ltstatusrow(nlat,nmos,icc), &
                        burnfracrow(nlat,nmos), smfuncvegrow(nlat,nmos,icc), &
                        lucemcomrow(nlat,nmos), lucltrinrow(nlat,nmos), &
                        lucsocinrow(nlat,nmos), &
                        dstcemls3row(nlat,nmos), tracerCO2rot(nlat,nmos)
    !
    !     fire variables
    real, intent(in) :: emit_co2row(nlat,nmos,icc), emit_corow(nlat,nmos,icc), &
                        emit_ch4row(nlat,nmos,icc), emit_nmhcrow(nlat,nmos,icc), &
                        emit_h2row(nlat,nmos,icc), emit_noxrow(nlat,nmos,icc), &
                        emit_n2orow(nlat,nmos,icc), emit_pm25row(nlat,nmos,icc), &
                        emit_tpmrow(nlat,nmos,icc), emit_tcrow(nlat,nmos,icc), &
                        emit_ocrow(nlat,nmos,icc), emit_bcrow(nlat,nmos,icc), &
                        burnvegfrow(nlat,nmos,icc), btermrow(nlat,nmos,icc), &
                        ltermrow(nlat,nmos), mtermrow(nlat,nmos,icc), &
                        popdinrow(nlat,nmos)
    !
    real, intent(in) :: farerow(nlat,nmos)
    !
    real, intent(in) :: gavgscmsrow(nlat,nmos)
    !
    real, intent(in) :: rmlvegaccrow(nlat,nmos,icc), rmsvegrow(nlat,nmos,icc), &
                        rmrvegrow(nlat,nmos,icc), rgvegrow(nlat,nmos,icc)
    !
    real, intent(in) :: vgbiomas_vegrow(nlat,nmos,icc)
    !
    real, intent(in) :: gppvegrow(nlat,nmos,icc), nepvegrow(nlat,nmos,iccp1), &
                        nbpvegrow(nlat,nmos,iccp1), hetroresvegrow(nlat,nmos,iccp1), &
                        autoresvegrow(nlat,nmos,icc), &
                        ! COMBAK PERLAY
                        litresvegrow(nlat,nmos,iccp2), &
                        soilcresvegrow(nlat,nmos,iccp2)
    ! &      litresvegrow(nlat,nmos,iccp2,ignd), &
    ! &      soilcresvegrow(nlat,nmos,iccp2,ignd)
    ! COMBAK PERLAY

    !   Methane related variables
    real, intent(in)  :: slopefracrow(nlat,nmos,8), ch4WetSpecrow(nlat,nmos), wetfdynrow(nlat,nmos), &
                         ch4WetDynrow(nlat,nmos), ch4soillsrow(nlat,nmos)

    real, intent(out) :: slopefracgat(ilg,8), ch4WetSpecgat(ilg), wetfdyngat(ilg), &
                         ch4WetDyngat(ilg), ch4soillsgat(ilg)

    real, intent(in)  :: twarmmrow(nlat,nmos), tcoldmrow(nlat,nmos), gdd5row(nlat,nmos), &
                         aridityrow(nlat,nmos), srplsmonrow(nlat,nmos), defctmonrow(nlat,nmos), &
                         anndefctrow(nlat,nmos), annsrplsrow(nlat,nmos), annpcprow(nlat,nmos), &
                         dry_season_lengthrow(nlat,nmos)

    real, intent(out) :: twarmmgat(ilg), tcoldmgat(ilg), gdd5gat(ilg), &
                         ariditygat(ilg), srplsmongat(ilg), defctmongat(ilg), &
                         anndefctgat(ilg), annsrplsgat(ilg), annpcpgat(ilg), &
                         dry_season_lengthgat(ilg)

    !    --peatland variables
    real, intent(in)  :: anmossrow(nlat,nmos), rmlmossrow(nlat,nmos), gppmossrow(nlat,nmos), &
                         armossrow(nlat,nmos), nppmossrow(nlat,nmos), peatdeprow(nlat,nmos), &
                         litrmsmossrow(nlat,nmos), Cmossmasrow(nlat,nmos), dmossrow(nlat,nmos)

    real, intent(out) :: anmossgat(ilg), rmlmossgat(ilg), gppmossgat(ilg), &
                         armossgat(ilg), nppmossgat(ilg), peatdepgat(ilg), &
                         litrmsmossgat(ilg), Cmossmasgat(ilg), dmossgat(ilg), &
                         !    9                   thlqaccrow_m(nlat,nmos,ignd), thlqaccgat_m(ilg,ignd), &
                         !    1                   thicaccrow_m(nlat,nmos,ignd), thicaccgat_m(ilg,ignd), &
                         pddrow(nlat,nmos), pddgat(ilg)

    integer, intent(out) :: ipeatlandgat(ilg), ipeatlandrow(nlat,nmos)


    !----------------------------------------------------------------------
    do k = 1,nml ! loop 100
      sdepgat(k)      = sdeprow(ilmos(k),jlmos(k))
      co2concgat(k)   = co2concrow(ilmos(k),jlmos(k))
      ch4concgat(k)   = ch4concrow(ilmos(k),jlmos(k))
      cfluxcggat(k)   = cfluxcgrow(ilmos(k),jlmos(k))
      cfluxcsgat(k)   = cfluxcsrow(ilmos(k),jlmos(k))
      canresgat(k)    = canresrow(ilmos(k),jlmos(k))
      xdiffusgat(k)   = xdiffus(ilmos(k))
      prbfrhucgat(k)  = prbfrhucrow(ilmos(k),jlmos(k))
      extnprobgat(k)  = extnprobrow(ilmos(k),jlmos(k))
      daylgat(k)      = daylrow(ilmos(k))
      dayl_maxgat(k)  = dayl_maxrow(ilmos(k))
      vgbiomasgat(k)  = vgbiomasrow(ilmos(k),jlmos(k))
      gavgltmsgat(k)  = gavgltmsrow(ilmos(k),jlmos(k))
      gavglaigat(k)   = gavglairow(ilmos(k),jlmos(k))
      nppgat(k)       = npprow(ilmos(k),jlmos(k))
      nepgat(k)       = neprow(ilmos(k),jlmos(k))
      hetroresgat(k)  = hetroresrow(ilmos(k),jlmos(k))
      autoresgat(k)   = autoresrow(ilmos(k),jlmos(k))
      soilrespgat(k)  = soilresprow(ilmos(k),jlmos(k))
      rmgat(k)        = rmrow(ilmos(k),jlmos(k))
      rggat(k)        = rgrow(ilmos(k),jlmos(k))
      nbpgat(k)       = nbprow(ilmos(k),jlmos(k))
      litresgat(k)    = litresrow(ilmos(k),jlmos(k))
      socresgat(k)    = socresrow(ilmos(k),jlmos(k))
      gppgat(k)       = gpprow(ilmos(k),jlmos(k))
      dstcemlsgat(k)  = dstcemlsrow(ilmos(k),jlmos(k))
      litrfallgat(k)  = litrfallrow(ilmos(k),jlmos(k))
      humiftrsgat(k)  = humiftrsrow(ilmos(k),jlmos(k))
      rmlgat(k)       = rmlrow(ilmos(k),jlmos(k))
      rmsgat(k)       = rmsrow(ilmos(k),jlmos(k))
      rmrgat(k)       = rmrrow(ilmos(k),jlmos(k))
      burnfracgat(k)  = burnfracrow(ilmos(k),jlmos(k))
      lucemcomgat(k)  = lucemcomrow(ilmos(k),jlmos(k))
      lucltringat(k)  = lucltrinrow(ilmos(k),jlmos(k))
      lucsocingat(k)  = lucsocinrow(ilmos(k),jlmos(k))
      dstcemls3gat(k) = dstcemls3row(ilmos(k),jlmos(k))
      tracerCO2gat(k) = tracerCO2rot(ilmos(k),jlmos(k))
      faregat(k)      = farerow(ilmos(k),jlmos(k))
      gavgscmsgat(k)  = gavgscmsrow(ilmos(k),jlmos(k))
      do n = 1,8
        slopefracgat(k,n) = slopefracrow(ilmos(k),jlmos(k),n)
      end do
      ch4WetSpecgat(k)   = ch4WetSpecrow(ilmos(k),jlmos(k))
      wetfdyngat(k)   = wetfdynrow(ilmos(k),jlmos(k))
      ch4WetDyngat(k)   = ch4WetDynrow(ilmos(k),jlmos(k))
      ch4soillsgat(k) = ch4soillsrow(ilmos(k),jlmos(k))

      twarmmgat(k)    = twarmmrow(ilmos(k),jlmos(k))
      tcoldmgat(k)    = tcoldmrow(ilmos(k),jlmos(k))
      gdd5gat(k)      = gdd5row(ilmos(k),jlmos(k))
      ariditygat(k)   = aridityrow(ilmos(k),jlmos(k))
      srplsmongat(k)  = srplsmonrow(ilmos(k),jlmos(k))
      defctmongat(k)  = defctmonrow(ilmos(k),jlmos(k))
      anndefctgat(k)  = anndefctrow(ilmos(k),jlmos(k))
      annsrplsgat(k)  = annsrplsrow(ilmos(k),jlmos(k))
      annpcpgat(k)    = annpcprow(ilmos(k),jlmos(k))
      popdingat(k)    = popdinrow(ilmos(k),jlmos(k))
      dry_season_lengthgat(k) = &
                                dry_season_lengthrow(ilmos(k),jlmos(k))
      !
    end do ! loop 100
    !
    do l = 1,icc ! loop 101
      do k = 1,nml
        btermgat(k,l)    = btermrow(ilmos(k),jlmos(k),l)
        ltermgat(k)      = ltermrow(ilmos(k),jlmos(k))  ! not per pft but keeping in with the other term vars.
        mtermgat(k,l)    = mtermrow(ilmos(k),jlmos(k),l)
        smfuncveggat(k,l) = smfuncvegrow(ilmos(k),jlmos(k),l)
        fcancmxgat(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
        !          ailcggat(k,l)    = ailcgrow(ilmos(k),jlmos(k),l)
        ailcgsgat(k,l)   = ailcgsrow(ilmos(k),jlmos(k),l)
        fcancsgat(k,l)   = fcancsrow(ilmos(k),jlmos(k),l)
        fcancgat(k,l)    = fcancrow(ilmos(k),jlmos(k),l)
        co2i1cggat(k,l)  = co2i1cgrow(ilmos(k),jlmos(k),l)
        co2i1csgat(k,l)  = co2i1csrow(ilmos(k),jlmos(k),l)
        co2i2cggat(k,l)  = co2i2cgrow(ilmos(k),jlmos(k),l)
        co2i2csgat(k,l)  = co2i2csrow(ilmos(k),jlmos(k),l)
        !         slaigat(k,l)     = slairow(ilmos(k),jlmos(k),l)
        ancsveggat(k,l)  = ancsvegrow(ilmos(k),jlmos(k),l)
        ancgveggat(k,l)  = ancgvegrow(ilmos(k),jlmos(k),l)
        rmlcsveggat(k,l) = rmlcsvegrow(ilmos(k),jlmos(k),l)
        rmlcgveggat(k,l) = rmlcgvegrow(ilmos(k),jlmos(k),l)
        anveggat(k,l)    = anvegrow(ilmos(k),jlmos(k),l)
        rmlveggat(k,l)   = rmlvegrow(ilmos(k),jlmos(k),l)
        pfcancmxgat(k,l) = pfcancmxrow(ilmos(k),jlmos(k),l)
        nfcancmxgat(k,l) = nfcancmxrow(ilmos(k),jlmos(k),l)
        stemmassgat(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
        rootmassgat(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
        pstemmassgat(k,l) = pstemmassrow(ilmos(k),jlmos(k),l)
        pgleafmassgat(k,l) = pgleafmassrow(ilmos(k),jlmos(k),l)
        gleafmasgat(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
        bleafmasgat(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
        !          ailcbgat(k,l)    = ailcbrow(ilmos(k),jlmos(k),l)
        flhrlossgat(k,l) = flhrlossrow(ilmos(k),jlmos(k),l)
        pandaysgat(k,l)  = pandaysrow(ilmos(k),jlmos(k),l)
        lfstatusgat(k,l) = lfstatusrow(ilmos(k),jlmos(k),l)
        grwtheffgat(k,l) = grwtheffrow(ilmos(k),jlmos(k),l)
        lystmmasgat(k,l) = lystmmasrow(ilmos(k),jlmos(k),l)
        lyrotmasgat(k,l) = lyrotmasrow(ilmos(k),jlmos(k),l)
        tymaxlaigat(k,l) = tymaxlairow(ilmos(k),jlmos(k),l)
        stmhrlosgat(k,l) = stmhrlosrow(ilmos(k),jlmos(k),l)
        !         bmasveggat(k,l)  = bmasvegrow(ilmos(k),jlmos(k),l)
        rothrlosgat(k,l) = rothrlosrow(ilmos(k),jlmos(k),l)
        !         veghghtgat(k,l)  = veghghtrow(ilmos(k),jlmos(k),l)
        !         rootdpthgat(k,l) = rootdpthrow(ilmos(k),jlmos(k),l)
        tltrleafgat(k,l) = tltrleafrow(ilmos(k),jlmos(k),l)
        tltrstemgat(k,l) = tltrstemrow(ilmos(k),jlmos(k),l)
        tltrrootgat(k,l) = tltrrootrow(ilmos(k),jlmos(k),l)
        leaflitrgat(k,l) = leaflitrrow(ilmos(k),jlmos(k),l)
        roottempgat(k,l) = roottemprow(ilmos(k),jlmos(k),l)
        afrleafgat(k,l)  = afrleafrow(ilmos(k),jlmos(k),l)
        afrstemgat(k,l)  = afrstemrow(ilmos(k),jlmos(k),l)
        afrrootgat(k,l)  = afrrootrow(ilmos(k),jlmos(k),l)
        wtstatusgat(k,l) = wtstatusrow(ilmos(k),jlmos(k),l)
        ltstatusgat(k,l) = ltstatusrow(ilmos(k),jlmos(k),l)
        rmlvegaccgat(k,l) = rmlvegaccrow(ilmos(k),jlmos(k),l)
        rmsveggat(k,l)   = rmsvegrow(ilmos(k),jlmos(k),l)
        rmrveggat(k,l)   = rmrvegrow(ilmos(k),jlmos(k),l)
        rgveggat(k,l)    = rgvegrow(ilmos(k),jlmos(k),l)
        gppveggat(k,l)   = gppvegrow(ilmos(k),jlmos(k),l)
        autoresveggat(k,l) = autoresvegrow(ilmos(k),jlmos(k),l)
        vgbiomas_veggat(k,l) = vgbiomas_vegrow(ilmos(k),jlmos(k),l)
        pftexistgat(k,l) = pftexistrow(ilmos(k),jlmos(k),l)

        !         fire emission variables
        emit_co2gat(k,l)  = emit_co2row(ilmos(k),jlmos(k),l)
        emit_cogat(k,l)   = emit_corow(ilmos(k),jlmos(k),l)
        emit_ch4gat(k,l)  = emit_ch4row(ilmos(k),jlmos(k),l)
        emit_nmhcgat(k,l) = emit_nmhcrow(ilmos(k),jlmos(k),l)
        emit_h2gat(k,l)   = emit_h2row(ilmos(k),jlmos(k),l)
        emit_noxgat(k,l)  = emit_noxrow(ilmos(k),jlmos(k),l)
        emit_n2ogat(k,l)  = emit_n2orow(ilmos(k),jlmos(k),l)
        emit_pm25gat(k,l) = emit_pm25row(ilmos(k),jlmos(k),l)
        emit_tpmgat(k,l)  = emit_tpmrow(ilmos(k),jlmos(k),l)
        emit_tcgat(k,l)   = emit_tcrow(ilmos(k),jlmos(k),l)
        emit_ocgat(k,l)   = emit_ocrow(ilmos(k),jlmos(k),l)
        emit_bcgat(k,l)   = emit_bcrow(ilmos(k),jlmos(k),l)
        burnvegfgat(k,l)  = burnvegfrow(ilmos(k),jlmos(k),l)
        !
      end do
    end do ! loop 101
    !
    do l = 1,iccp1 ! loop 102
      do k = 1,nml
        hetroresveggat(k,l) = hetroresvegrow(ilmos(k),jlmos(k),l)
        nepveggat(k,l)   = nepvegrow(ilmos(k),jlmos(k),l)
        nbpveggat(k,l)   = nbpvegrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 102

    ! COMBAK PERLAY
    do l = 1,iccp2 ! loop 103
      do k = 1,nml
        litrmassgat(k,l) = litrmassrow(ilmos(k),jlmos(k),l)
        soilcmasgat(k,l) = soilcmasrow(ilmos(k),jlmos(k),l)
        litresveggat(k,l) = litresvegrow(ilmos(k),jlmos(k),l)
        soilcresveggat(k,l) = soilcresvegrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 103
    !       do l=1,iccp2
    !        do k=1,nml
    !         do m = 1,ignd
    !          litrmassgat(k,l,m)=litrmassrow(ilmos(k),jlmos(k),l,m)
    !          soilcmasgat(k,l,m)=soilcmasrow(ilmos(k),jlmos(k),l,m)
    !          litresveggat(k,l,m)= litresvegrow(ilmos(k),jlmos(k),l,m)
    !          soilcresveggat(k,l,m)= soilcresvegrow(ilmos(k),jlmos(k),l,m)
    !         end do
    ! 103   continue
    ! COMBAK PERLAY
    !
    do l = 1,2 ! loop 106
      do k = 1,nml
        colddaysgat(k,l) = colddaysrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 106
    !
    do l = 1,ican ! loop 201
      do k = 1,nml
        !         ailcgat(k,l)     = ailcrow(ilmos(k),jlmos(k),l)
        !         zolncgat(k,l)    = zolncrow(ilmos(k),jlmos(k),l)
        !         cmasvegcgat(k,l) = cmasvegcrow(ilmos(k),jlmos(k),l)
        !         paicgat(k,l)     = paicrow(ilmos(k),jlmos(k),l)
        !         slaicgat(k,l)    = slaicrow(ilmos(k),jlmos(k),l)
        !
        !         alvsctmgat(k,l)  = alvsctmrow(ilmos(k),jlmos(k),l)
        !         alirctmgat(k,l)  = alirctmrow(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 201
    !
    do l = 1,ignd ! loop 250
      do k = 1,nml
        sandgat(k,l) = sandrow(ilmos(k),jlmos(k),l)
        claygat(k,l) = clayrow(ilmos(k),jlmos(k),l)
        orgmgat(k,l) = orgmrow(ilmos(k),jlmos(k),l)
        tbaraccgat_m(k,l) = tbaraccrow_m(ilmos(k),jlmos(k),l)
        ! thlqaccgat_m(k,l) = thlqaccrow_m(ilmos(k),jlmos(k),l)
        ! thicaccgat_m(k,l) = thicaccrow_m(ilmos(k),jlmos(k),l)
      end do
    end do ! loop 250
    !
    do l = 1,icc ! loop 280
      do m = 1,ignd
        do k = 1,nml
          !         rmatctemgat(k,l,m) = rmatctemrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 280
    !
    do l = 1,ican ! loop 290
      do m = 1,ignd
        do k = 1,nml
          !         rmatcgat(k,l,m) = rmatcrow(ilmos(k),jlmos(k),l,m)
        end do
      end do
    end do ! loop 290
    !
    do k = 1,nml ! loop 300
      ! If peatlands are not present these just gather regardless.
      anmossgat(k) =  anmossrow(ilmos(k),jlmos(k))
      rmlmossgat(k) = rmlmossrow(ilmos(k),jlmos(k))
      gppmossgat(k) = gppmossrow(ilmos(k),jlmos(k))
      armossgat(k) =  armossrow(ilmos(k),jlmos(k))
      nppmossgat(k) = nppmossrow(ilmos(k),jlmos(k))
      peatdepgat(k)    =  peatdeprow(ilmos(k),jlmos(k))
      litrmsmossgat(k) =  litrmsmossrow(ilmos(k),jlmos(k))
      Cmossmasgat(k) = Cmossmasrow(ilmos(k),jlmos(k))
      dmossgat(k) = dmossrow(ilmos(k),jlmos(k))
      ipeatlandgat(k) = ipeatlandrow(ilmos(k),jlmos(k))
      pddgat(k)   = pddrow(ilmos(k),jlmos(k))
    end do ! loop 300

    return

  end subroutine ctemg2
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \namespace ctemgatherscatter
  !! Transfers information between the 'gathered' and 'scattered' form of the CTEM data arrays.
  !!

end module ctemGatherScatter
