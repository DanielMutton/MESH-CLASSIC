!> \file
!> Calculates the litter generated from stem and root turnover
!! Also updates the vegetation pools for changes due to the turnover calculations.
!! @author Vivek Arora, Joe Melton
!!
module turnover

  implicit none

  public :: turnoverStemRoot
  public :: updatePoolsTurnover

contains

  !> \ingroup turnover
  !! @{
  !> Calculates the litter generated from stem and root turnover
  !> @author Vivek Arora and Joe Melton
  subroutine turnoverStemRoot (stemmass, rootmass, lfstatus, ailcg, & ! in
                               il1, il2, ilg, leapnow, useTracer, & ! In
                               sort, fcancmx, tracerStemMass, tracerRootMass, & ! In
                               stmhrlos, rothrlos, & ! In / Out
                               stemlitr, rootlitr, tracerStemLitr, tracerRootLitr) ! Out
    !
    !     06  Dec 2018  - Pass ilg back in as an argument
    !     V. Arora
    !
    !     17  Jan 2014  - Moved parameters to global file (ctem_params.f90)
    !     J. Melton
    !
    !     22  Jul 2013  - Add in module for parameters
    !     J. Melton
    !
    !     24  Sep 2012  - add in checks to prevent calculation of non-present
    !     J. Melton       pfts
    !
    !     07  May 2003  - this subroutine calculates the litter generated
    !     V. Arora        from stem and root turnover

    use classicParams, only : icc, ican, kk, zero, stemlife, &
                                rootlife, stmhrspn, classpfts, &
                                nol2pfts, reindexPFTs

    implicit none

    integer, intent(in) :: ilg !<
    integer, intent(in) :: il1 !< il1=1
    integer, intent(in) :: il2 !< il2=ilg
    integer, intent(in) :: lfstatus(ilg,icc) !< leaf status. an integer :: indicating if leaves are in "max.
    !< growth", "normal growth", "fall/harvest", or "no leaves" mode. see phenolgy subroutine for more details.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    logical, intent(in) :: leapnow        !< true if this year is a leap year. Only used if the switch 'leap' is true.
    integer, intent(in) :: sort(icc)      !< index for correspondence between ctem pfts and size of parameter vectors
    real, intent(in) :: stemmass(ilg,icc) !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: rootmass(ilg,icc) !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: tracerStemMass(:,:)  !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerRootMass(:,:)  !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: ailcg(ilg,icc)    !< green or live lai
    real, intent(in) :: fcancmx(ilg,icc)  !< max. fractional coverage of CTEM pfts, but this can be modified by
    !! land-use change,and competition between pfts
    real, intent(inout) :: rothrlos(ilg,icc) !< root death for crops. when in "harvest" mode for crops, root is assumed
    !! to die in a similar way as stem is harvested.
    real, intent(inout) :: stmhrlos(ilg,icc) !< stem harvest loss for crops. when in "harvest" mode for crops, stem is
    !! also assumed to be harvested and this generates litter.
    real, intent(out) :: stemlitr(ilg,icc) !< stem litter \f$(kg C/m^2)\f$
    real, intent(out) :: rootlitr(ilg,icc) !< root litter \f$(kg C/m^2)\f$
    real, intent(out) :: tracerStemLitr(ilg,icc) !< Tracer stem litter \f$tracer C units/m^2\f$
    real, intent(out) :: tracerRootLitr(ilg,icc) !< Tracer root litter \f$tracer C units/m^2\f$

    integer :: i, j, k, n, m
    real :: nrmlsmlr(ilg, icc) ! stem litter from normal turnover
    real :: nrmlrtlr(ilg, icc) ! root litter from normal turnover
    real :: frac !< temp. var.

    ! Initialize required arrays to zero
    stemlitr = 0.0
    rootlitr = 0.0
    nrmlsmlr = 0.0
    nrmlrtlr = 0.0

    !------------------------------------------------------------------

    !> Calculate normal stem and root litter using the amount of stem and
    !! root biomass and their turnover time scales.
    do j = 1,icc
      n = sort(j)
      do i = il1,il2
        if (fcancmx(i,j) > 0.0) then

          if (stemlife(n) > zero) then
            if (leapnow) then
              nrmlsmlr(i,j) = stemmass(i,j) * (1.0 - exp( - 1.0 / (366.0 * stemlife(n))))
            else
              nrmlsmlr(i,j) = stemmass(i,j) * (1.0 - exp( - 1.0 / (365.0 * stemlife(n))))
            end if
          end if

          if (rootlife(n) > zero) then
            if (leapnow) then
              nrmlrtlr(i,j) = rootmass(i,j) * (1.0 - exp( - 1.0 / (366.0 * rootlife(n))))
            else
              nrmlrtlr(i,j) = rootmass(i,j) * (1.0 - exp( - 1.0 / (365.0 * rootlife(n))))
            end if
          end if
        end if
      end do ! loop 210
    end do ! loop 200

    !> If crops are in harvest mode then we start harvesting stem as well.
    !! If stem has already been harvested then we set the stem harvest
    !! loss equal to zero. the roots of the crop die in a similar way.
    do j = 1,ican
      do m = reindexPFTs(j,1),reindexPFTs(j,2)
        do i = il1,il2
          if (fcancmx(i,m) > 0.0) then
            select case (classpfts(j))

            case ('NdlTr' , 'BdlTr', 'Grass', 'BdlSh')

              stmhrlos(i,m) = 0.0
              rothrlos(i,m) = 0.0

            case ('Crops')

              if (lfstatus(i,m) == 3) then
                if (stmhrlos(i,m) <= zero .and. stemmass(i,m) > zero) stmhrlos(i,m) = stemmass(i,m) * (1.0 / stmhrspn)
                if (rothrlos(i,m) <= zero .and. rootmass(i,m) > zero) rothrlos(i,m) = rootmass(i,m) * (1.0 / stmhrspn)
              end if

              if (stemmass(i,m) <= zero .or. lfstatus(i,m) == 1 .or. lfstatus(i,m) == 2) stmhrlos(i,m) = 0.0
              if (rootmass(i,m) <= zero .or. lfstatus(i,m) == 1 .or. lfstatus(i,m) == 2) rothrlos(i,m) = 0.0

            case default

              print * ,'Unknown CLASS PFT in turnover ',classpfts(j)
              call errorHandler('turnover', - 1)

            end select
          end if
        end do ! loop 260
      end do ! loop 255
    end do ! loop 250

    !> add stem and root litter from all sources
    !>
    do j = 1,icc ! loop 350
      do i = il1,il2 ! loop 360
        if (fcancmx(i,j) > 0.0) then
          stemlitr(i,j) = nrmlsmlr(i,j) + stmhrlos(i,j)
          rootlitr(i,j) = nrmlrtlr(i,j) + rothrlos(i,j)
          if (useTracer > 0) then
            if (stemmass(i,j) > 0.) then
              frac = tracerStemMass(i,j) / stemmass(i,j)
              tracerStemLitr(i,j) = stemlitr(i,j) * frac
            end if

            if (rootmass(i,j) > 0.) then
              frac = tracerRootMass(i,j) / rootmass(i,j)
              tracerRootLitr(i,j) = rootlitr(i,j) * frac
            end if

          end if
        end if
      end do ! loop 360
    end do ! loop 350

    return

  end subroutine turnoverStemRoot
  !! @}

  !> \ingroup turnover
  !! @{
  !> Update green leaf biomass for trees and crops, brown leaf biomass for grasses,
  !! stem and root biomass for litter deductions, and update litter pool with leaf
  !! litter calculated in the phenology subroutine and stem and root litter
  !! calculated in the turnover subroutine. Also add the reproduction
  !!  carbon directly to the litter pool.
  !> @author Vivek Arora and Joe Melton
  subroutine updatePoolsTurnover (il1, il2, ilg, reprocost, rmatctem, useTracer, tracerReproCost, & ! In
                                  stemmass, rootmass, litrmass, rootlitr, & ! In/Out
                                  gleafmas, bleafmas, leaflitr, stemlitr, & ! In/Out
                                  tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Out
                                  tracerRootMass, tracerStemMass, & ! In/Out
                                  tracerLeafLitr, tracerStemLitr, tracerRootLitr) ! In/Out

    use classicParams, only : ican, nol2pfts, classpfts, deltat, icc, ignd, reindexPFTs

    implicit none

    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg             !< no. of grid cells/tiles in latitude circle
    real, intent(in) :: reprocost(:,:)     !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: tracerReproCost(:,:)  !< Tracer cost of making reproductive tissues, only non-zero when
    !! NPP is positive (\f$tracer C units m^{-2} s^{-1}\f$)
    real, intent(in)    :: rmatctem(:,:,:) !< fraction of roots for each of ctem's 9 pfts in each soil layer
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    real, intent(inout) :: rootmass(:,:)   !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: gleafmas(:,:)   !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: bleafmas(:,:)   !< brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(inout) :: stemmass(:,:)   !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    real, intent(inout) :: litrmass(:,:) !< litter mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! real, intent(inout) :: litrmass(:,:,:) !< litter mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    real, intent(inout) :: leaflitr(:,:)   !< leaf litter \f$(kg C/m^2)\f$
    real, intent(inout) :: rootlitr(:,:)   !< root litter \f$(kg C/m^2)\f$
    real, intent(inout) :: stemlitr(:,:)   !< stem litter \f$(kg C/m^2)\f$
    real, intent(inout) :: tracerStemLitr(:,:) !< Tracer stem litter \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootLitr(:,:) !< Tracer root litter \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLeafLitr(ilg,icc)  !< Tracer leaf litter generated by normal turnover, cold
    !< and drought stress, and leaf fall/harvest, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$tracer C units/m^2\f$

    ! Local.
    integer :: j, m, i, k

    ! ------------------------------

    do j = 1,ican ! loop 700
      do m = reindexPFTs(j,1),reindexPFTs(j,2)
        do i = il1,il2
          select case (classpfts(j))

          case ('NdlTr' , 'BdlTr', 'Crops', 'BdlSh') ! tree,crops and shrub

            gleafmas(i,m) = gleafmas(i,m) - leaflitr(i,m)
            if (useTracer > 0) tracerGLeafMass(i,m) = tracerGLeafMass(i,m) - tracerLeafLitr(i,m)

            if (gleafmas(i,m) < 0.0) then
              leaflitr(i,m) = leaflitr(i,m) + gleafmas(i,m)
              gleafmas(i,m) = 0.0
              if (useTracer > 0) then
                tracerLeafLitr(i,m) = tracerLeafLitr(i,m) + tracerGLeafMass(i,m)
                tracerGLeafMass(i,m) = 0.0
              end if
            end if

          case ('Grass')

            bleafmas(i,m) = bleafmas(i,m) - leaflitr(i,m)
            if (useTracer > 0) tracerBLeafMass(i,m) = tracerBLeafMass(i,m) - tracerLeafLitr(i,m)

            if (bleafmas(i,m) < 0.0) then
              leaflitr(i,m) = leaflitr(i,m) + bleafmas(i,m)
              bleafmas(i,m) = 0.0
              if (useTracer > 0) then
                tracerLeafLitr(i,m) = tracerLeafLitr(i,m) + tracerBLeafMass(i,m)
                tracerBLeafMass(i,m) = 0.0
              end if
            end if

          case default

            print * ,'Unknown CLASS PFT in ctem ',classpfts(j)
            call errorHandler('ctem', - 5)

          end select
        end do ! loop 710
      end do ! loop 705
    end do ! loop 700

    !> Update stem and root biomass for litter deductions
    do j = 1,icc ! loop 780
      do i = il1,il2 ! loop 790
        stemmass(i,j) = stemmass(i,j) - stemlitr(i,j)
        if (useTracer > 0) tracerStemMass(i,j) = tracerStemMass(i,j) - tracerStemLitr(i,j)

        if (stemmass(i,j) < 0.0) then
          stemlitr(i,j) = stemlitr(i,j) + stemmass(i,j)
          stemmass(i,j) = 0.0
          if (useTracer > 0) then
            tracerStemLitr(i,j) = tracerStemLitr(i,j) + tracerStemMass(i,j)
            tracerStemMass(i,j) = 0.0
          end if
        end if

        rootmass(i,j) = rootmass(i,j) - rootlitr(i,j)
        if (useTracer > 0) tracerRootMass(i,j) = tracerRootMass(i,j) - tracerRootLitr(i,j)

        if (rootmass(i,j) < 0.0) then
          rootlitr(i,j) = rootlitr(i,j) + rootmass(i,j)
          rootmass(i,j) = 0.0
          if (useTracer > 0) then
            tracerRootLitr(i,j) = tracerRootLitr(i,j) + tracerRootMass(i,j)
            tracerRootMass(i,j) = 0.0
          end if
        end if

      end do ! loop 790
    end do ! loop 780

    !> Update litter pool with leaf litter calculated in the phenology
    !! subroutine and stem and root litter calculated in the turnover
    !! subroutine. Also add the reproduction carbon directly to the litter pool
    do i = il1,il2 ! loop 800
      do j = 1,icc ! loop 805
        ! COMBAK PERLAY
        litrmass(i,j) = litrmass(i,j) + leaflitr(i,j) + stemlitr(i,j) &
                        + rootlitr(i,j) &
                        + reprocost(i,j) * deltat / 963.62

        !       do k = 1,ignd
        !         if (k == 1) then
        !           ! The first layer gets the leaf and stem litter as well as the reprocost,
        !           ! which is assumed to be cones/seeds. The root litter is given in proportion
        !           ! to the root distribution
        !           litrmass(i,j,k) = litrmass(i,j,k) + leaflitr(i,j) + stemlitr(i,j) &
        !                             + rootlitr(i,j) * rmatctem(i,j,k) &
        !                             + reprocost(i,j) * deltat / 963.62
        !
        !           if (useTracer > 0) tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) + tracerLeafLitr(i,j) + tracerStemLitr(i,j) &
        !                               + tracerRootLitr(i,j) * rmatctem(i,j,k) &
        !                               + tracerReproCost(i,j) * deltat / 963.62
        !         else ! the lower soil layers get the roots,in the proportion that they
        !              ! are in the unfrozen soil column.
        !           litrmass(i,j,k)=litrmass(i,j,k) + rootlitr(i,j) * rmatctem(i,j,k)
        !
        !           if (useTracer > 0) tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) + tracerRootLitr(i,j) * rmatctem(i,j,k)
        !         end if
        !
        ! 810     continue
        ! COMBAK PERLAY
      end do ! loop 805
    end do ! loop 800


  end subroutine updatePoolsTurnover
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \namespace turnover
  !!
  !! The turnover of stem and root components is modelled via their PFT-dependent specified lifetimes.
  !! The litter generation (\f$kg\, C\, m^{-2}\f$ \f$day^{-1}\f$) associated with turnover of stem
  !! (\f$D_\mathrm{S}\f$) and root (\f$D_\mathrm{R}\f$) components is calculated based on the
  !! amount of biomass in the respective components (\f$C_\mathrm{S},
  !! C_\mathrm{R}\f$; \f$kg\, C\, m^{-2}\f$) and their respective turnover timescales
  !! (\f$\tau_\mathrm{S}\f$ and \f$\tau_\mathrm{R}\f$; \f$yr\f$; see also classicParams.f90) as
  !!
  !! \f[ D_{i} = C_{i}\left[1 - \exp\left(-\frac{1}{365\, \tau_{i}}\right)\right], \quad
  !! i = S, R.\f]
  !!
  !! Litter contributions are either put in the first soil layer (leaf and stem litter) whereas
  !! root litter is added in proportion to the root distribution to non-perennially frozen soil layers.
  !! For defining which layers are frozen, we use the active layer depth.
  !!

end module turnover
