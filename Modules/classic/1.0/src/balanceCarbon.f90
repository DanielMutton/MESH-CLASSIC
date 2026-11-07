!> \file
!> Checks carbon pools for conservation
module balanceCarbon

  implicit none

  ! Subroutines
  public :: balcar
  public :: prepBalanceC

contains

  ! ------------------------------------------------------------------
  !> \ingroup balancecarbon_balcar
  !! @{
  !> Check C budget by going through each pool for each vegetation type.
  !! Unless mentioned all pools are in kg C/m2 and all fluxes are in units
  !! of u-mol CO2/m2.sec
  !> @author Vivek Arora
  subroutine balcar (gleafmas, stemmass, rootmass, bleafmas, &
                     litrmass, soilcmas, ntchlveg, ntchsveg, &
                     ntchrveg, tltrleaf, tltrstem, tltrroot, &
                     glcaemls, blcaemls, stcaemls, rtcaemls, &
                     ltrcemls, ltresveg, scresveg, humtrsvg, &
                     pglfmass, pblfmass, pstemass, protmass, &
                     plitmass, psocmass, vgbiomas, repro_cost, &
                     pvgbioms, gavgltms, pgavltms, gavgscms, &
                     pgavscms, galtcels, repro_cost_g, autores, &
                     hetrores, gpp, litres, socres, dstcemls, &
                     litrfall, humiftrs, il1, il2, ilg, &
                     ipeatland, Cmossmas, pCmossmas, &
                     nppmosstep, litrfallmoss, litrmsmoss, &
                     plitrmsmoss, ltrestepmoss, humicmosstep)

    !     -----------------------------------------------------------------
    !
    !     07  Dec 2018  - Pass ilg back in as an argument
    !     V. Arora
    !
    !     22  Nov 2012  - calling this version 1.1 since a fair bit of ctem
    !     V. Arora        subroutines were changed for compatibility with class
    !                     version 3.6 including the capability to run ctem in
    !                     mosaic/tile version along with class.
    !
    !     24  Sep 2012  - add in checks to prevent calculation of non-present
    !     J. Melton       pfts
    !
    !     27  May 2003  - this subroutine checks if the various c fluxes
    !     V. Arora        between the different pools balance properly to
    !                     make sure that conservation of mass is achieved
    !                     with in a specified tolerance.
    !
    use classicParams,      only : tolrance, icc, deltat, &
                                   iccp2, ignd, iccp1
    !
    implicit none
    !
    integer, intent(in) :: ilg !<
    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg
    integer :: i, j, k
    !
    real, intent(in) :: stemmass(ilg,icc)  !< pools (after being updated): stem mass for each of the 9 ctem pfts
    real, intent(in) :: rootmass(ilg,icc)  !< pools (after being updated): root mass for each of the 9 ctem pfts
    real, intent(in) :: gleafmas(ilg,icc)  !< pools (after being updated): green leaf mass for each of the 9 ctem pfts
    real, intent(in) :: bleafmas(ilg,icc)  !< pools (after being updated): brown leaf mass for each of the 9 ctem pfts
    ! real :: litrmass(ilg,iccp2,ignd)      !< pools (after being updated): litter mass over the 9 pfts and the bare fraction of the grid cell ! COMBAK PERLAY
    ! real :: soilcmas(ilg,iccp2,ignd)      !< pools (after being updated): soil carbon mass over the 9 pfts and the bare fraction of the grid cell ! COMBAK PERLAY
    real, intent(in) :: litrmass(ilg,iccp2)!< pools (after being updated): litter mass over the 9 pfts and the bare fraction of the grid cell
    real, intent(in) :: soilcmas(ilg,iccp2)!< pools (after being updated): soil carbon mass over the 9 pfts and the bare fraction of the grid cell
    real, intent(in) :: ntchlveg(ilg,icc)  !< fluxes for each pft: net change in leaf biomass
    real, intent(in) :: ntchsveg(ilg,icc)  !< fluxes for each pft: net change in stem biomass
    real, intent(in) :: ntchrveg(ilg,icc)  !< fluxes for each pft: net change in root biomass the net change is the difference
    !< between allocation and autotrophic respiratory fluxes
    real, intent(in) :: tltrleaf(ilg,icc)  !< fluxes for each pft: total leaf litter falling rate
    real, intent(in) :: tltrstem(ilg,icc)  !< fluxes for each pft: total stem litter falling rate
    real, intent(in) :: tltrroot(ilg,icc)  !< fluxes for each pft: total root litter falling rate
    real, intent(in) :: glcaemls(ilg,icc)  !< fluxes for each pft: carbon emission losses mainly due to fire: green leaf carbon emission losses
    real, intent(in) :: blcaemls(ilg,icc)  !< fluxes for each pft: carbon emission losses mainly due to fire: brown leaf carbon emission losses
    real, intent(in) :: stcaemls(ilg,icc)  !< fluxes for each pft: carbon emission losses mainly due to fire: stem carbon emission losses
    real, intent(in) :: rtcaemls(ilg,icc)  !< fluxes for each pft: carbon emission losses mainly due to fire: root carbon emission losses
    real, intent(in) :: ltrcemls(ilg,icc)  !< fluxes for each pft: carbon emission losses mainly due to fire: litter carbon emission losses
    ! COMBAK PERLAY
    real, intent(in) :: ltresveg(ilg,iccp2)!< fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(in) :: scresveg(ilg,iccp2)!< fluxes for each pft: soil c respiration for each pft + bare fraction
    real, intent(in) :: humtrsvg(ilg,iccp2)!< fluxes for each pft: humification for each pft + bare fraction
    ! real :: ltresveg(ilg,iccp2,ignd)!< fluxes for each pft: litter respiration for each pft + bare fraction
    ! real :: scresveg(ilg,iccp2,ignd)!< fluxes for each pft: soil c respiration for each pft + bare fraction
    ! real :: humtrsvg(ilg,iccp2,ignd)!< fluxes for each pft: humification for each pft + bare fraction
    ! COMBAK PERLAY
    real, intent(in) :: pglfmass(ilg,icc)  !< pools (before being updated): previous green leaf mass
    real, intent(in) :: pblfmass(ilg,icc)  !< pools (before being updated): previous brown leaf mass
    real, intent(in) :: pstemass(ilg,icc)  !< pools (before being updated): previous stem mass
    real, intent(in) :: protmass(ilg,icc)  !< pools (before being updated): previous root mass
    real, intent(in) :: plitmass(ilg,iccp2)!< pools (before being updated): previous litter mass
    real, intent(in) :: psocmass(ilg,iccp2)!< pools (before being updated): previous soil c mass
    real, intent(in) :: vgbiomas(ilg)      !< pools (after being updated): grid averaged pools: vegetation biomass
    real, intent(in) :: pvgbioms(ilg)      !< pools (before being updated): grid average pools: previous vegetation biomass
    real, intent(in) :: gavgltms(ilg)      !< pools (after being updated): grid averaged pools: litter mass
    real, intent(in) :: pgavltms(ilg)      !< pools (before being updated): grid average pools: previous litter mass
    real, intent(in) :: gavgscms(ilg)      !< pools (after being updated): grid averaged pools: soil carbon mass
    real, intent(in) :: pgavscms(ilg)      !< pools (before being updated): grid average pools: previous soil c mass
    real, intent(in) :: autores(ilg)       !< grid averaged flux: autotrophic respiration
    real, intent(in) :: hetrores(ilg)      !< grid averaged flux: heterotrophic respiration
    real, intent(in) :: gpp(ilg)           !< grid averaged flux: gross primary productivity
    real, intent(in) :: litres(ilg)        !< grid averaged flux: litter respiration
    real, intent(in) :: socres(ilg)        !< grid averaged flux: soil carbon respiration
    real, intent(in) :: dstcemls(ilg)      !< grid averaged flux: carbon emission losses due to disturbance, mainly fire
    real, intent(in) :: litrfall(ilg)      !< grid averaged flux: combined (leaves, stem, and root) total litter fall rate
    real, intent(in) :: humiftrs(ilg)      !< grid averaged flux: humification
    real, intent(in) :: repro_cost(ilg,icc)!< pools (after being updated): amount of C transferred to litter due to reproductive tissues
    integer, intent(in) ::  ipeatland (ilg)!< Peatland flag, non-peatlands = 0
    real, intent(in) :: Cmossmas(ilg)      !< moss biomass C (kgC/m2)
    real, intent(in) :: pCmossmas(ilg)     !< moss biomass C at the previous time step (kgC/m2)
    real, intent(in) :: nppmosstep(ilg)    !< moss npp (kgC/m2/timestep)
    real, intent(in) :: litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
    real, intent(in) :: litrmsmoss(ilg)    !< moss litter C (kgC/m2)
    real, intent(in) :: plitrmsmoss(ilg)   !< moss litter C at the previous time step (kgC/m2)
    real, intent(in) :: ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
    real, intent(in) :: humicmosstep(ilg)  !< moss humification (kgC/m2/timestep)
    real, intent(in) :: galtcels(ilg)      !< grid averaged flux: carbon emission losses from litter
    real, intent(in) :: repro_cost_g(ilg)  !< grid averaged flux: amount of C used to generate reproductive tissues
    !
    real ::  soiltempor
    real ::  littempor
    real ::  humtrstemp
    real :: scresveg_temp
    real :: litrestemp
    real :: diff1
    real :: diff2
    !
    !>
    !! To check C budget we go through each pool for each vegetation type.
    !!
    !! Green and brown leaves
    !!
    do j = 1,icc ! loop 100
      do i = il1,il2 ! loop 110
        diff1 = (gleafmas(i,j) + bleafmas(i,j) - pglfmass(i,j) - &
                pblfmass(i,j))
        diff2 = (ntchlveg(i,j) - tltrleaf(i,j) - glcaemls(i,j) - &
                blcaemls(i,j)) * (deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2000)i,j,abs(diff1 - diff2),tolrance
2000      format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
  than our tolerance of ',f12.6,' for leaves')
          call errorHandler('balcar', - 1)
        end if
        !         end if
      end do ! loop 110
    end do ! loop 100
    !
    ! Stem
    !
    do j = 1,icc ! loop 150
      do i = il1,il2 ! loop 160
        diff1 = stemmass(i,j) - pstemass(i,j)
        diff2 = (ntchsveg(i,j) - tltrstem(i,j) - &
                stcaemls(i,j)) * (deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2001)i,j,abs(diff1 - diff2),tolrance
2001      format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
  than our tolerance of ',f12.6,' for stem')
          call errorHandler('balcar', - 2)
        end if
        !         end if
      end do ! loop 160
    end do ! loop 150
    !
    ! Root
    !
    do j = 1,icc ! loop 200
      do i = il1,il2 ! loop 210
        diff1 = rootmass(i,j) - protmass(i,j)
        diff2 = (ntchrveg(i,j) - tltrroot(i,j) - &
                rtcaemls(i,j)) * (deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2002)i,j,abs(diff1 - diff2),tolrance
2002      format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
  than our tolerance of ',f12.6,' for root')
          call errorHandler('balcar', - 3)
        end if
        !         end if
      end do ! loop 210
    end do ! loop 200
    !
    ! Litter over all pfts
    !
    do j = 1,icc ! loop 250
      do i = il1,il2 ! loop 260
        littempor = 0.
        litrestemp = 0.
        humtrstemp = 0.
        ! COMBAK PERLAY
        littempor = littempor + litrmass(i,j)
        litrestemp = litrestemp + ltresveg(i,j)
        humtrstemp = humtrstemp + humtrsvg(i,j)
        ! do k = 1,ignd
        !  littempor = littempor + litrmass(i,j,k)
        !  litrestemp = litrestemp + ltresveg(i,j,k)
        !  humtrstemp = humtrstemp + humtrsvg(i,j,k)
        ! end do
        ! COMBAK PERLAY
        diff1 = littempor - plitmass(i,j)
        diff2 = ( tltrleaf(i,j) + tltrstem(i,j) + tltrroot(i,j) &
                - litrestemp - humtrstemp - ltrcemls(i,j) &
                + repro_cost(i,j)) * (deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2003)i,j,abs(diff1 - diff2),tolrance
2003      format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
   than our tolerance of ',f12.6,' for litter')
          call errorHandler('balcar', - 4)
        end if
      end do ! loop 260
    end do ! loop 250
    !
    ! Litter over the bare fraction
    !
    do i = il1,il2 ! loop 280
      if (ipeatland(i) == 0) then ! Over the non-peatland areas.
        littempor = 0.
        litrestemp = 0.
        humtrstemp = 0.
        ! COMBAK PERLAY
        littempor = littempor + litrmass(i,iccp1)
        litrestemp = litrestemp + ltresveg(i,iccp1)
        humtrstemp = humtrstemp + humtrsvg(i,iccp1)
        ! do k = 1,ignd
        !   littempor = littempor + litrmass(i,iccp1,k)
        !   litrestemp = litrestemp + ltresveg(i,iccp1,k)
        !   humtrstemp = humtrstemp + humtrsvg(i,iccp1,k)
        ! end do
        ! COMBAK PERLAY
        diff1 = littempor - plitmass(i,iccp1)
        diff2 = ( - litrestemp - humtrstemp) * ( deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2003)i,iccp1,abs(diff1 - diff2),tolrance
          call errorHandler('balcar', - 5)
        end if
      end if
    end do ! loop 280

    !
    ! Litter over the LUC pool
    !
    do i = il1,il2 ! loop 290
      littempor = 0.
      litrestemp = 0.
      humtrstemp = 0.
      ! COMBAK PERLAY
      littempor = littempor + litrmass(i,iccp2)
      litrestemp = litrestemp + ltresveg(i,iccp2)
      humtrstemp = humtrstemp + humtrsvg(i,iccp2)
      ! do k = 1,ignd
      !   littempor = littempor + litrmass(i,iccp2,k)
      !   litrestemp = litrestemp + ltresveg(i,iccp2,k)
      !   humtrstemp = humtrstemp + humtrsvg(i,iccp2,k)
      ! end do
      ! COMBAK PERLAY
      diff1 = littempor - plitmass(i,iccp2)
      diff2 = ( - litrestemp - humtrstemp) * ( deltat/963.62)

      if ((abs(diff1 - diff2)) > tolrance) then
        write(6,2003)i,iccp2,abs(diff1 - diff2),tolrance
        call errorHandler('balcar', - 6)
      end if
    end do ! loop 290

    !
    ! Soil carbon for the vegetated areas
    !
    do j = 1,icc ! loop 300
      do i = il1,il2 ! loop 310
        if (ipeatland(i) == 0) then ! Over the non-peatland regions
          soiltempor = 0.
          scresveg_temp = 0.
          humtrstemp = 0.
          ! COMBAK PERLAY
          soiltempor = soilcmas(i,j)
          scresveg_temp = scresveg(i,j)
          humtrstemp = humtrsvg(i,j)
          ! do k = 1,ignd
          !   soiltempor = soiltempor + soilcmas(i,j,k)
          !   scresveg_temp = scresveg_temp + scresveg(i,j,k)
          !   humtrstemp = humtrstemp + humtrsvg(i,j,k)
          ! end do
          ! COMBAK PERLAY
          diff1 = soiltempor - psocmass(i,j)
          diff2 = ( humtrstemp - scresveg_temp) * (deltat/963.62)
          if ((abs(diff1 - diff2)) > tolrance) then
            ! write(6,3001)'soilCmas(',i,')=',soilcmas(i,j)
            ! write(6,3001)'psocmass(',i,')=',psocmass(i,j)
            write(6,3001)'humtr(',i,') = ',humtrstemp * (deltat/963.62)
            ! write(6,3001)'scres(',i,')=',scresveg(i,j)*(deltat/963.62)
            write(6,2004)i,j,abs(diff1 - diff2),tolrance
2004        format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
  than our tolerance of ',f12.6,' for soil c')
            call errorHandler('balcar', - 7)
          end if
        end if
      end do ! loop 310
    end do ! loop 300

    !
    ! Soil carbon over the bare area and LUC product pool
    !
    do j = iccp1,iccp2 ! loop 320
      do i = il1,il2 ! loop 330
        soiltempor = 0.
        scresveg_temp = 0.
        humtrstemp = 0.
        ! COMBAK PERLAY
        soiltempor = soilcmas(i,j)
        scresveg_temp = scresveg(i,j)
        humtrstemp = humtrsvg(i,j)
        ! do k = 1,ignd
        !   soiltempor = soiltempor + soilcmas(i,j,k)
        !   scresveg_temp = scresveg_temp + scresveg(i,j,k)
        !   humtrstemp = humtrstemp + humtrsvg(i,j,k)
        ! end do
        ! COMBAK PERLAY
        diff1 = soiltempor - psocmass(i,j)
        diff2 = ( humtrstemp - scresveg_temp) * (deltat/963.62)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,2004)i,j,abs(diff1 - diff2),tolrance
2014      format('at (i) = (',i3,'),pft = ',i2,', ',f12.6,' is greater &
  than our tolerance of ',f12.6,' for soil c')
          call errorHandler('balcar', - 8)
        end if
      end do ! loop 330
    end do ! loop 320

    !
    ! Moss C balance
    !
    do i = il1,il2 ! loop 400
      if (ipeatland(i) > 0) then ! Peatlands only
        diff1 = Cmossmas(i) - pCmossmas(i)
        diff2 = nppmosstep(i) - litrfallmoss(i)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,3001)'Cmossmas(',i,') = ',Cmossmas(i)
          write(6,3001)'pCmossmas(',i,') = ',pCmossmas(i)
          write(6,3001)'nppmosstep(',i,') = ',nppmosstep(i)
          write(6,3001)' litrfallmoss(',i,') = ',litrfallmoss(i)
          write(6,2008)i,abs(diff1 - diff2),tolrance
2008      format('at (i) = (',i3,'),',f12.6,' is greater' &
          'than our tolerance of ',f12.6,' for moss carbon')
          call errorHandler('balcar', - 9)
        end if
        !
        ! Moss litter pool C balance
        !
        diff1 = litrmsmoss(i) - plitrmsmoss(i)
        diff2 = litrfallmoss(i) - ltrestepmoss(i) - humicmosstep(i)
        if ((abs(diff1 - diff2)) > tolrance) then
          write(6,3001)'litrmsmoss(',i,') = ',litrmsmoss(i)
          write(6,3001)'plitrmsmoss(',i,') = ',plitrmsmoss(i)
          write(6,3001)'litrfallmoss(',i,') = ',litrfallmoss(i)
          write(6,3001)' ltrestepmoss(',i,') = ',ltrestepmoss(i)
          write(6,3001)' humicmosstep(',i,') = ',humicmosstep(i)
          write(6,2009)i,abs(diff1 - diff2),tolrance
2009      format('at (i) = (',i3,'),',f12.6,' is greater &
 than our tolerance of ',f12.6,' for moss litter')
          call errorHandler('balcar', - 10)
        end if
      end if
    end do ! loop 400
    !
    ! Grid averaged fluxes must also balance
    !
    ! Vegetation biomass
    !
    do i = il1,il2 ! loop 350
      diff1 = vgbiomas(i) - pvgbioms(i)
      diff2 = (gpp(i) - autores(i) - litrfall(i) - &
              dstcemls(i) - repro_cost_g(i)) * (deltat/963.62)
      if ((abs(diff1 - diff2)) > (tolrance + 0.00003)) then !   then ! YW the difference for moss litter
        ! and biomass can go a bit over the tolerance
        write(6,3001)'vgbiomas(',i,') = ',vgbiomas(i)
        write(6,3001)'pvgbioms(',i,') = ',pvgbioms(i)
        write(6,3001)'     gpp(',i,') = ',gpp(i)
        write(6,3001)' autores(',i,') = ',autores(i)
        write(6,3001)'litrfall(',i,') = ',litrfall(i)
        write(6,3001)'dstcemls(',i,') = ',dstcemls(i)
        write(6,3001)'repro_cost_g(',i,') = ',repro_cost_g(i)
3001    format(a9,i2,a2,f14.9)
        write(6,2005)i,abs(diff1 - diff2),tolrance
2005    format('at (i) = (',i3,'),',f12.6,' is greater &
  than our tolerance of ',f12.6,' for vegetation biomass')
        write(90, * )    abs(diff1 - diff2),tolrance
        call errorHandler('balcar', - 11)
      end if
    end do ! loop 350
    !
    ! Litter
    !
    do i = il1,il2 ! loop 380
      diff1 = gavgltms(i) - pgavltms(i)
      diff2 = (litrfall(i) - litres(i) - humiftrs(i) - galtcels(i) &
              + repro_cost_g(i)) * &
              (deltat/963.62)
      if ((abs(diff1 - diff2)) > (tolrance + 0.00003)) then    ! YW the difference for moss litter
        ! and biomass can go a bit over the tolerance
        write(6,3001)'pgavltms(',i,') = ',pgavltms(i)
        write(6,3001)'gavgltms(',i,') = ',gavgltms(i)
        write(6,3001)'litrfall(',i,') = ',litrfall(i)
        write(6,3001)'  litres(',i,') = ',litres(i)
        write(6,3001)'humiftrs(',i,') = ',humiftrs(i)
        write(6,3001)'galtcels(',i,') = ',galtcels(i)
        write(6,2006)i,abs(diff1 - diff2),tolrance
        write( * , * )i,abs(diff1 - diff2),tolrance
2006    format('at (i) = (',i3,'),',f12.6,' is greater &
  than our tolerance of ',f12.6,' for litter mass')
        call errorHandler('balcar', - 12)
      end if
    end do ! loop 380
    !
    ! Soil carbon
    !
    do i = il1,il2 ! loop 390
      diff1 = gavgscms(i) - pgavscms(i)
      diff2 = (humiftrs(i) - socres(i)) * (deltat/963.62)
      if ((abs(diff1 - diff2)) > tolrance) then
        write(6,3001)'pgavscms(',i,') = ',pgavscms(i)
        write(6,3001)'gavgscms(',i,') = ',gavgscms(i)
        write(6,3001)'humiftrs(',i,') = ',humiftrs(i)
        write(6,3001)'  socres(',i,') = ',socres(i)
        write(6,2007)i,abs(diff1 - diff2),tolrance
2007    format('at (i) = (',i3,'),',f12.6,' is greater &
  than our tolerance of ',f12.6,' for soil c mass')
        call errorHandler('balcar', - 13)
      end if
    end do ! loop 390

    return
  end subroutine balcar
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup balancecarbon_prepBalanceC
  !! @{
  !> Prepare for the carbon balance check. Calculate total litter fall from each
  !! component (leaves, stem, and root) from all causes (normal turnover, drought
  !! and cold stress for leaves, mortality, and disturbance), calculate grid-average
  !! vegetation biomass, litter mass, and soil carbon mass, and litter fall rate.
  !! Also add the bare ground values to the grid-average. If a peatland, we assume no bareground and
  !! add the moss values instead. Note: peatland soil C is not aggregated from plants but updated
  !! by humification and respiration from the previous stored value
  !> @author J. Melton

  subroutine prepBalanceC (il1, il2, ilg, fcancmx, glealtrm, glfltrdt, & ! In
                           blfltrdt, stemltrm, stemltdt, rootltrm, rootltdt, & ! In
                           ipeatland, nppmosstep, pgavscms, humstepmoss, & ! In
                           ltrestepmoss, stemlitr, rootlitr, rootmass, & ! In
                           litrmass, soilCmas, hutrstep_g, stemmass, bleafmas, & ! In
                           gleafmas, socrestep, fg, litrfallmoss, & ! In
                           leaflitr, Cmossmas, litrmsmoss, & ! In/Out
                           tltrleaf, tltrstem, tltrroot, vgbiomas, litrfall, & ! Out
                           gavgltms, litrfallveg, gavgscms, vgbiomas_veg)        ! Out

    use classicParams,     only : icc, deltat, iccp1, ignd

    implicit none

    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg             !< il2=ilg (no. of grid cells in latitude circle)
    real, intent(in)    :: fcancmx(:,:)    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !! modified by land-use change,and competition between pfts
    real, intent(in) :: glealtrm(:,:) !< green leaf litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in) :: glfltrdt(:,:) !< green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(in) :: blfltrdt(:,:) !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(in) :: stemlitr(:,:) !< stem litter \f$(kg C/m^2)\f$
    real, intent(in) :: stemltrm(:,:) !< stem litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in) :: stemltdt(:,:) !< stem litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(in) :: rootlitr(:,:) !< root litter \f$(kg C/m^2)\f$
    real, intent(in) :: rootltrm(:,:) !< root litter generated due to mortality \f$(kg C/m^2)\f$
    real, intent(in) :: rootltdt(:,:) !< root litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(in) :: bleafmas(:,:) !< brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: gleafmas(:,:) !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: rootmass(:,:) !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: stemmass(:,:) !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! real, intent(in) :: litrmass(:,:,:) ! COMBAK PERLAY
    ! real, intent(in) :: soilCmas(:,:,:) ! COMBAK PERLAY
    real, intent(in) :: litrmass(:,:) !< litter mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: soilCmas(:,:) !<
    real, intent(in) :: hutrstep_g(:) !< grid sum of humification from vascualr litter \f$(kg C/m^2/timestep)\f$
    real, intent(in) :: socrestep(:) !< heterotrophic respiration from soil \f$(kg C/m^2/timestep)\f$
    integer, intent(in) :: ipeatland(:) !< Peatland switch: 0 = not a peatland, 1 = bog, 2 = fen
    real, intent(in) :: nppmosstep(:) !< moss npp \f$(kg C/m^2/timestep)\f$
    real, intent(in) :: humstepmoss(:) !< moss humification \f$(kg C/m^2/timestep)\f$
    real, intent(in) :: ltrestepmoss(:) !< litter respiration from moss \f$(kg C/m^2/timestep)\f$
    real, intent(in) :: pgavscms(:) !<
    real, intent(in) :: fg(:) !< Fraction of grid cell that is bare ground.
    real, intent(in) :: litrfallmoss(:) !< moss litter fall \f$(kg C/m^2/timestep)\f$

    real, intent(inout) :: leaflitr(:,:) !<
    real, intent(inout) :: Cmossmas(:) !< C in moss biomass, \f$kg C/m^2\f$
    real, intent(inout) :: litrmsmoss(:) !< moss litter mass, \f$kg C/m^2\f$

    real, intent(out) :: tltrleaf(ilg,icc) !< total leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tltrstem(ilg,icc) !< total stem litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tltrroot(ilg,icc) !< total root litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: litrfallveg(ilg,icc) !< litter fall for each pft (\f$kg c/m^2\f$)
    real, intent(out) :: vgbiomas(ilg) !< grid averaged vegetation biomass (\f$kg C/m^2\f$)
    real, intent(out) :: vgbiomas_veg(ilg,icc) !< vegetation biomass for each pft
    real, intent(out) :: litrfall(ilg) !< total litter fall (from leaves, stem, and root) due to
    !! all causes (mortality,turnover,and disturbance)
    real, intent(out) :: gavgltms(ilg) !< grid averaged litter mass including the LUC product pool, (\f$kg C/m^2\f$)
    real, intent(out) :: gavgscms(ilg) !< grid averaged soil c mass, (\f$kg C/m^2\f$)

    ! Local
    integer :: i, j, k

    !>
    !! Calculate total litter fall from each component (leaves, stem, and root)
    !!  from all causes (normal turnover, drought and cold stress for leaves, mortality,
    !! and disturbance) for use in balcar subroutine
    !!
    do j = 1,icc ! loop 1050
      do i = il1,il2 ! loop 1060

        ! units here are \f$kg c/m^2 .day\f$
        tltrleaf(i,j) = leaflitr(i,j) + glealtrm(i,j) + glfltrdt(i,j) &
                        + blfltrdt(i,j)
        tltrstem(i,j) = stemlitr(i,j) + stemltrm(i,j) + stemltdt(i,j)
        tltrroot(i,j) = rootlitr(i,j) + rootltrm(i,j) + rootltdt(i,j)
        !>
        !> convert units to u-mol co2/m2.sec
        leaflitr(i,j) = leaflitr(i,j) * (963.62 / deltat)
        tltrleaf(i,j) = tltrleaf(i,j) * (963.62 / deltat)
        tltrstem(i,j) = tltrstem(i,j) * (963.62 / deltat)
        tltrroot(i,j) = tltrroot(i,j) * (963.62 / deltat)
      end do ! loop 1060
    end do ! loop 1050
    !>
    !> calculate grid-average vegetation biomass, litter mass, and soil carbon mass, and litter fall rate
    !>
    litrfall(:) = 0.0
    vgbiomas(:) = 0.0
    gavgltms(:) = 0.0
    gavgscms(:) = 0.0
    do j = 1,icc ! loop 1100
      do i = il1,il2 ! loop 1110

        vgbiomas(i) = vgbiomas(i) + fcancmx(i,j) * (gleafmas(i,j) &
                      + bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))

        litrfall(i) = litrfall(i) + fcancmx(i,j) * (tltrleaf(i,j) &
                      + tltrstem(i,j) + tltrroot(i,j))

        ! Store the per PFT litterfall for outputting.
        litrfallveg(i,j) = (tltrleaf(i,j) + tltrstem(i,j) + tltrroot(i,j))

        ! COMBAK PERLAY
        gavgltms(i) = gavgltms(i) + fcancmx(i,j) * litrmass(i,j)
        if (ipeatland(i) == 0) then ! Non-peatlands
          gavgscms(i) = gavgscms(i) + fcancmx(i,j) * soilcmas(i,j)
          ! else
          ! Peatland soil C is calculated from peat depth (peatdep) in the peatland
        end if
        ! do k = 1,ignd
        !   gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j,k)
        !   if (ipeatland(i)==0) then ! Non-peatlands
        !     gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j,k)
        !   ! else
        !     ! Peatland soil C is calculated from peat depth (peatdep) in the peatland
        !   end if
        ! end do
        ! COMBAK PERLAY

        vgbiomas_veg(i,j) = gleafmas(i,j) + bleafmas(i,j) + stemmass(i,j) &
                            + rootmass(i,j)
      end do ! loop 1110
    end do ! loop 1100
    !
    !> Add the bare ground values to the grid-average. If a peatland, we assume no bareground and
    !! add the moss values instead.
    !! Note: peatland soil C is not aggregated from plants but updated
    !! by humification and respiration from the previous stored value

    do i = il1,il2 ! loop 1020
      if (ipeatland(i) == 0) then
        ! COMBAK PERLAY
        ! Add the bare fraction dead C
        gavgltms(i) = gavgltms(i) + fg(i) * litrmass(i,iccp1)
        gavgscms(i) = gavgscms(i) + fg(i) * soilcmas(i,iccp1)
        ! do k  = 1,ignd
        !   ! Add the bare fraction dead C
        !   gavgltms(i)=gavgltms(i)+ fg(i) *litrmass(i,iccp1,k)
        !   gavgscms(i)=gavgscms(i)+ fg(i) * soilcmas(i,iccp1,k)
        ! 1025    continue
        ! COMBAK PERLAY
      else
        litrmsmoss(i) = litrmsmoss(i) + litrfallmoss(i) &  ! FLAG, I am thinking this needs to be reset to zero somewhere in ctemDriver, but I don't see it. JM.
                        - ltrestepmoss(i) - humstepmoss(i)     ! kg/m2
        Cmossmas(i) = Cmossmas(i) + nppmosstep(i) - litrfallmoss(i)
        vgbiomas(i) = vgbiomas(i) + Cmossmas(i)
        litrfall(i) = litrfall(i) + litrfallmoss(i) * (963.62 / deltat)! umolCO2/m2/s
        gavgltms(i) = gavgltms(i) + litrmsmoss(i)
        gavgscms(i) = pgavscms(i) + hutrstep_g(i) - socrestep(i)
      end if
    end do ! loop 1020

  end subroutine prepBalanceC
  !! @}
  !> \namespace balancecarbon
  !! Checks carbon pools for conservation

end module balanceCarbon
