!> \file
!> Autotrophic respiration. Leaf respiration is calculated in phtsyn subroutine, while stem
!! and root maintenance respiration are calculated here.
!!
module autotrophicRespiration

  implicit none

  public :: mainres
  public :: growthRespiration

contains

  !> \ingroup autotrophic_res_mainres
  !! @{
  !> Calculates maintenance respiration for roots and stems
  !> @author Vivek Arora and Joe Melton
  subroutine mainres (fcan, fct, stemmass, rootmass, il1, il2, ilg, & ! In
                      leapnow, ta, tbar, rmatctem, sort, isand, & ! In
                      useTracer, tracerStemMass, tracerRootMass, & ! In
                      rmsveg, rmrveg, roottemp, rmsTracer, rmrTracer) ! Out
    !
    !     20  sep. 2001 - this subroutine calculates maintenance respiration,
    !     V. Arora        over a given sub-area, for stem and root components.
    !                     leaf respiration is estimated within the phtsyn
    !                     subroutine.

    !     change history:

    !     J. Melton 22  Apr 2019 - Moved to own module. Added tracer.
    !     J. Melton 22  Jul 2015 - Loops 180 and 190 were not set up for > 3 soil layers. Fixed.
    !
    !     J. Melton 17  Jan 2014 - Moved parameters to global file (classicParams.f90)
    !
    !     J. Melton 22  Jul 2013 - Add in module for parameters
    !
    !     J. Melton 20 sep 2012 - made it so does not do calcs for pfts with
    !                             fcan = 0.
    !     J. Melton 23 aug 2012 - change sand to isand, converting sand to
    !                             int was missing some gridcells assigned
    !                             to bedrock in classb. isand is now passed
    !                             in.
    !
    use classicParams,     only : icc, ignd, ican, kk, zero, &
                                  bsrtstem, bsrtroot, minlvfr, &
                                  classpfts, nol2pfts, reindexPFTs

    implicit none

    integer, intent(in) :: ilg          !< Number of grid cells in latitude circle
    integer, intent(in) :: il1          !< il1=1
    integer, intent(in) :: il2          !< il2=ilg
    integer, intent(in) :: sort(icc)    !< index for correspondence between 9 pfts and 12 values in the parameter vectors
    integer, intent(in) :: isand(ilg,ignd) !< flag for bedrock or ice in a soil layer
    logical, intent(in) :: leapnow        !< true if this year is a leap year. Only used if the switch 'leap' is true.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: fcan(ilg,icc)     !< fractional coverage of ctem's 9 pfts over the given sub-area
    real, intent(in) :: fct(ilg)          !< sum of all fcan fcan & fct are not used at this time but could be used at some later stage
    real, intent(in) :: stemmass(ilg,icc) !< stem biomass for the 9 pfts in \f$kg c/m^2\f$
    real, intent(in) :: ta(ilg)           !< Air temperature, K
    real, intent(in) :: tbar(ilg,ignd)    !< Soil temperature, K
    real, intent(in) :: rootmass(ilg,icc) !< root biomass for the 9 pfts in \f$kg c/m^2\f$
    real, intent(in) :: rmatctem(ilg,icc,ignd) !< fraction of roots in each layer for each pft

    real, intent(inout) :: roottemp(ilg,icc) !< root temperature (k)
    real, intent(inout) :: rmsveg(ilg,icc)   !< Maintenance respiration for stem for the CTEM pfts (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(inout) :: rmrveg(ilg,icc)   !< Maintenance respiration for root for the CTEM pfts (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, intent(in) :: tracerStemMass(:,:)  !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerRootMass(:,:)  !< Tracer mass in the root for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(out) :: rmsTracer(ilg,icc)  !< Tracer maintenance respiration for stem for the CTEM pfts (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rmrTracer(ilg,icc)  !< Tracer maintenance respiration for root for the CTEM pfts both (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real :: tempq10r(ilg,icc)
    real :: tempq10s(ilg)
    integer :: i, j, k
    integer :: n,m
    real :: q10
    real :: q10funcStem, q10funcRoot
    real :: livstmfr(ilg,icc)
    real :: livrotfr(ilg,icc)
    real :: tot_rmat(ilg,icc)
    logical consq10

    !---------------------------------------------------

    !> Set the following switch to .true. for using constant temperature
    !! indepedent q10 specified below
    consq10 = .false.

    !> q10 - if using a constant temperature independent value, i.e.
    !> if consq10 is set to true
    q10 = 2.00

    ! ---------------------------------------------------

    ! Initialize required arrays to zero
    roottemp = 0.0        ! root temperature
    tot_rmat = 0.0
    rmsveg = 0.0          ! stem maintenance respiration
    rmrveg = 0.0          ! root maintenance respiration
    livstmfr = 0.0         ! live stem fraction
    livrotfr = 0.0         ! live root fraction
    rmsTracer = 0.
    rmrTracer = 0.

    ! initialization ends

    !> Based on root and stem biomass, find fraction which is live.
    !! for stem this would be the sapwood to total wood ratio.
    do j = 1,ican ! loop 120
      do m = reindexPFTs(j,1),reindexPFTs(j,2) ! loop 125
        do i = il1,il2 ! loop 130
          select case (classpfts(j))
          case ('Crops', 'Grass') ! crops and grass
            livstmfr(i,m) = 1.0
            livrotfr(i,m) = 1.0
          case ('NdlTr','BdlTr','BdlSh')
            livstmfr(i,m) = exp( - 0.2835 * stemmass(i,m))  ! following century model
            livstmfr(i,m) = max(minlvfr,min(livstmfr(i,m),1.0))
            livrotfr(i,m) = exp( - 0.2835 * rootmass(i,m))
            livrotfr(i,m) = max(minlvfr,min(livrotfr(i,m),1.0))
          case default
            print * ,'Unknown CLASS PFT in mainres ',classpfts(j)
            call errorHandler('mainres', - 1)
          end select
        end do ! loop 130
      end do ! loop 125
    end do ! loop 120

    !> Fraction of roots for each vegetation type, for each soil layer,
    !! in each grid cell is given by rmatctem (grid cell, veg type, soil layer)
    !! which bio2str subroutine calculates. rmatctem can thus be used
    !! to find average root temperature for each plant functional type
    do j = 1,icc
      do  i = il1,il2
        if (fcan(i,j) > 0.) then
          do  k = 1,ignd
            if (isand(i,k) >= - 2) then
              roottemp(i,j) = roottemp(i,j) + tbar(i,k) * rmatctem(i,j,k)
              tot_rmat(i,j) = tot_rmat(i,j) + rmatctem(i,j,k)
            end if
          end do
          roottemp(i,j) = roottemp(i,j) / tot_rmat(i,j)
        end if
      end do
    end do

    !> We assume that stem temperature is same as air temperature, ta.
    !! using stem and root temperatures we can find their maintenance respirations rates
    do i = il1,il2 ! loop 200

      !> first find the q10 response function to scale base respiration
      !! rate from 15 c to current temperature, we do the stem first.
      if (.not. consq10) then
        !> when finding temperature dependent q10, use temperature which
        !! is close to average of actual temperature and the temperature
        !! at which base rate is specified
        tempq10s(i) = (15.0 + 273.16 + ta(i)) / 1.9
        q10 = 3.22 - 0.046 * (tempq10s(i) - 273.16)
        q10 = min(4.0,max(1.5,q10))
      end if

      q10funcStem = q10 ** (0.1 * (ta(i) - 288.16))

      do j = 1,icc ! loop 210
        if (fcan(i,j) > 0.) then

          !> This q10 value is then used with the base rate of respiration
          !! (commonly taken at some reference temperature (15 deg c), see Tjoelker et
          !! al. 2009 New Phytologist or Atkin et al. 2000 New Phyto for
          !! an example.). Long-term acclimation to temperature could be occuring
          !! see King et al. 2006 Nature SOM for a possible approach. JM.
          if (leapnow) then
            rmsveg(i,j) = stemmass(i,j) * livstmfr(i,j) * q10funcStem &
                          * (bsrtstem(sort(j)) / 366.0)
          else
            rmsveg(i,j) = stemmass(i,j) * livstmfr(i,j) * q10funcStem &
                          * (bsrtstem(sort(j)) / 365.0)
          end if

          if (useTracer > 0) then
            ! If our tracers are present then calculate the tracer flux values,
            ! include here the conversion from  kg C/m2.day -> u mol CO2/m2.sec
            if (leapnow) then
              rmsTracer(i,j) = tracerStemMass(i,j) * livstmfr(i,j) * q10funcStem &
                               * (bsrtstem(sort(j)) / 366.0) * 963.62
            else
              rmsTracer(i,j) = tracerStemMass(i,j) * livstmfr(i,j) * q10funcStem &
                               * (bsrtstem(sort(j)) / 365.0) * 963.62
            end if
          end if

          !> convert kg c/m2.day -> u mol co2/m2.sec
          rmsveg(i,j) = rmsveg(i,j) * 963.62

          !> root respiration
          if (.not. consq10) then
            tempq10r(i,j) = (15.0 + 273.16 + roottemp(i,j)) / 1.9
            q10 = 3.22 - 0.046 * (tempq10r(i,j) - 273.16)
            q10 = min(4.0,max(1.5,q10))
          end if

          q10funcRoot = q10 ** (0.1 * (roottemp(i,j) - 288.16))

          if (leapnow) then
            rmrveg(i,j) = rootmass(i,j) * livrotfr(i,j) * q10funcRoot &
                          * (bsrtroot(sort(j))/366.0)
          else
            rmrveg(i,j) = rootmass(i,j) * livrotfr(i,j) * q10funcRoot &
                          * (bsrtroot(sort(j))/365.0)
          end if

          !> convert kg c/m2.day -> u mol co2/m2.sec
          rmrveg(i,j) = rmrveg(i,j) * 963.62

          if (useTracer > 0) then
            ! If our tracers are present then calculate the tracer flux values
            ! include here the conversion from  kg C/m2.day -> u mol CO2/m2.sec
            if (leapnow) then
              rmrTracer(i,j) = tracerRootMass(i,j) * livrotfr(i,j) * q10funcRoot &
                               * (bsrtroot(sort(j))/366.0) * 963.62
            else
              rmrTracer(i,j) = tracerRootMass(i,j) * livrotfr(i,j) * q10funcRoot &
                               * (bsrtroot(sort(j))/365.0) * 963.62
            end if
          end if

        end if ! fcan check.
      end do ! loop 210
    end do ! loop 200

    return

  end subroutine mainres
  !! @}
  !> \ingroup autotrophic_res_growthRespiration
  !! @{
  !> Calculates growth respiration for all PFTs
  !> @author Vivek Arora and Joe Melton
  subroutine growthRespiration(il1,il2,ilg,sort,useTracer,nppveg,tracerNPP,&
                                rgveg,tracerRG)

    use classicParams,     only : icc,grescoef,zero

    implicit none

    integer, intent(in) :: il1      !< il1=1
    integer, intent(in) :: il2      !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg      !< nlat * nmos 
    integer, intent(in) :: sort(:)  !< index for correspondence between biogeochem pfts and the number of values in parameters vectors in run params file.   
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.   
    real, intent(in) :: nppveg(:,:)  !< NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: tracerNPP(:,:) !< tracer NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerRG(ilg,icc)    !< Tracer growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rgveg(ilg,icc)   !< PFT level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    
    integer :: j,i

    do j = 1,icc
      do i = il1,il2
        if (nppveg(i,j) > zero) then
          rgveg(i,j) = grescoef(sort(j)) * nppveg(i,j)
          if (useTracer > 0) then          
            tracerRG(i,j) = 0.
            if (tracerNPP(i,j) > 0.) tracerRG(i,j) = grescoef(sort(j)) * tracerNPP(i,j)
          end if 
        else
          rgveg(i,j) = 0.0
        end if
      end do 
    end do 
    
  end subroutine growthRespiration
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \namespace autotrophic_res
  !> Calculates maintenance respiration, over a given sub-area, for stem and root components.
  !! leaf respiration is estimated within the phtsyn subroutine.
  !! @author V. Arora, J. Melton
  !!
  !! Autotrophic respiration (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$) is composed of maintenance, \f$R_\mathrm{m}\f$, and growth respirations, \f$R_\mathrm{g}\f$,
  !! \f[
  !! R_\mathrm{a} =R_\mathrm{m} + R_\mathrm{g}.
  !! \f]
  !! Maintenance respiration accounts for carbon consumed by processes that keep existing plant tissues alive and is a function of environmental stresses. Maintenance respiration is calculated on a half-hourly time step (with photosynthesis) for the leaves, \f$R_{mL}\f$, and at a daily time step for the stem, \f$R_{mS}\f$, and root, \f$R_{mR}\f$, components
  !! \f[
  !! \label{mainres_all} R_\mathrm{m} = R_{mL} + R_{mS} + R_{mR}.
  !! \f]
  !!
  !! Maintenance respiration is generally strongly correlated with nitrogen content (Reich et al. 1998; Ryan 1991) \cite Reich1998-zr \cite Ryan1991-ai. The current version of CTEM does not explicitly track nitrogen in its vegetation components. Therefore, we adopt the approach of Collatz et al. (1991,1992) \cite Collatz1991-5bc \cite Collatz1992-jf in which the close relation between maximum catalytic capacity of Rubisco, \f$V_\mathrm{m}\f$, and leaf nitrogen content is used as a proxy to estimate leaf maintenance respiration,
  !! \f[
  !! R_{mL} = \varsigma_\mathrm{L}V_\mathrm{m}\, f_{25}(Q_10d, n)f_{PAR},
  !! \f]
  !! where \f$\varsigma_\mathrm{L}\f$ is set to 0.015 and 0.025 for \f$C_3\f$ and \f$C_4\f$ plants, respectively, \f$f_{PAR}\f$ scales respiration from the leaf to the canopy level, similar to Eq. (\ref{G_canopy}), and the \f$f_{25}(Q_10d, n)\f$ function accounts for different temperature sensitivities of leaf respiration during day (\f$d\f$) and night (\f$n\f$). Pons and Welschen (2003) \cite Pons2003-f26 and Xu and Baldocchi (2003) \cite Xu2003-d75 suggest lower temperature sensitivity for leaf respiration during the day compared to night, and therefore we use values of \f$Q_10d=1.3\f$ and \f$Q_10n=2.0\f$ for day and night, respectively.
  !!
  !! Maintenance respiration from the stem and root components is estimated based on PFT-specific base respiration rates (\f$\varsigma_\mathrm{S}\f$ and \f$\varsigma_\mathrm{R}\f$ specified at \f$15\, C\f$, \f$kg\, C\, (kg\, C)^{-1}\, yr^{-1}\f$; see also classicParams.f90) that are modified to account for temperature response following a \f$Q_{10}\f$ function. Maintenance respiration from stem and root components, \f$R_{m\{S, R\}}\f$, is calculated as
  !! \f[
  !! \label{r_msr} R_{\mathrm{m}, i} = 2.64 \times 10^{-6}\varsigma_{i}l_{\mathrm{v}, i}C_{i}f_{15}(Q_{10}), \quad i = \mathrm{S}, \mathrm{R},
  !! \f]
  !! where \f$l_{v, i}\f$ is the live fraction of stem or root component, i.e. the sapwood, and \f$C_i\f$ is the stem or root carbon mass (\f$kg\, C\, m^{-2}\f$). The constant \f$2.64 \times 10^{-6}\f$ converts units from \f$kg\, C\, m^{-2}\, yr^{-1}\f$ to \f$mol\, CO_2\, m^{-2}\, s^{-1}\f$. The live sapwood fraction, \f$l_{\mathrm{v}, i}\f$, for stem or root component is calculated following the CENTURY model (Parton et al. 1996) \cite Parton1996-zv as
  !! \f[
  !! l_{\mathrm{v}, i} = \max(0.05, \min[1.0, \exp^{-0.2835 C_i} ]), \quad i = \mathrm{S}, \mathrm{R}.
  !! \f]
  !!
  !! The \f$Q_{10}\f$ value used in Eq. (\ref{r_msr}) is not assumed to be constant but modelled as a function of temperature following Tjoelker et al. (2001) \cite Tjoelker2001-uz as
  !! \f[
  !! Q_{10} = 3.22 - 0.046 \left(\frac{15.0 + T_{\{S, R\}}}{1.9}\right),
  !! \f]
  !! where \f$T_{\{S, R\}}\f$ is stem or root temperature (\f$C\f$). Stem temperature is assumed to be the same as air temperature while root temperature is based on the soil temperature weighted by the fraction of roots present in each soil layer (Arora and Boer, 2003) \cite Arora2003838. The calculated \f$Q_{10}\f$ value is additionally constrained to be between 1.5 and 4.0.
  !!
  !! Growth respiration, \f$R_\mathrm{g}\f$ (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), is estimated as a fraction (\f$\epsilon_\mathrm{g}=0.15\f$) of the positive gross canopy photosynthetic rate after maintenance respiration has been accounted for
  !! \f[
  !! \label{growth_res} R_\mathrm{g}=\epsilon_\mathrm{g}\max[0, (G_{canopy} - R_\mathrm{m})].
  !! \f]
  !! Finally, net primary productivity (\f$NPP\f$) is calculated as
  !! \f[
  !! NPP = G_{canopy} - R_\mathrm{m} - R_\mathrm{g}.
  !! \f]
  !!
end module autotrophicRespiration
