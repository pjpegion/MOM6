!> Provides the ocean stochastic equation of state
module MOM_stoch_eos
! This file is part of MOM6. See LICENSE.md for the license.
use MOM_grid,                  only : ocean_grid_type
use MOM_file_parser,          only : get_param, param_file_type
use MOM_domains, only : root_PE
use mpp_mod, only : mpp_pe,mpp_root_pe,mpp_broadcast
use random_numbers_mod, only : getRandomNumbers,initializeRandomNumberStream,randomNumberStream

implicit none

public MOM_stoch_eos_init
public MOM_stoch_eos_run

real, allocatable, private :: l2_inv(:,:)
real, allocatable, private :: rgauss(:,:)
logical                    :: global_index_logic,first_time
real, parameter,private :: tfac=0.27
real, parameter,private :: amplitude=0.624499 !0.39 is variance
type(randomNumberStream)  ::  rns
integer(8) :: npts
real :: pi

contains
  subroutine MOM_stoch_eos_init(G,param_file,global_indexing)
! initialization subroutine called my MOM.F90,
  type(param_file_type),   intent(in)   :: param_file  !< structure indicating parameter file to parse
  type(ocean_grid_type),   intent(in) :: G
  logical,                 intent(in)  :: global_indexing
  integer :: seed
  integer :: i,j
  integer(8) :: count, count_rate, count_max, count_trunc
  integer :: count4
  integer(8) :: iscale = 10000000000
  ! contants
  pi=2*acos(0.0)
  allocate(l2_inv(G%isd:G%ied,G%jsd:G%jed)) ! currently allocating for data domain, should I which to compute domain?
  allocate(rgauss(1-G%Domain%nihalo:G%Domain%niglobal+G%Domain%nihalo,1-G%Domain%njhalo:G%Domain%njglobal+G%Domain%njhalo))
  global_index_logic=global_indexing
  first_time=.true.
  call get_param(param_file, "MOM", "SEED_STOCH_EOS", seed, &
                 "Specfied seed for random number sequence ", default=0)
  ! set random number seed on root task
  if (mpp_pe()==mpp_root_pe()) then
     if (seed == 0) then
       ! generate a random seed from system clock and ens member number
       call system_clock(count, count_rate, count_max)
       ! seed is elapsed time since unix epoch began (secs)
       ! truncate to 4 byte integer
       count_trunc = iscale*(count/iscale)
       count4 = count - count_trunc !+ member_id
       print *,'using seed',count4
     else
       ! don't rely on compiler to truncate integer(8) to integer(4) on
       ! overflow, do wrap around explicitly.
       !count4 = mod(seed + 2147483648, 4294967296) - 2147483648
       count4 = seed
       print *,'using seed',count4,seed!,member_id
     endif
  endif
  call mpp_broadcast(count4, root_PE())
  rns=initializeRandomNumberStream(count4)
  npts=(G%Domain%niglobal+G%Domain%nihalo*2)*(G%Domain%njglobal+G%Domain%njhalo*2)
  ! fill array with approximation of grid area needed for decorrelation
  ! time-scale calculation
  if (global_index_logic) then
     do j=G%jsd,G%jed
        do i=G%isd,G%ied
           l2_inv(i,j)=1.0/(G%dxT(i,j)**2+G%dyT(i,j)**2)
        enddo
     enddo
  else
     do j=G%jsd,G%jed
        do i=G%isd,G%ied
           l2_inv(i,j)=1.0/(G%dxT(i,j)**2+G%dyT(i,j)**2)
        enddo
     enddo
  endif

  
  end subroutine MOM_stoch_eos_init

  subroutine MOM_stoch_eos_run(G,u,v,random_pattern,phi_out,delt)
  type(ocean_grid_type),   intent(in)    :: G
  real,                    intent(in)    :: u(:,:,:)  ! zonal current
  real,                    intent(in)    :: v(:,:,:)  ! meridional current
  real,                    intent(inout) :: random_pattern(:,:)  
  real,                    intent(inout) :: phi_out(:,:)     ! autocorrelation output for diagnoistic purposes
  real,                    intent(in)    :: delt 
! locals
  integer                                ::  i,j
  real                                   :: phi,ubar,vbar

  call random_gauss(rns,rgauss)  ! get random white noise on global grid so each PE is marching through the same random sequence for reproducibility
  if (global_index_logic) then
     if (first_time) then  ! populate pattern with random white noise
        do j=G%jsd,G%jed
           do i=G%isd,G%ied
              random_pattern(i,j)=amplitude*rgauss(i,j)
           enddo
        enddo
     else  ! advance AR(1)
        do j=G%jsd,G%jed
           do i=G%isd,G%ied
              ubar=0.5*(u(i,j,1)+u(i+1,j,1))
              vbar=0.5*(v(i,j,1)+v(i,j+1,1))
              phi=exp(-1*delt*tfac*sqrt((ubar**2+vbar**2)*l2_inv(i,j)))
              random_pattern(i,j)=phi*random_pattern(i,j) + amplitude*sqrt(1-phi**2)*rgauss(i,j)
              phi_out(i,j)=phi
           enddo
        enddo
     endif
  else  
     if (first_time) then ! populate pattern with random white noise
        do j=G%jsd,G%jed
           do i=G%isd,G%ied
              random_pattern(i,j)=amplitude*rgauss(i+G%idg_offset,j+G%jdg_offset)
           enddo
        enddo
     else  ! advance AR(1)
        do j=G%jsd,G%jed
           do i=G%isd,G%ied
              ubar=0.5*(u(i,j,1)+u(i+1,j,1))
              vbar=0.5*(v(i,j,1)+v(i,j+1,1))
              phi=exp(-1*delt*tfac*sqrt((ubar**2+vbar**2)*l2_inv(i,j)))
              random_pattern(i,j)=phi*random_pattern(i,j) + amplitude*sqrt(1-phi**2)*rgauss(i+G%idg_offset,j+G%jdg_offset)
              phi_out(i,j)=phi
           enddo
        enddo
     endif
  endif
  first_time=.false.
  !print*,'stoch_run ubar',minval(u),maxval(u),minval(v),maxval(v)
  !print*,'stoch_run rp',minval(random_pattern),maxval(random_pattern),minval(rgauss),maxval(rgauss)
  !print*,'stoch_run phi',minval(phi_out),maxval(phi_out),minval(rgauss),maxval(rgauss)
  !print*,'indicies',G%isd,G%ied,G%idg_offset,G%jsd,G%jed,G%jdg_offset
  !print*,'rgauss bounds',lbound(rgauss,1),lbound(rgauss,2),ubound(rgauss,1),ubound(rgauss,2)
  !print*,'rp bounds',lbound(random_pattern,1),lbound(random_pattern,2),ubound(random_pattern,1),ubound(random_pattern,2)
  !print*,'subset',minval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd)),&
  !                maxval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd))
  
  end subroutine MOM_stoch_eos_run

  subroutine random_gauss(rns,harvest)
! use box-muller transform to calculate gaussian values from two uniform
! distributions
  implicit none
  real,intent(inout):: harvest(npts)
  type(randomNumberStream),intent(inout)  ::  rns
  integer(8)            ::nxe
  real:: random_arr(npts)
  real r2(2)
  call getRandomNumbers(rns,random_arr)
  nxe=(npts*2)/2
  harvest(1:nxe/2-1)=sqrt(-2*log(random_arr(1:nxe/2-1)))*cos(2*pi*random_arr(nxe/2:nxe))
  harvest(nxe/2:nxe)=sqrt(-2*log(random_arr(nxe/2:nxe)))*cos(2*pi*random_arr(1:nxe/2-1))
  if (npts.GT.nxe) then ! odd grid size
      call getRandomNumbers(rns,r2)
      harvest(npts)=sqrt(-2*log(r2(1)))*cos(2*pi*r2(2))
  endif
  end subroutine random_gauss

end module MOM_stoch_eos

