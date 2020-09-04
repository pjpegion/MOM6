!> Provides the ocean stochastic equation of state
module MOM_stoch_eos
! This file is part of MOM6. See LICENSE.md for the license.
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_file_parser,     only : get_param, param_file_type
use MOM_random,          only : PRNG,random_2d_constructor,random_2d_norm
use MOM_time_manager,    only : time_type
!use random_numbers_mod, only : getRandomNumbers,initializeRandomNumberStream,randomNumberStream

implicit none

public MOM_stoch_eos_init
public MOM_stoch_eos_run

real, allocatable, private :: l2_inv(:,:)
real, allocatable, private :: rgauss(:,:)
real, parameter,private :: tfac=0.27
real, parameter,private :: amplitude=0.624499 !0.39 is variance
integer        ,private :: seed
type(PRNG)  ::  CS
!integer(8) :: npts
!real :: pi
logical first_time

contains
  subroutine MOM_stoch_eos_init(G,Time,param_file,random_pattern,new_run)
! initialization subroutine called my MOM.F90,
  type(param_file_type), intent(in)    :: param_file  !< structure indicating parameter file to parse
  type(ocean_grid_type), intent(in)    :: G
  type(time_type),       intent(in)    :: Time
  real,                  intent(inout) :: random_pattern(:,:)
  logical,               intent(in)    :: new_run
  integer :: i,j
  seed=0
  ! contants
  !pi=2*acos(0.0)
  allocate(l2_inv(G%isd:G%ied,G%jsd:G%jed)) ! currently allocating for data domain, should I which to compute domain?
  allocate(rgauss(G%isd:G%ied,G%jsd:G%jed)) ! currently allocating for data domain, should I which to compute domain?
  call get_param(param_file, "MOM", "SEED_STOCH_EOS", seed, &
                 "Specfied seed for random number sequence ", default=0)
  call random_2d_constructor(CS, G%HI, Time, seed)
  call random_2d_norm(CS, G%HI, rgauss)
  ! fill array with approximation of grid area needed for decorrelation
  ! time-scale calculation
  do j=G%jsd,G%jed
     do i=G%isd,G%ied
        l2_inv(i,j)=1.0/(G%dxT(i,j)**2+G%dyT(i,j)**2)
     enddo
  enddo 
  if (new_run) then
     do j=G%jsd,G%jed
        do i=G%isd,G%ied
           random_pattern(i,j)=amplitude*rgauss(i,j)
        enddo
     enddo
  endif

  
  end subroutine MOM_stoch_eos_init

  subroutine MOM_stoch_eos_run(G,u,v,random_pattern,phi_out,delt,Time)
  type(ocean_grid_type), intent(in)    :: G
  real,                  intent(in)    :: u(:,:,:)  ! zonal current
  real,                  intent(in)    :: v(:,:,:)  ! meridional current
  real,                  intent(inout) :: random_pattern(:,:)  
  real,                  intent(inout) :: phi_out(:,:)     ! autocorrelation output for diagnoistic purposes
  real,                  intent(in)    :: delt 
  type(time_type),       intent(in)    :: Time
! locals
  integer                                ::  i,j
  integer :: yr,mo,dy,hr,mn,sc
  real                                   :: phi,ubar,vbar

  call random_2d_constructor(CS, G%HI, Time, seed)
  call random_2d_norm(CS, G%HI, rgauss)
  ! advance AR(1)
  do j=G%jsd,G%jed
     do i=G%isd,G%ied
        ubar=0.5*(u(i,j,1)+u(i+1,j,1))
        vbar=0.5*(v(i,j,1)+v(i,j+1,1))
        phi=exp(-1*delt*tfac*sqrt((ubar**2+vbar**2)*l2_inv(i,j)))
        random_pattern(i,j)=phi*random_pattern(i,j) + amplitude*sqrt(1-phi**2)*rgauss(i,j)
        phi_out(i,j)=phi
     enddo
  enddo
  !print*,'stoch_run ubar',minval(u),maxval(u),minval(v),maxval(v)
  !print*,'stoch_run rp',minval(random_pattern),maxval(random_pattern),minval(rgauss),maxval(rgauss)
  !print*,'stoch_run phi',minval(phi_out),maxval(phi_out),minval(rgauss),maxval(rgauss)
  !print*,'indicies',G%isd,G%ied,G%idg_offset,G%jsd,G%jed,G%jdg_offset
  !print*,'rgauss bounds',lbound(rgauss,1),lbound(rgauss,2),ubound(rgauss,1),ubound(rgauss,2)
  !print*,'rp bounds',lbound(random_pattern,1),lbound(random_pattern,2),ubound(random_pattern,1),ubound(random_pattern,2)
  !print*,'subset',minval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd)),&
  !                maxval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd))
  
  end subroutine MOM_stoch_eos_run

!  subroutine random_gauss(rns,harvest)
!! use box-muller transform to calculate gaussian values from two uniform
!! distributions
!  implicit none
!  real,intent(inout):: harvest(npts)
!  type(randomNumberStream),intent(inout)  ::  rns
!  integer(8)            ::nxe
!  real:: random_arr(npts)
!  real r2(2)
!  call getRandomNumbers(rns,random_arr)
!  nxe=(npts*2)/2
!  harvest(1:nxe/2-1)=sqrt(-2*log(random_arr(1:nxe/2-1)))*cos(2*pi*random_arr(nxe/2:nxe))
!  harvest(nxe/2:nxe)=sqrt(-2*log(random_arr(nxe/2:nxe)))*cos(2*pi*random_arr(1:nxe/2-1))
!  if (npts.GT.nxe) then ! odd grid size
!      call getRandomNumbers(rns,r2)
!      harvest(npts)=sqrt(-2*log(r2(1)))*cos(2*pi*r2(2))
!  endif
!  end subroutine random_gauss

end module MOM_stoch_eos

