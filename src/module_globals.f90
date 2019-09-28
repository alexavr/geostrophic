module module_globals
use,intrinsic :: iso_fortran_env,only:real32,real64
implicit none

  real(kind=real64), parameter :: dFillValue = -(HUGE(1.d0))
  real(kind=real32), parameter :: fFillValue = 1E-35
  real(kind=real32), parameter :: missing = fFillValue
  integer          , parameter :: iFillValue = -(HUGE(1))
  real(kind=real64), parameter :: g = 9.80665D0, rg = 1.D0/g
  real(kind=real32), parameter :: r_d = 287      ! [J/kg/K] gas const. dry air
  real(kind=real32), parameter :: cp = 1004.5    ! [J/kg/K] specific heat of dry air 
  real(kind=real32), parameter :: p1000mb = 100000
  real(kind=real32), parameter :: pi = 3.14159
  character(len=80), parameter :: namelist_name = "namelist.geo"

  character(len=120)  :: file_in, file_out
  integer  :: ncid_in, ncid_out
  ! integer  :: dims_w(4), dims_u(4), dims_v(4), dims_m(4)
  integer  :: ntimes

  real(kind=real32),dimension(:),allocatable :: lon1d, lat1d, levels
  real(kind=real32),dimension(:,:),allocatable :: lon2d, lat2d


  ! real(kind=real32),dimension(:,:,:),allocatable :: P, PB, PT
  ! real(kind=real32),dimension(:,:,:),allocatable :: PH, PHB, PHT
  ! real(kind=real32),dimension(:,:,:),allocatable :: alp_d, alt_d, alt_m
  ! real(kind=real32),dimension(:,:)  ,allocatable :: MU, MUB, MUT
  ! real(kind=real32),dimension(:,:)  ,allocatable :: F, MAPFAC_MX, MAPFAC_MY 
  ! real(kind=real32),dimension(:,:,:),allocatable :: Ug, Vg, Wg, Ua, Va, Wa, U, V, W 

  ! real(kind=real32),dimension(:,:,:),allocatable :: DPDX, DPDY 
  ! real(kind=real32),dimension(:)    ,allocatable :: RDNW,RDN                           ! inverse vertical steps 
  ! real(kind=real32),dimension(:)    ,allocatable :: FNM, FNP                           ! upper/lower weight for vertical stretching

  integer :: vg_id, ug_id, va_id, ua_id, v_id, u_id 
  integer :: var_dimIDs(3) 
  integer :: xDimIDout=-1, yDimIDout=-1, zDimIDout=-1, tDimIDout=-1,    &
             xDimIDout_stag=-1, yDimIDout_stag=-1, zDimIDout_stag=-1
  integer :: xDimIDin=-1, yDimIDin=-1, zDimIDin=-1, tDimIDin=-1,        &
             xDimIDin_stag=-1, yDimIDin_stag=-1, zDimIDin_stag=-1
  integer :: xDimSIZE=-1, yDimSIZE=-1, zDimSIZE=-1, tDimSIZE=-1,        &
             xDimSIZE_stag=-1, yDimSIZE_stag=-1, zDimSIZE_stag=-1

  integer           :: nb_ticks_sec, nb_ticks_max, nb_ticks_initial
  real(kind=real32) :: elapsed_time

! NAMELIST SCHEME
  logical :: iswrf=.true.
  logical :: regional=.true.
  integer :: scheme=.true.
  logical :: debug=.false.

  namelist /settings/ iswrf,regional,scheme,debug

end module module_globals
