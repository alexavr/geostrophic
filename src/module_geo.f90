module module_geo
use module_globals
use module_io, only:  get_nc_var, wrf_user_unstagger, put_var_time, finalize, &
                      time_start, time_end

implicit none

private

public :: run_computations 

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine run_computations ! potentially parallel
  integer :: itime, nb_ticks_run_computations 
  real(kind=real32) :: elapsed_time_run_computations

  do itime = 1, tDimSIZE

    if (debug) call time_start(nb_ticks_run_computations)
    
    if(iswrf) then

      if(scheme.eq.0) then
        
        call run_computations_wrf_simple(itime)
      
      else if(scheme.eq.1) then
        
        call run_computations_wrf_pert(itime)
      
      else

        print*,"UNKNOW SCHEME OPTIONS. STOP."
        call finalize

      end if
    
    else
      
      call run_computations_rean_simple(itime)
    
    end if

    if (debug) then
      call time_end(nb_ticks_run_computations,elapsed_time_run_computations)
      print'(a,i3,a,i5,a,f7.2,a)',"run_computations: time step ",itime," of ", tDimSIZE, " took ", elapsed_time_run_computations, " sec"
    end if


  end do


end subroutine run_computations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine run_computations_wrf_pert(itime) ! potentially parallel
  integer, intent(in) :: itime 
  real(kind=real32)   :: RDX, RDY
  real(kind=real32),dimension(:)    ,allocatable :: FNM, FNP, RDNW, RDN
  real(kind=real32),dimension(:,:)  ,allocatable :: MUB, MU, MUT
  real(kind=real32),dimension(:,:,:),allocatable :: PT, P, PB
  real(kind=real32),dimension(:,:,:),allocatable :: PH, PHB, PHT
  real(kind=real32),dimension(:,:,:),allocatable :: tmp3d1,tmp3d2,tmp3d3
  real(kind=real32),dimension(:,:,:),allocatable :: dPHTdnu, dPHdnu
  real(kind=real32),dimension(:,:,:),allocatable :: alp_d, alt_d, alt_m
  real(kind=real32),dimension(:,:,:),allocatable :: DPDX, DPDY
  real(kind=real32),dimension(:,:)  ,allocatable :: F, MAPFAC_MX, MAPFAC_MY
  real(kind=real32),dimension(:,:,:),allocatable :: U,V,Ug,Vg,Ua,Va
  real(kind=real32),dimension(:,:,:),allocatable :: W,Wg,Wa


  integer :: dimIDs_m(4),dimIDs_u(4),dimIDs_v(4)
  integer :: ii, jj, kk

  dimIDs_m = (/ xDimIDout, yDimIDout, zDimIDout, tDimIDout /)
  dimIDs_u = (/ xDimIDout_stag, yDimIDout     , zDimIDout, tDimIDout /)
  dimIDs_v = (/ xDimIDout     , yDimIDout_stag, zDimIDout, tDimIDout /)

  call get_nc_var(ncid_in, "RDX", itime, RDX)
  call get_nc_var(ncid_in, "RDY", itime, RDY)

  allocate( FNM (zDimSIZE) )
  allocate( FNP (zDimSIZE) )
  allocate( RDNW(zDimSIZE) )
  allocate(  RDN(zDimSIZE) )

  call get_nc_var(ncid_in, "FNM", itime, FNM)   ! upper weight for vertical stretching
  call get_nc_var(ncid_in, "FNP", itime, FNP)   ! lower weight for vertical stretching
  call get_nc_var(ncid_in, "RDNW",itime, RDNW)  ! inverse d(eta) values between full (w) levels
  call get_nc_var(ncid_in, "RDN" ,itime, RDN)   ! inverse d(eta) values between half (mass) levels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE dry air mass
  allocate( MU  ( xDimSIZE, yDimSIZE ) )
  allocate( MUB ( xDimSIZE, yDimSIZE ) )
  allocate( MUT ( xDimSIZE, yDimSIZE ) )
  call get_nc_var(ncid_in, "MU" , itime, MU)   ! perturbation dry air mass in column [Pa]
  call get_nc_var(ncid_in, "MUB", itime, MUB)  ! base state dry air mass in column [Pa]
  MUT = MU + MUB                               ! total dry air mass in column [Pa]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE TOTAL PRESSURE
  allocate( P  ( xDimSIZE, yDimSIZE, zDimSIZE ) ) 
  allocate( PB ( xDimSIZE, yDimSIZE, zDimSIZE ) ) 

  call get_nc_var(ncid_in, "P",  itime, P )
  call get_nc_var(ncid_in, "PB", itime, PB)

  allocate( PT ( xDimSIZE, yDimSIZE, zDimSIZE ) ) 
  PT = P + PB         ! Total pressure [Pa]
  call put_var_time(ncid_out,PT,"PT","Total pressure","Pa", itime, dimIDs_m )
  deallocate( PT )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE GEOPOTENTIAL (originally on w-grid, but we need to save on m-grid)
  allocate( PH     ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  allocate( PHB    ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  allocate( PHT    ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  allocate( tmp3d1 ( xDimSIZE, yDimSIZE, zDimSIZE ) )

  call get_nc_var(ncid_in, "PH",  itime, PH )
  call get_nc_var(ncid_in, "PHB", itime, PHB)

  PHT = PH + PHB                          ! Total geopotential [Pa] on w-grid
  
  call wrf_user_unstagger(PHT,tmp3d1)     ! Ustagering for output only
  
  call put_var_time(ncid_out, tmp3d1,"PHT","Total geopotential","m2 s-2",itime, dimIDs_m )        
  
  tmp3d1 = tmp3d1/g
  ! call put_var_time(ncid_out, tmp3d1,"Z","Total geopotential surface hight","m",itime, dimIDs_m ) 
  deallocate(tmp3d1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! VERTICAL GRADIENTS FOR PHI AND PHI'
  allocate( dPHTdnu( xDimSIZE,yDimSIZE,zDimSIZE ) )  
  allocate(  dPHdnu( xDimSIZE,yDimSIZE,zDimSIZE ) )  

  ! vertical gradients of phi and phi'
  do jj = 1, yDimSIZE
  do ii = 1, xDimSIZE
  do kk = 1, zDimSIZE
      dPHTdnu(ii,jj,kk) = ( PHT(ii,jj,kk+1) - PHT(ii,jj,kk) ) * RDNW(kk)
       dPHdnu(ii,jj,kk) = (  PH(ii,jj,kk+1) -  PH(ii,jj,kk) ) * RDNW(kk)
  end do
  end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE DRY INVERSE DENSITY (ALPHA)
! Using diagnostic relation (eq 28 and 32)
  allocate( alt_d( xDimSIZE,yDimSIZE,zDimSIZE ) ) ! total dry air inverse density
  allocate( alp_d( xDimSIZE,yDimSIZE,zDimSIZE ) ) ! perturbation

  do jj = 1, yDimSIZE
  do ii = 1, xDimSIZE
  do kk = 1, zDimSIZE
      alt_d(ii,jj,kk) = (-1.)*dPHTdnu(ii,jj,kk)/MUT(ii,jj)
      alp_d(ii,jj,kk) = (-1.)*(dPHdnu(ii,jj,kk)+alt_d(ii,jj,kk)*MU(ii,jj))/MUB(ii,jj)
  end do
  end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE WET DENSITY
  allocate( tmp3d1 ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( tmp3d2 ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( tmp3d3 ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( alt_m  ( xDimSIZE,yDimSIZE,zDimSIZE ) ) ! total wet air inverse density 
  call get_nc_var(ncid_in, "QVAPOR", itime, tmp3d1)
  call get_nc_var(ncid_in, "QCLOUD", itime, tmp3d2)
  tmp3d1 = tmp3d1 + tmp3d2
  call get_nc_var(ncid_in, "QRAIN",  itime, tmp3d2)
  tmp3d3 = tmp3d1 + tmp3d2
  
  ! We used WSM3 mp scheme in coarse run, so no QICE and QGRAUP in output
  if( (1000./RDX) .LT. 30 ) then
      call get_nc_var(ncid_in, "QICE",   itime, tmp3d2)
      tmp3d1 = tmp3d3 + tmp3d2
      call get_nc_var(ncid_in, "QSNOW",  itime, tmp3d2)
      tmp3d1 = tmp3d1 + tmp3d2
      call get_nc_var(ncid_in, "QGRAUP", itime, tmp3d2)
      tmp3d1 = tmp3d1 + tmp3d2
      alt_m = alt_d*(1.+tmp3d1)**(-1)   ! From WRF Tech Note (p12)
      deallocate(tmp3d1, tmp3d2)
  else
      alt_m = alt_d*( 1.+tmp3d3 )**(-1) ! From WRF Tech Note (p12)
      deallocate(tmp3d1, tmp3d2, tmp3d3 )
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! RETRIVE HORIZONTAL PRESSURE GRADIENT
! The idea is to use different appoach: full variables, pertrubations with no 
! analitical asumptions and wrf formulation
  allocate( DPDX( xDimSIZE_stag,yDimSIZE     ,zDimSIZE ) )
  allocate( DPDY( xDimSIZE     ,yDimSIZE_stag,zDimSIZE ) )

  call horizontal_pressure_gradient_pertrubations_wrf(P,PB,PH,PHT,PHB, &
                                                    MU,MUT,alt_d,alt_m,alp_d,    &
                                                    FNM,FNP,RDX,RDY,RDNW, &
                                                    DPDX, DPDY,itime)

  ! call put_var_time(ncid_out,DPDX ,"DPDX ","DPDX","kg m-2 s-2",itime,dimIDs_u)
  ! call put_var_time(ncid_out,DPDY ,"DPDY ","DPDY","kg m-2 s-2",itime,dimIDs_v)

  deallocate( P, PB )
  deallocate( PHT, PH, PHB )
  deallocate( dPHTdnu, dPHdnu )
  deallocate( alp_d, alt_d, alt_m )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! GEOSTROPHIC BALANCE
! VERCION ON MASS GRID POINTS
  allocate( Ug ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( Vg ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( F  ( xDimSIZE,yDimSIZE          ) )
  allocate( MAPFAC_MX ( xDimSIZE,yDimSIZE ) )
  allocate( MAPFAC_MY ( xDimSIZE,yDimSIZE ) )
  Ug = fFillValue
  Vg = fFillValue

  call get_nc_var(ncid_in, "F" , itime, F)                                ! Coriolis sine latitude term [2*omega*sin(lat)]
  call get_nc_var(ncid_in, "MAPFAC_MX" , itime, MAPFAC_MX)                ! map scale factors: 
  call get_nc_var(ncid_in, "MAPFAC_MY" , itime, MAPFAC_MY)                       
  ! call get_nc_var(ncid_in, "MAPFAC_VX" , itime, MAPFAC_VX)                       
  ! call get_nc_var(ncid_in, "MAPFAC_VY" , itime, MAPFAC_VY)                       

  do kk = 1, zDimSIZE
  do ii = 1, xDimSIZE
  do jj = 1, yDimSIZE
      if( DPDX(ii,jj,kk).NE.fFillValue .AND. DPDX(ii,jj,kk).NE.fFillValue ) then
          Vg(ii,jj,kk) =   0.5 * 1./F(ii,jj) * MAPFAC_MX(ii,jj)/MAPFAC_MY(ii,jj)*MAPFAC_MX(ii,jj)/MUT(ii,jj)*(DPDX(ii+1,jj  ,kk)+DPDX(ii,jj,kk))
          Ug(ii,jj,kk) =  -0.5 * 1./F(ii,jj) * MAPFAC_MY(ii,jj)/MAPFAC_MX(ii,jj)*MAPFAC_MY(ii,jj)/MUT(ii,jj)*(DPDY(ii  ,jj+1,kk)+DPDY(ii,jj,kk))
      end if
  end do
  end do
  end do

  ! where(Ug.NE.missing .AND. Vg.NE.missing) Wg = sqrt(Ug**2+Vg**2)
  ! write(*,'(10x, "WG overspeed: ", i6 , " Max value at 2 level: ", f7.1)') count(Wg.GT.100), maxval(Wg(:,:,2)) 

  deallocate( DPDX, DPDY )
  deallocate( F, MAPFAC_MX, MAPFAC_MY )
  deallocate( MUT, MU, MUB )
  deallocate( FNM, FNP, RDNW, RDN)

  ! AGEOSTROPHIC PORTION 
  allocate( Ua ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( Va ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  Ua = fFillValue

  allocate( U  ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( V  ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  Va = fFillValue

  allocate( tmp3d1 ( xDimSIZE_stag,yDimSIZE     ,zDimSIZE ) )
  call get_nc_var(ncid_in, "U", itime, tmp3d1)
  call wrf_user_unstagger(tmp3d1,U)      
  deallocate( tmp3d1 )

  allocate( tmp3d1 ( xDimSIZE     ,yDimSIZE_stag,zDimSIZE ) )
  call get_nc_var(ncid_in, "V", itime, tmp3d1)
  call wrf_user_unstagger(tmp3d1,V)              
  deallocate( tmp3d1 )

  do kk = 1, zDimSIZE
  do ii = 1, xDimSIZE
  do jj = 1, yDimSIZE
      if( Ug(ii,jj,kk).NE.fFillValue .AND. Vg(ii,jj,kk).NE.fFillValue ) then
          Ua(ii,jj,kk) = U(ii,jj,kk) - Ug(ii,jj,kk)
          Va(ii,jj,kk) = V(ii,jj,kk) - Vg(ii,jj,kk)
      end if
  end do
  end do
  end do

  call put_var_time(ncid_out,Ug,"Ug","U-component of geostrophic wind","m s-1" ,itime,dimIDs_m)
  call put_var_time(ncid_out,Vg,"Vg","V-component of geostrophic wind","m s-1" ,itime,dimIDs_m)
  call put_var_time(ncid_out,Ua,"Ua","U-component of ageostrophic wind","m s-1",itime,dimIDs_m)
  call put_var_time(ncid_out,Va,"Va","V-component of ageostrophic wind","m s-1",itime,dimIDs_m)

  ! Fast diag
    allocate( W  ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    allocate( Wg ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    allocate( Wa ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    W  = fFillValue
    Wg = fFillValue
    Wa = fFillValue
    W  = sqrt ( U**2  + V**2  )
    Wg = sqrt ( Ug**2 + Vg**2 )
    Wa = sqrt ( Ua**2 + Va**2 )
    call put_var_time(ncid_out,W  ,"W" ,"W" ,"m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,Wg ,"WG","WG","m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,Wa ,"WA","WA","m s-1",itime,dimIDs_m)
    deallocate( W,Wg,Wa )
  ! end fast diag


  deallocate( Ua, Va )
  deallocate( Ug, Vg )
  deallocate( U, V )

end subroutine run_computations_wrf_pert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine run_computations_wrf_simple(itime) ! potentially parallel
  integer, intent(in) :: itime 


end subroutine run_computations_wrf_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Assume reanalyze is on the p-levels => calculates Vg using geopotential (PHT)
! (no rho needed in this case). 
! Vg = 1/f * grad(PHT)
! * Forward scheme for grad(PHT) applied
!  
! GAVR: rewrite more accurate! 
subroutine run_computations_rean_simple(itime) 
  integer, intent(in) :: itime 
  real(kind=real32),dimension(:,:,:),allocatable :: U, V, Ug, Vg, Ua, Va, Z, PHT
  real(kind=real32),dimension(:,:,:),allocatable :: W, Wg, Wa
  real(kind=real32),dimension(:,:,:),allocatable :: tmp3d
  integer :: dimIDs_m(4)
  integer :: ii, jj, iip1, kk
  real(kind=real32) :: f, rf, rad, dx, dy
  real(kind=real32) :: dlon, dlat, degree
  real(kind=real32) :: omega = 7.292E-5

  dimIDs_m = (/ xDimIDout, yDimIDout, zDimIDout, tDimIDout /)

  allocate( PHT ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  call get_nc_var(ncid_in, "z", itime, PHT)
  ! tmp3d = PHT/g
  ! call put_var_time(ncid_out,tmp3d,"Z","Z","m" ,itime, dimIDs_m)
  ! deallocate(tmp3d)

  allocate( V ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( U ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  call get_nc_var(ncid_in, "v", itime, V)
  call get_nc_var(ncid_in, "u", itime, U)

  allocate( Vg ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( Va ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( Ug ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  allocate( Ua ( xDimSIZE,yDimSIZE,zDimSIZE ) )
  Vg = fFillValue
  Va = fFillValue
  Ug = fFillValue
  Ua = fFillValue

  dlat = abs(lat1d(2)-lat1d(1))
  dlon = abs(lon1d(2)-lon1d(1))

  do kk = 1, zDimSIZE
  do jj = 1, yDimSIZE-1

    rad = lat1d(jj)*pi/180.
    f = 2. * omega * sin( rad ) ! [ f == 2 * omega * sin(lat) ]

    rf = 1./f

    degree = 111300.*cos( rad )
    dy = degree * dlat

    do ii = 1, xDimSIZE
      
      iip1 = ii + 1
      if(ii.eq.xDimSIZE) iip1 = 1

      dx = degree * dlon

      Vg(ii,jj,kk) =   rf * (PHT(iip1,jj,kk) - PHT(ii,jj,kk) )/dx
      Ug(ii,jj,kk) = - rf * (PHT(ii,jj+1,kk) - PHT(ii,jj,kk) )/dy

    end do
  
  end do
  end do

  do kk = 1, zDimSIZE
  do ii = 1, xDimSIZE
  do jj = 1, yDimSIZE
      if( Ug(ii,jj,kk).NE.fFillValue .AND. Vg(ii,jj,kk).NE.fFillValue ) then
          Ua(ii,jj,kk) = U(ii,jj,kk) - Ug(ii,jj,kk)
          Va(ii,jj,kk) = V(ii,jj,kk) - Vg(ii,jj,kk)
      end if
  end do
  end do
  end do

  call put_var_time(ncid_out,Ug,"Ug","U-component of geostrophic wind","m s-1" ,itime,dimIDs_m)
  call put_var_time(ncid_out,Vg,"Vg","V-component of geostrophic wind","m s-1" ,itime,dimIDs_m)
  call put_var_time(ncid_out,Ua,"Ua","U-component of ageostrophic wind","m s-1",itime,dimIDs_m)
  call put_var_time(ncid_out,Va,"Va","V-component of ageostrophic wind","m s-1",itime,dimIDs_m)

  ! Fast diag
    allocate( W  ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    allocate( Wg ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    allocate( Wa ( xDimSIZE,yDimSIZE,zDimSIZE ) )
    W  = fFillValue
    Wg = fFillValue
    Wa = fFillValue
    W  = sqrt ( U**2  + V**2  )
    Wg = sqrt ( Ug**2 + Vg**2 )
    Wa = sqrt ( Ua**2 + Va**2 )
    call put_var_time(ncid_out,W  ,"W" ,"W" ,"m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,Wg ,"WG","WG","m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,Wa ,"WA","WA","m s-1",itime,dimIDs_m)
    deallocate( W,Wg,Wa )
  ! end fast diag

  deallocate(PHT)
  deallocate(U, V)
  deallocate(Ug, Vg, Ua, Va)

end subroutine run_computations_rean_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine horizontal_pressure_gradient_pertrubations_wrf(P,PB,PH,PHT,PHB,            &
                                                          MU,MUT,alt_d,alt_m,alp_d,   &
                                                          FNM,FNP,RDX,RDY,RDNW,       &
                                                          gradPdx, gradPdy,itime)
implicit none
real(kind=real32),intent(in)  :: P(:,:,:),PB(:,:,:),PH(:,:,:),PHT(:,:,:),PHB(:,:,:),  &
                      alt_d(:,:,:),alt_m(:,:,:),alp_d(:,:,:),         &
                      MU(:,:),MUT(:,:),FNM(:),FNP(:), RDNW(:),        &
                      RDX,RDY
integer,intent(in)  :: itime
real(kind=real32),intent(out)  :: gradPdx(:,:,:), gradPdy(:,:,:)
integer              :: ii, jj, kk
real(kind=real32)    :: gMU, gMUT, gALTd, gALTm, gALPd
real(kind=real32)    :: gdPHTdx, gdPHdx, gdPBdx, gdPdx, gdPdnu
real(kind=real32)    :: gdPHTdy, gdPHdy, gdPBdy, gdPdy
real(kind=real32)    :: gdPHBdx, gdPHBdy 
real(kind=real32)    :: cfn, cfn1, cf1, cf2, cf3, gdPnum,gdPnup
logical     :: verbose = .false.

gradPdx = fFillValue
gradPdy = fFillValue

call get_nc_var(ncid_in, "CFN",  itime, cfn)                            ! extrapolation constants
call get_nc_var(ncid_in, "CFN1", itime, cfn1)                           ! extrapolation constants
call get_nc_var(ncid_in, "CF1",  itime, cf1)                            ! extrapolation constants
call get_nc_var(ncid_in, "CF2",  itime, cf2)                            ! extrapolation constants
call get_nc_var(ncid_in, "CF3",  itime, cf3)                            ! extrapolation constants

! start with the east-west (x) pressure gradient
if(verbose) then
  write(*,'(" PDPX....")')
  write(*,'("    1      |     2      |     3      |     4      |     5      |     6      |     7      |     8      |")')
end if

    do jj = 1, yDimSIZE
    do ii = 2, xDimSIZE
    do kk = 1, zDimSIZE
      gMU  =  0.5*(   MU (ii,jj)    +   MU (ii-1,jj)    )
      gMUT =  0.5*(   MUT(ii,jj)    +   MUT(ii-1,jj)    )
      gALTd = 0.5*( alt_d(ii,jj,kk) + alt_d(ii-1,jj,kk) )
      gALTm = 0.5*( alt_m(ii,jj,kk) + alt_m(ii-1,jj,kk) )
      gALPd = 0.5*( alp_d(ii,jj,kk) + alp_d(ii-1,jj,kk) )
      gdPHTdx = 0.5*RDX*( PHT(ii,jj,kk+1) + PHT(ii,jj,kk) - PHT(ii-1,jj,kk+1) - PHT(ii-1,jj,kk) )
      gdPHBdx = 0.5*RDX*( PHB(ii,jj,kk+1) + PHB(ii,jj,kk) - PHB(ii-1,jj,kk+1) - PHB(ii-1,jj,kk) )
      gdPHdx  = 0.5*RDX*(  PH(ii,jj,kk+1) +  PH(ii,jj,kk) -  PH(ii-1,jj,kk+1) -  PH(ii-1,jj,kk) )
      gdPdx  = RDX*(  P(ii,jj,kk) -  P(ii-1,jj,kk) )
      gdPBdx = RDX*( PB(ii,jj,kk) - PB(ii-1,jj,kk) )

      if(kk.EQ.1) then
        gdPnum = 0.5*( cf1*(P(ii-1,jj,kk  )+P(ii,jj,kk  ))      &
                      +cf2*(P(ii-1,jj,kk+1)+P(ii,jj,kk+1))      &
                      +cf3*(P(ii-1,jj,kk+2)+P(ii,jj,kk+2))  )
        gdPnup = 0.5*( FNM(2)*(P(ii,jj,2  )+P(ii-1,jj,2  ))     & 
                      +FNP(2)*(P(ii,jj,2-1)+P(ii-1,jj,2-1)) )   

      elseif(kk.EQ.zDimSIZE) then
        gdPnup = 0.5*( cfn *(P(ii-1,jj,kk-1)+P(ii,jj,kk-1))     &
                      +cfn1*(P(ii-1,jj,kk-2)+P(ii,jj,kk-2)) )
        gdPnum = 0.5*( FNM(kk-1)*(P(ii,jj,kk-1  )+P(ii-1,jj,kk-1  ))     & 
                      +FNP(kk-1)*(P(ii,jj,kk-1-1)+P(ii-1,jj,kk-1-1)) )     

      else
        gdPnup = 0.5*( FNM(kk+1)*(P(ii,jj,kk+1  )+P(ii-1,jj,kk+1  ))     &  ! 1. mass -> u-grid -> w-grid
                      +FNP(kk+1)*(P(ii,jj,kk+1-1)+P(ii-1,jj,kk+1-1)) )      ! 
        gdPnum = 0.5*( FNM(kk  )*(P(ii,jj,kk    )+P(ii-1,jj,kk    ))     &  ! 1. mass -> u-grid -> w-grid     
                      +FNP(kk  )*(P(ii,jj,kk-1  )+P(ii-1,jj,kk-1  )) )      !      
      end if

      gdPdnu = RDNW(kk) * (gdPnup-gdPnum)

      gradPdx(ii,jj,kk) = gALTm/gALTd * ( gMUT*( gdPHdx + gALTd*gdPdx + gALPd*gdPBdx ) + gdPHTdx*(gdPdnu-gMU)  )

      if(verbose) then
        if(kk.EQ.2 .AND. jj.EQ.480 .AND. MOD(ii,10).EQ.0)  write(*,'(8(f15.7," "))') &
          ! gMUT*gdPHdx, gMUT*gALTd*gdPdx, gMUT*gALPd*gdPBdx, gdPHBdx*gdPdnu, -gdPHBdx*gMU, gdPHdx*gdPdnu, -gdPHdx*gMU, gradPdx(ii,jj,kk)
          abs(gMUT*gdPHdx), abs(gMUT*gALTd*gdPdx), abs(gMUT*gALPd*gdPBdx), abs(gdPHBdx*gdPdnu), abs(gdPHBdx*gMU), abs(gdPHdx*gdPdnu), abs(gdPHdx*gMU), abs(gradPdx(ii,jj,kk))
      endif 
    end do
    end do
    end do

! now the north-south (y) pressure gradient

    if(verbose) then
      write(*,'(" DPDY....")')
      write(*,'("    1      |     2      |     3      |     4      |     5      |     6      |     7      |     8      |")')
    end if

    do jj = 2, yDimSIZE
    do ii = 1, xDimSIZE
    do kk = 1, zDimSIZE
      gMU  = 0.5 * ( MU (ii,jj)    + MU (ii,jj-1)    )
      gMUT = 0.5 * ( MUT(ii,jj)    + MUT(ii,jj-1)    )
      gALTd = 0.5 * ( alt_d(ii,jj,kk) + alt_d(ii,jj-1,kk) )
      gALTm = 0.5 * ( alt_m(ii,jj,kk) + alt_m(ii,jj-1,kk) )
      gALPd = 0.5 * ( alp_d(ii,jj,kk) + alp_d(ii,jj-1,kk) )
      gdPHTdy = 0.5*RDY*( PHT(ii,jj,kk+1) + PHT(ii,jj,kk) - PHT(ii,jj-1,kk+1) - PHT(ii,jj-1,kk) )
      gdPHBdy = 0.5*RDY*( PHB(ii,jj,kk+1) + PHB(ii,jj,kk) - PHB(ii,jj-1,kk+1) - PHB(ii,jj-1,kk) )
      gdPHdy  = 0.5*RDY*(  PH(ii,jj,kk+1) +  PH(ii,jj,kk) -  PH(ii,jj-1,kk+1) -  PH(ii,jj-1,kk) )
      gdPdy  = RDY * (  P(ii,jj,kk) -  P(ii,jj-1,kk) )
      gdPBdy = RDY * ( PB(ii,jj,kk) - PB(ii,jj-1,kk) )

      if(kk.EQ.1) then
        gdPnum = 0.5*( cf1*(P(ii,jj-1,kk  )+P(ii,jj,kk  ))   &
                      +cf2*(P(ii,jj-1,kk+1)+P(ii,jj,kk+1))   &
                      +cf3*(P(ii,jj-1,kk+2)+P(ii,jj,kk+2))  )
        gdPnup = 0.5*( FNM(2)*(P(ii,jj,2  )+P(ii,jj-1,2  ))    & 
                      +FNP(2)*(P(ii,jj,2-1)+P(ii,jj-1,2-1)) )   

      elseif(kk.EQ.zDimSIZE) then
        gdPnup = 0.5*( cfn *(P(ii,jj-1,kk-1)+P(ii,jj,kk-1))   &
                      +cfn1*(P(ii,jj-1,kk-2)+P(ii,jj,kk-2)) )
        gdPnum = 0.5*( FNM(kk-1)*(P(ii,jj,kk-1  )+P(ii,jj-1,kk-1  ))     & 
                      +FNP(kk-1)*(P(ii,jj,kk-1-1)+P(ii,jj-1,kk-1-1)) )     

      else
        gdPnup = 0.5*( FNM(kk+1)*(P(ii,jj,kk+1  )+P(ii,jj-1,kk+1  ))     & ! 1. mass -> u-grid -> w-grid
                      +FNP(kk+1)*(P(ii,jj,kk+1-1)+P(ii,jj-1,kk+1-1)) )     ! 
        gdPnum = 0.5*( FNM(kk  )*(P(ii,jj,kk    )+P(ii,jj-1,kk    ))     & ! 1. mass -> u-grid -> w-grid
                      +FNP(kk  )*(P(ii,jj,kk-1  )+P(ii,jj-1,kk-1  )) )     ! 
      end if

      gdPdnu = RDNW(kk) * (gdPnup-gdPnum)

      gradPdy(ii,jj,kk) = gALTm/gALTd * ( gMUT*( gdPHdy + gALTd*gdPdy + gALPd*gdPBdy ) + gdPHTdy*(gdPdnu-gMU)  )

    if(verbose) then
      if(kk.EQ.2 .AND. jj.EQ.480 .AND. MOD(ii,10).EQ.0) write(*,'(8(f15.7," "))') &
      ! gMUT*gdPHdy, gMUT*gALTd*gdPdy, gMUT*gALPd*gdPBdy, gdPHBdy*gdPdnu, -gdPHBdy*gMU, gdPHdy*gdPdnu, -gdPHdy*gMU, gradPdy(ii,jj,kk)
      abs(gMUT*gdPHdy), abs(gMUT*gALTd*gdPdy), abs(gMUT*gALPd*gdPBdy), abs(gdPHBdy*gdPdnu), abs(gdPHBdy*gMU), abs(gdPHdy*gdPdnu), abs(gdPHdy*gMU), abs(gradPdy(ii,jj,kk))
    end if
 
    end do
    end do
    end do

end subroutine horizontal_pressure_gradient_pertrubations_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

end module module_geo