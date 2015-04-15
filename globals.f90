!=======================================================================
! globals.f90
! YONGFEI WANG Feb, 2015
! Defined global variabls
!-----------------------------------------------------------------------


module m_globals
implicit none

real, parameter :: pi = 3.1415926, &
                   d2r = pi/180.0
                   
integer, parameter :: sizereal=4,sizeint=4

integer, parameter :: master = 0
 
real, dimension(3) :: &
	  origin
	  
integer, dimension(3) :: &
	  ihypo 
                   
real :: &
        lx,  	&
        ly,  	&
        lz,  	&
        dx,  	&
        dy,  	& 
        dz,  	& 
        dt,  	&
        t,   	&
        tt,  	&
        rho, 	&
        vp,  	&
        vs,  	&
        mu,  	&
        lam, 	&
        degint,	&
        falloff1in, &
        falloff2in, &
        df,       &
        dist,     &
        deg0,     &
        ifrq1f,   &
        ifrq2f,   &
        psfrqf(4),   &
        cfc,      &
        fc1,      &
        fc2
        
        
integer :: &
	  nst,  	&
	  nx,		&
	  ny, 	&
	  nz,		&
	  nt,		&
	  ntt,	&
	  nfreq,	&
	  ifrq1,	&
	  ifrq2,	&
	  getfc,    &
	  soff,     &
	  yorp      ! yongfei or peter

integer :: &
	  myid, numprocs, intervsta, stnum
	  
	  
real, allocatable, target, dimension(:,:,:,:) :: &
	  timeseries,     &
	  velocity,       &
	  timedot,        &
	  tmp_time    
	  
real, allocatable, target, dimension(:,:,:) :: &
   	  displacement,	&
   	  tmp_disp,       &
   	  disspectrum,	&
   	  normvector,     &
   	  unit_norm
   	  
real, allocatable, target, dimension(:,:) :: &
	  xyzstation,	&
	  xx,			&
	  yy,			&
	  zz,			&
	  slrx,           &
	  slry,           &
	  slrz,           &
	  slrrx,          &
	  slrry,          &
	  slrrz,          &
	  area,           &
	  trup
	  
real, allocatable, target, dimension(:) :: &
        subarea
	  
double complex, allocatable, target, dimension(:) :: &
	  in,  		&
	  out	
real :: radius, rref !(rref set in simulation; radius is measured from trup)
! spectra portion
real :: sig0, fcorn, falloffbest, fitrms
real,allocatable,dimension(:,:) :: &
	asig0, afcorn, afalloffbest, afitrms
real,dimension(2) :: safc,safr,samo, energy!per processor, for reduction
! spherical average corner freq, fall-off rate, spectral level(moment),energy
real,dimension(2) :: safc0,safr0,samo0,energy0! for receiving reduce result in master node
real :: saar,saar0

! spherical area(discrete)
real,allocatable,dimension(:,:,:) :: &
      fitspectrum

real, allocatable, target, dimension(:) :: &
	  spec,		&
	  specfit
	  	  		  	  
end module

