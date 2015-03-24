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
        getfc,    &
        cfc
        
        
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
	  soff 

integer :: &
	  myid, numprocs, intervsta, stnum
	  
	  
real, allocatable, target, dimension(:,:,:,:) :: &
	  timeseries
	  
real, allocatable, target, dimension(:,:,:) :: &
   	  displacement,	&
   	  disspectrum,	&
   	  normvector
   	  
real, allocatable, target, dimension(:,:) :: &
	  xyzstation,	&
	  xx,			&
	  yy,			&
	  zz,			&
	  slrx,           &
	  slry,           &
	  slrz,           &
	  area,           &
	  trup
	  
double complex, allocatable, target, dimension(:) :: &
	  in,  		&
	  out,            &
	  subarea	
	  
! spectra portion
real :: sig0, fcorn, falloffbest, fitrms
real,allocatable,dimension(:,:) :: &
	asig0, afcorn, afalloffbest, afitrms
real,allocatable,dimension(:,:,:) :: &
      fitspectrum

real, allocatable, target, dimension(:) :: &
	  spec,		&
	  specfit
	  	  		  	  
end module

