module m_arrays
implicit none
contains

subroutine arrays
use m_globals
use mpi
integer :: l,k,j
real :: theta, phi
!print *,myid,'of',numprocs
! allocate arrays
if(mod(nst,numprocs)==0) then
intervsta = nst/numprocs
else
intervsta = nst/numprocs+1
end if
do l=1,intervsta
stnum=l
if(l+myid*intervsta > nst ) then
stnum=stnum-1
exit
end if
end do
!!
! allocate all dimension-dependent arrays
!!
allocate ( &
	timeseries(2,3,ntt,stnum),   &
	displacement(2,ntt,stnum),   &
	disspectrum(2,nfreq,stnum), &
	xyzstation(stnum,3),	   &	
	asig0(2,stnum),      &
	afcorn(2,stnum),  	&
	afalloffbest(2,stnum),	&
	afitrms(2,stnum),	&
	fitspectrum(2,nfreq,stnum),&
	slrx(nx,ny),            &
	slry(nx,ny),		&
	slrz(nx,ny),		&
	vrup(nx,3),             &
	in(ntt),			&
	out(ntt),			&
	spec(nfreq),		&
	specfit(nfreq)		&
)
timeseries = 0.0
displacement = 0.0
xyzstation = 0.0
slrx = 0.0
slry = 0.0
slrz = 0.0
vrup = 0.0
disspectrum = 0.0

      
!! read station distributions
!     reading station information
!     nst is number of stations
!     spherical domain
do l=1,stnum

j = mod(l+intervsta*myid,floor(360/degint))
k = (l+intervsta*myid-j)/floor(360/degint) + 1
phi=deg0+(j-1)*degint
theta=-90+(k-1)*degint

xyzstation(l,1)=origin(1)+dist*cos(theta*d2r)*cos(phi*d2r)
xyzstation(l,2)=origin(2)-dist*sin(theta*d2r)
xyzstation(l,3)=origin(3)+dist*cos(theta*d2r)*sin(phi*d2r)

end do

end subroutine

end module