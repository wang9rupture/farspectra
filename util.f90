module m_util
implicit none
contains

integer function deltaf(p, q)
integer,intent(in) :: p,q
if (p == q) then
deltaf=1
else
deltaf=0
end if
end function

real function cjkpq(j,k,p,q)
use m_globals
integer,intent(in) :: j,k,p,q
cjkpq=lam*deltaf(j,k)*deltaf(p,q)+mu*(deltaf(j,p) &
  *deltaf(k,q)+deltaf(j,q)*deltaf(k,p))
end function

real function coefp(i,j,m0,gamma,nv)
use m_globals
real,intent(in) :: m0
real,intent(in) :: gamma(:),nv(:)
integer,intent(in) :: i,j
integer :: k,p,q
coefp=0;
do k=1,3
do p=1,3
do q=1,3
coefp=coefp+cjkpq(j,k,p,q)*gamma(i)*gamma(p)* &
 gamma(q)*nv(k)/(4.0*pi*rho*vp**3*m0)
end do
end do
end do
end function
      
      
real function coefs(i,j,m0,gamma,nv)
use m_globals
real, intent(in) :: m0
real, intent(in) :: gamma(:),nv(:)
integer,intent(in) :: i,j
integer :: k,p,q
coefs=0;
do k=1,3
do p=1,3
do q=1,3
coefs=coefs-cjkpq(j,k,p,q)*(gamma(i)*gamma(p)- &
 deltaf(i,p))*gamma(q)*nv(k) &
 /(4.0*pi*rho*vs**3*m0)
end do
end do
end do
end function  

subroutine rupradius    
use m_globals
integer :: i,j
real,allocatable,dimension(:,:) :: ruparea

allocate(ruparea(nx,ny))
ruparea = 0.0

open(101,file='in/trup',recl=nx*ny*4,access='direct')
read(101,rec=1) ((trup(i,j),i=1,nx),j=1,ny)
close(101)

where(trup < 1e8 .and. trup > 0) 
ruparea = area
end where
radius = sqrt(sum(ruparea)/pi)
write(0,*) 'Rupture Radius = ', radius ,' m'
end subroutine 

subroutine differ(intmsers,outmsers)
real,intent(in) :: intmsers(:)
real,dimension(size(intmsers)),intent(out):: outmsers
integer :: n,i
n=size(intmsers)
Loop: do i=2,n-1
outmsers(i)=(intmsers(i+1)-intmsers(i-1))/2
end do Loop
outmsers(1)=outmsers(2)
outmsers(n)=outmsers(n-1)
end subroutine
end module
