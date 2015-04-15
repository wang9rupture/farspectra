module m_timeseries
implicit none
contains


subroutine compute_timeseries_1(it)
use m_globals
use m_util
use mpi 
integer,intent(in) :: it
integer :: nsta,nc,n,inx,iny
integer :: i,j,k,l,ierr
real  ::  delay,tmp1(3),tmp2(3),tmp3(3),subfault
real  ::  r(3),m0,gamma(3),xi(3),nv(3),unv(3)
real  ::  the,phi
real  :: locord(3,3),x(3)
real  :: fact,a_fp(3),a_fs(3)
real  :: up(3),us(3),u2p(3),u2s(3)
real  ::   sliprate(3),moslipr
real  :: tnow,tp,ts
if (myid == master) then
read(29,rec=it) ((slrx(j,k),j=1,nx),k=1,ny)
read(30,rec=it) ((slry(j,k),j=1,nx),k=1,ny)
read(31,rec=it) ((slrz(j,k),j=1,nx),k=1,ny)
end if
!! Broadcast slrx,slry,slrz
call MPI_BCAST(slrx(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slry(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slrz(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!if (mod(it-1,500)==0) print *,myid,origin,normvector(:,25,25)

do nsta=1,stnum
!if (mod(it-1,500)==0) print *,myid,xyzstation(nsta,:)
!     compute each station	
!write(0,*) 'staion is ',xyzstation(nsta,1),xyzstation(nsta,2) &
!      ,xyzstation(nsta,3)  
!     compute each component x-y-z
do iny=1,ny
!     compute y-direction subfault      
do inx=1,nx
!     compute x-direction subfault

xi(1)=xx(inx,iny)
xi(2)=yy(inx,iny)
xi(3)=zz(inx,iny)

r(:)=xyzstation(nsta,:)-xi(:)
m0=sqrt(sum(r*r))
gamma = r/m0

nv(1:3)=normvector(1:3,inx,iny)
unv(1:3)=unit_norm(1:3,inx,iny)

!if (mod(it-1,500)==0) print *,xi,nv  
sliprate(1)=slrx(inx,iny)    
sliprate(2)=slry(inx,iny)
sliprate(3)=slrz(inx,iny)
moslipr = sqrt(dot_product(sliprate,sliprate))
subfault = area(inx,iny) 

do nc=1,3  
!     for P wave seismogram      
delay=m0/vp      
timeseries(1,nc,it+floor(delay/dt),nsta)= &
    timeseries(1,nc,it+floor(delay/dt),nsta)+ &
    coefp(nc,1,m0,gamma,unv)*sliprate(1)*subfault+      &
    coefp(nc,2,m0,gamma,unv)*sliprate(2)*subfault+      &
    coefp(nc,3,m0,gamma,unv)*sliprate(3)*subfault
    
!if (mod(it-1,500)==0) print *,myid,m0,gamma,coefp(nc,1,m0,gamma,nv)
!     for S wave seismogram      
delay=m0/vs      
timeseries(2,nc,it+floor(delay/dt),nsta)= &
    timeseries(2,nc,it+floor(delay/dt),nsta)+  &
    coefs(nc,1,m0,gamma,unv)*sliprate(1)*subfault+       &
    coefs(nc,2,m0,gamma,unv)*sliprate(2)*subfault+       &
    coefs(nc,3,m0,gamma,unv)*sliprate(3)*subfault                 

end do
end do
end do
end do

end subroutine

subroutine compute_timeseries_2(it)
use m_globals
use m_util
use mpi 
integer,intent(in) :: it
integer :: nsta,nc,n,inx,iny
integer :: i,j,k,l,ierr
real  ::  delay,tmp1(3),tmp2(3),tmp3(3),subfault
real  ::  r(3),m0,gamma(3),xi(3),nv(3),unv(3)
real  ::  the,phi
real  :: locord(3,3),x(3)
real  :: fact,a_fp(3),a_fs(3)
real  :: up(3),us(3),u2p(3),u2s(3)
real  ::   sliprate(3),moslipr
real  :: tnow,tp,ts
if (myid == master) then
read(29,rec=it) ((slrx(j,k),j=1,nx),k=1,ny)
read(30,rec=it) ((slry(j,k),j=1,nx),k=1,ny)
read(31,rec=it) ((slrz(j,k),j=1,nx),k=1,ny)
end if
!! Broadcast slrx,slry,slrz
call MPI_BCAST(slrx(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slry(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slrz(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!if (mod(it-1,500)==0) print *,myid,origin,normvector(:,25,25)

do nsta=1,stnum
!if (mod(it-1,500)==0) print *,myid,xyzstation(nsta,:)
!     compute each station	
!write(0,*) 'staion is ',xyzstation(nsta,1),xyzstation(nsta,2) &
!      ,xyzstation(nsta,3)  
!     compute each component x-y-z
do iny=1,ny
!     compute y-direction subfault      
do inx=1,nx
!     compute x-direction subfault

up = 0.0
us = 0.0
u2p = 0.0
u2s = 0.0

xi(1)=xx(inx,iny)
xi(2)=yy(inx,iny)
xi(3)=zz(inx,iny)

r(:)=xyzstation(nsta,:)-xi(:)
m0=sqrt(sum(r*r))
gamma = r/m0

unv(1:3)=unit_norm(1:3,inx,iny)

!if (mod(it-1,500)==0) print *,xi,nv  
sliprate(1)=slrx(inx,iny)    
sliprate(2)=slry(inx,iny)
sliprate(3)=slrz(inx,iny)
moslipr = sqrt(dot_product(sliprate,sliprate))
subfault = area(inx,iny) 
locord =0.0
if(moslipr<1e-20) then
locord(1,1)=1
locord(2,2)=1
locord(3,3)=1
else
!convert to spherical coordinate (use unit_norm)
! local coordinate
locord(3,:) = unv
locord(1,:) = sliprate/moslipr
tmp1(1:3)=locord(3,1:3)
tmp2(1:3)=locord(1,1:3)
tmp3 = cross_product(tmp1,tmp2)
locord(2,1:3)=tmp3(1:3)
end if
! convert local to spherical coordinate 
x(1) = dot_product(gamma,locord(1,:))
x(2) = dot_product(gamma,locord(2,:))
x(3) = dot_product(gamma,locord(3,:))

the = acos(x(3))
phi = atan2(x(2),x(1))

   a_fp(1) =  sin(2.*the)*cos(phi)
   a_fp(2) =  0.
   a_fp(3) =  0.
!   print *, 'a_fp = ', a_fp   
    
   a_fs(1) =  0.
   a_fs(2) =  cos(2.*the)*cos(phi)
   a_fs(3) = -cos(the)*sin(phi)
!   print *, 'a_fs = ', a_fs  

fact = 1./(4.*pi*rho*vp**3*m0)
up(1) = fact*a_fp(1)*mu*subfault*moslipr
up(2) = fact*a_fp(2)*mu*subfault*moslipr
up(3) = fact*a_fp(3)*mu*subfault*moslipr

fact = 1./(4.*pi*rho*vs**3*m0)
us(1) = fact*a_fs(1)*mu*subfault*moslipr
us(2) = fact*a_fs(2)*mu*subfault*moslipr
us(3) = fact*a_fs(3)*mu*subfault*moslipr

!convert spherical to local cartesian 
u2p(1) = up(1)*sin(the)*cos(phi) + up(2)*cos(the)*cos(phi) - up(3)*sin(phi)
u2p(2) = up(1)*sin(the)*sin(phi) + up(2)*cos(the)*sin(phi) + up(3)*cos(phi)
u2p(3) = up(1)*cos(the) - up(2)*sin(the)

u2s(1) = us(1)*sin(the)*cos(phi) + us(2)*cos(the)*cos(phi) - us(3)*sin(phi)
u2s(2) = us(1)*sin(the)*sin(phi) + us(2)*cos(the)*sin(phi) + us(3)*cos(phi)
u2s(3) = us(1)*cos(the) - us(2)*sin(the)

!convert local cartesian to global cartesian
tp = m0/vp
ts = m0/vs
tnow=real(it-1)*dt
! p wave
k = nint((tnow+tp)/dt)
timeseries(1,1,k,nsta) = timeseries(1,1,k,nsta) + &
    u2p(1)*locord(1,1)+u2p(2)*locord(2,1)+u2p(3)*locord(3,1)
timeseries(1,2,k,nsta) = timeseries(1,2,k,nsta) + &
    u2p(1)*locord(1,2)+u2p(2)*locord(2,2)+u2p(3)*locord(3,2)
timeseries(1,3,k,nsta) = timeseries(1,3,k,nsta) + &
    u2p(1)*locord(1,3)+u2p(2)*locord(2,3)+u2p(3)*locord(3,3)
! s wave
k = nint((tnow+ts)/dt)
timeseries(2,1,k,nsta) = timeseries(2,1,k,nsta) + &
    u2s(1)*locord(1,1)+u2s(2)*locord(2,1)+u2s(3)*locord(3,1)
timeseries(2,2,k,nsta) = timeseries(2,2,k,nsta) + &
    u2s(1)*locord(1,2)+u2s(2)*locord(2,2)+u2s(3)*locord(3,2)
timeseries(2,3,k,nsta) = timeseries(2,3,k,nsta) + &
    u2s(1)*locord(1,3)+u2s(2)*locord(2,3)+u2s(3)*locord(3,3)

end do
end do
end do

end subroutine

subroutine compute_timeseries_3(it)
use m_globals
use m_util
use mpi 
integer,intent(in) :: it
integer :: nsta,nc,n,inx,iny
integer :: i,j,k,l,ierr
real  ::  delay,tmp1(3),tmp2(3),tmp3(3),subfault
real  ::  r(3),m0,gamma(3),xi(3),nv(3),unv(3)
real  ::  the,phi
real  :: locord(3,3),x(3)
real  :: fact
real  ::   sliprate(3),moslipr
real  :: tnow,tp,ts
if (myid == master) then
read(29,rec=it) ((slrx(j,k),j=1,nx),k=1,ny)
read(30,rec=it) ((slry(j,k),j=1,nx),k=1,ny)
read(31,rec=it) ((slrz(j,k),j=1,nx),k=1,ny)
end if
!! Broadcast slrx,slry,slrz
call MPI_BCAST(slrx(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slry(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(slrz(1,1),nx*ny,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!if (mod(it-1,500)==0) print *,myid,origin,normvector(:,25,25)
do nsta=1,stnum
!if (mod(it-1,500)==0) print *,myid,xyzstation(nsta,:)
!     compute each station	
!write(0,*) 'staion is ',xyzstation(nsta,1),xyzstation(nsta,2) &
!      ,xyzstation(nsta,3)  
!     compute each component x-y-z
do iny=1,ny
!     compute y-direction subfault      
do inx=1,nx
!     compute x-direction subfault

xi(1)=xx(inx,iny)
xi(2)=yy(inx,iny)
xi(3)=zz(inx,iny)

r(:)=xyzstation(nsta,:)-xi(:)
m0=sqrt(sum(r*r))
gamma = r/m0
!gamma(1) = sin(60*d2r)*cos(60*d2r)
!gamma(2) = sin(60*d2r)*sin(60*d2r)
!gamma(3) = cos(60*d2r)
!m0 = dist
nv(1:3)=normvector(1:3,inx,iny)
unv(1:3)=unit_norm(1:3,inx,iny)

!if (mod(it-1,500)==0) print *,xi,nv  
sliprate(1)=slrx(inx,iny)    
sliprate(2)=slry(inx,iny)
sliprate(3)=slrz(inx,iny)
moslipr = sqrt(dot_product(sliprate,sliprate))
subfault = area(inx,iny) 

  
!     for P wave seismogram      
tp = m0/vp
ts = m0/vs
tnow=real(it-1)*dt
! p wave
k = nint((tnow+tp)/dt)
fact = 2.*subfault*mu*sum(gamma*unv)/(4*pi*rho*vp**3*m0)      
timeseries(1,1,k,nsta)= &
    timeseries(1,1,k,nsta)+ &
    fact*gamma(1)*gamma(1)*sliprate(1) + &
    fact*gamma(1)*gamma(2)*sliprate(2) + &
    fact*gamma(1)*gamma(3)*sliprate(3)

timeseries(1,2,k,nsta)= &
    timeseries(1,2,k,nsta)+ &
    fact*gamma(2)*gamma(1)*sliprate(1) + &
    fact*gamma(2)*gamma(2)*sliprate(2) + &
    fact*gamma(2)*gamma(3)*sliprate(3)
    
timeseries(1,3,k,nsta)= &
    timeseries(1,3,k,nsta)+ &
    fact*gamma(3)*gamma(1)*sliprate(1) + &
    fact*gamma(3)*gamma(2)*sliprate(2) + &
    fact*gamma(3)*gamma(3)*sliprate(3)    

velocity(1,1,k,nsta) = &
    velocity(1,1,k,nsta)+ &
    fact*gamma(1)*gamma(1)*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*gamma(1)*gamma(2)*(sliprate(2)-slrry(inx,iny))/dt + &
    fact*gamma(1)*gamma(3)*(sliprate(3)-slrrz(inx,iny))/dt
    
velocity(1,2,k,nsta) = &
    velocity(1,2,k,nsta)+ &
    fact*gamma(2)*gamma(1)*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*gamma(2)*gamma(2)*(sliprate(2)-slrry(inx,iny))/dt + &
    fact*gamma(2)*gamma(3)*(sliprate(3)-slrrz(inx,iny))/dt
    
velocity(1,3,k,nsta) = &
    velocity(1,3,k,nsta)+ &
    fact*gamma(3)*gamma(1)*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*gamma(3)*gamma(2)*(sliprate(2)-slrry(inx,iny))/dt + &
    fact*gamma(3)*gamma(3)*(sliprate(3)-slrrz(inx,iny))/dt

!if (mod(it-1,500)==0) print *,myid,m0,gamma,coefp(nc,1,m0,gamma,nv)
!     for S wave seismogram      
k = nint((tnow+ts)/dt) 
fact = mu*subfault/(4*pi*rho*vs**3*m0)      
timeseries(2,1,k,nsta)= &
    timeseries(2,1,k,nsta)+  &
    fact*(deltaf(1,1)*sum(gamma*unv)+unv(1)*gamma(1)-2*gamma(1) &
    *gamma(1)*sum(gamma*unv))*sliprate(1) + &
    fact*(deltaf(1,2)*sum(gamma*unv)+unv(1)*gamma(2)-2*gamma(1) &
    *gamma(2)*sum(gamma*unv))*sliprate(2) + & 
    fact*(deltaf(1,3)*sum(gamma*unv)+unv(1)*gamma(3)-2*gamma(1) &
    *gamma(3)*sum(gamma*unv))*sliprate(3)  
      
timeseries(2,2,k,nsta)= &
    timeseries(2,2,k,nsta)+  &
    fact*(deltaf(2,1)*sum(gamma*unv)+unv(2)*gamma(1)-2*gamma(2) &
    *gamma(1)*sum(gamma*unv))*sliprate(1) + &
    fact*(deltaf(2,2)*sum(gamma*unv)+unv(2)*gamma(2)-2*gamma(2) &
    *gamma(2)*sum(gamma*unv))*sliprate(2) + & 
    fact*(deltaf(2,3)*sum(gamma*unv)+unv(2)*gamma(3)-2*gamma(2) &
    *gamma(3)*sum(gamma*unv))*sliprate(3) 
    
timeseries(2,3,k,nsta)= &
    timeseries(2,3,k,nsta)+  &
    fact*(deltaf(3,1)*sum(gamma*unv)+unv(3)*gamma(1)-2*gamma(3) &
    *gamma(1)*sum(gamma*unv))*sliprate(1) + &
    fact*(deltaf(3,2)*sum(gamma*unv)+unv(3)*gamma(2)-2*gamma(3) &
    *gamma(2)*sum(gamma*unv))*sliprate(2) + & 
    fact*(deltaf(3,3)*sum(gamma*unv)+unv(3)*gamma(3)-2*gamma(3) &
    *gamma(3)*sum(gamma*unv))*sliprate(3)  
    
velocity(2,1,k,nsta)= &
    velocity(2,1,k,nsta)+  &
    fact*(deltaf(1,1)*sum(gamma*unv)+unv(1)*gamma(1)-2*gamma(1) &
    *gamma(1)*sum(gamma*unv))*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*(deltaf(1,2)*sum(gamma*unv)+unv(1)*gamma(2)-2*gamma(1) &
    *gamma(2)*sum(gamma*unv))*(sliprate(2)-slrry(inx,iny))/dt + & 
    fact*(deltaf(1,3)*sum(gamma*unv)+unv(1)*gamma(3)-2*gamma(1) &
    *gamma(3)*sum(gamma*unv))*(sliprate(3)-slrrz(inx,iny))/dt    

velocity(2,2,k,nsta)= &
    velocity(2,2,k,nsta)+  &
    fact*(deltaf(2,1)*sum(gamma*unv)+unv(2)*gamma(1)-2*gamma(2) &
    *gamma(1)*sum(gamma*unv))*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*(deltaf(2,2)*sum(gamma*unv)+unv(2)*gamma(2)-2*gamma(2) &
    *gamma(2)*sum(gamma*unv))*(sliprate(2)-slrry(inx,iny))/dt + & 
    fact*(deltaf(2,3)*sum(gamma*unv)+unv(2)*gamma(3)-2*gamma(2) &
    *gamma(3)*sum(gamma*unv))*(sliprate(3)-slrrz(inx,iny))/dt   
    
velocity(2,3,k,nsta)= &
    velocity(2,3,k,nsta)+  &
    fact*(deltaf(3,1)*sum(gamma*unv)+unv(3)*gamma(1)-2*gamma(3) &
    *gamma(1)*sum(gamma*unv))*(sliprate(1)-slrrx(inx,iny))/dt + &
    fact*(deltaf(3,2)*sum(gamma*unv)+unv(3)*gamma(2)-2*gamma(3) &
    *gamma(2)*sum(gamma*unv))*(sliprate(2)-slrry(inx,iny))/dt + & 
    fact*(deltaf(3,3)*sum(gamma*unv)+unv(3)*gamma(3)-2*gamma(3) &
    *gamma(3)*sum(gamma*unv))*(sliprate(3)-slrrz(inx,iny))/dt                    
end do
end do
end do
slrrx=slrx
slrry=slry
slrrz=slrz
end subroutine


subroutine compute_displacement
use m_globals
use m_util
use mpi
integer :: info,fh,ierr,i,j,k
integer(kind=mpi_offset_kind) :: offset
real :: tmp_timedot(2,ntt,stnum),tmp_energy(2)
real,dimension(ntt) :: array
!      peak=0.
energy = 0.0
! assign to 0 if small enough
!where(timeseries < 1e-7) 
!    timeseries = 0.0
!end where

!print *,myid,sum(timeseries)

displacement = sqrt(sum(timeseries*timeseries,2))
tmp_timedot = sqrt(sum(timedot*timedot,2))
tmp_timedot = sqrt(sum(velocity*velocity,2))
do j=1,stnum
       energy(1) = energy(1) + subarea(j)* &
                   sum(tmp_timedot(1,:,j)*tmp_timedot(1,:,j))*dt*vp*rho
       energy(2) = energy(2) + subarea(j)* &
                   sum(tmp_timedot(2,:,j)*tmp_timedot(2,:,j))*dt*vs*rho       
end do
!output time series records if necessary

!print *,maxval(timeseries(1,1,:,1))
!call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/timedot',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
!offset = myid*stnum*2*3*ntt
!call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
!call mpi_file_write(fh,timedot(1,1,1,1),stnum*2*3*ntt,mpi_real,mpi_status_ignore,ierr)
!call MPI_FILE_CLOSE(fh,ierr)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/velocity',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
!offset = myid*stnum*2*3*ntt
!call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
!call mpi_file_write(fh,velocity(1,1,1,1),stnum*2*3*ntt,mpi_real,mpi_status_ignore,ierr)
!call MPI_FILE_CLOSE(fh,ierr)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/displacement',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2*ntt
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,displacement(1,1,1),stnum*2*ntt,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!
!call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/tmp_disp',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
!offset = myid*stnum*2*ntt
!call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
!call mpi_file_write(fh,tmp_disp(1,1,1),stnum*2*ntt,mpi_real,mpi_status_ignore,ierr)
!call MPI_FILE_CLOSE(fh,ierr)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!
!      normalized
!      if(displacement(nsta,nps,mntt).gt.peak) then
!      peak=displacement(nsta,nps,mntt)
!      endif	
!	do mntt=1,ntt
!	displacement(nsta,nps,mntt)=displacement(nsta,nps,
!     $ mntt)/peak
!      enddo
	
!write(0,*) 'Calculation of Displacement Ends'
end subroutine

end module
