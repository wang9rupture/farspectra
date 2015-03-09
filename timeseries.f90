module m_timeseries
implicit none
contains


subroutine compute_timeseries(it)
use m_globals
use m_util
use mpi 
integer,intent(in) :: it
integer :: nsta,nc,n,inx,iny
integer :: i,j,k,l
real  ::  delay,tmp
real  ::  r(3),m0,gamma(3),xi(3),nv(3)
real  ::   sliprate(3)
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
do nc=1,3
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
!if (mod(it-1,500)==0) print *,xi,nv  
sliprate(1)=slrx(inx,iny)    
sliprate(2)=slry(inx,iny)
sliprate(3)=slrz(inx,iny)
   
!     for P wave seismogram      
delay=m0/vp      
timeseries(1,nc,it+floor(delay/dt),nsta)= &
    timeseries(1,nc,it+floor(delay/dt),nsta)+ &
    coefp(nc,1,m0,gamma,nv)*sliprate(1)+      &
    coefp(nc,2,m0,gamma,nv)*sliprate(2)+      &
    coefp(nc,3,m0,gamma,nv)*sliprate(3)
!if (mod(it-1,500)==0) print *,myid,m0,gamma,coefp(nc,1,m0,gamma,nv)
!     for S wave seismogram      
delay=m0/vs      
timeseries(2,nc,it+floor(delay/dt),nsta)= &
    timeseries(2,nc,it+floor(delay/dt),nsta)+  &
    coefs(nc,1,m0,gamma,nv)*sliprate(1)+       &
    coefs(nc,2,m0,gamma,nv)*sliprate(2)+       &
    coefs(nc,3,m0,gamma,nv)*sliprate(3)                 

end do
end do
end do
end do

end subroutine

subroutine compute_displacement
use m_globals
use mpi
integer :: info,fh,offset
!      peak=0.
displacement = sqrt(sum(timeseries*timeseries,2))
print *,maxval(timeseries(1,1,:,1))
!call MPI_FILE_OPEN(MPI_COMM_WORLD,'timeseries',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
!offset = myid*stnum*2*3*ntt
!call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
!call mpi_file_write(fh,timeseries(1,1,1,1),stnum*2*3*ntt,mpi_real,mpi_status_ignore,ierr)
!call MPI_FILE_CLOSE(fh,ierr)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!      normalized
!      if(displacement(nsta,nps,mntt).gt.peak) then
!      peak=displacement(nsta,nps,mntt)
!      endif	
!	do mntt=1,ntt
!	displacement(nsta,nps,mntt)=displacement(nsta,nps,
!     $ mntt)/peak
!      enddo
	
write(0,*) 'Calculation of Displacement Ends'
end subroutine

end module