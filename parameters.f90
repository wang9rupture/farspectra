! Read model parameters
module m_parameters
implicit none
contains

subroutine read_parameters
use m_globals
use m_arrays
use mpi
integer :: j,l,k,ierr
open( 22, file='parameters', status='old',iostat = k)

read(22,*) yorp
read(22,*) lx
read(22,*) ly
read(22,*) lz
read(22,*) dx
read(22,*) dy
read(22,*) dz
read(22,*) t
read(22,*) tt
read(22,*) rho
read(22,*) vp
read(22,*) vs
read(22,*) soff
read(22,*) getfc
read(22,*) psfrqf
read(22,*) fc1
read(22,*) fc2
read(22,*) cfc
read(22,*) falloff1in
read(22,*) falloff2in
read(22,*) dist
read(22,*) degint
read(22,*) deg0
read(22,*) rref


close(22)
if(myid == master) then
write(0,*) 'Use Yongfei(1) or Peter(2)', yorp
write(0,*) 'Dimension x meter ', lx
write(0,*) 'Dimension y meter ',ly
write(0,*) 'Dimension z meter ',lz
write(0,*) 'dx meter',dx
write(0,*) 'dy meter',dy
write(0,*) 'dz meter',dz
write(0,*) 'src period sec ',t
write(0,*) 'rceiv period sec ',tt
write(0,*) 'density kg/m^3 ',rho
write(0,*) 'P wave velocity m/s ',vp
write(0,*) 'S wave velocity m/s ',vs
write(0,*) 'Only measure spect switch 1-no 2-yes ',soff
write(0,*) 'Use getcorn or getcorn2 to measure fc', getfc
write(0,*) 'Fitting P freq band Hz',psfrqf(1),psfrqf(2)
write(0,*) 'Fitting S freq band Hz',psfrqf(3),psfrqf(4)
write(0,*) 'Min fc',fc1
write(0,*) 'Max fc',fc2
write(0,*) 'weight frequency Hz',cfc
write(0,*) 'Minimum falloff ',abs(falloff1in)
write(0,*) 'Maximum falloff ',falloff2in
write(0,*) 'Distance meter ',dist
write(0,*) 'Degree interval ',degint
write(0,*) 'Original azimuth degree',deg0
write(0,*) 'Set rupture radius m', rref
write(0,*) ''
end if
   
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
mu=vs**2*rho
lam=vp**2*rho-2*mu

! compute dimensions
! make preparation for array.f90
nx=floor(lx/dx+1.5)
ny=floor(ly/dy+1.5)
nz=floor(lz/dz+2.5)
dt=dx/12500.0
nt=floor(t/dt+1.5)
ntt=floor(tt/dt+1.5)
df=1./(dt*(ntt-1))
nst=floor(360.0/degint)*(floor(180.0/degint))
nfreq=(ntt-1)/2
allocate( &
      xx(nx,ny),			&
	yy(nx,ny),			&
	zz(nx,ny),			&
	normvector(3,nx,ny),	&
	unit_norm(3,nx,ny)      &
	)
xx = 0.0
yy = 0.0
zz = 0.0
normvector = 0.0
if (myid==master) then
write(0,*) 'Dimension check:'
write(0,*) 'nx ny nz ',nx,ny,nz
write(0,*) 'dt nt ntt ',dt,nt,ntt
write(0,*) 'df nst nfreq ',df,nst,nfreq
end if
!     set ihypo at domain center
ihypo(1)=(nx+1)/2
ihypo(2)=(ny+1)/2
ihypo(3)=(nz)/2

!! root node reading files
if (myid == master) then
!     reading in slip coordinate function
  open(26,file='in/xx',access='direct', &
    recl=nx*ny*sizereal,form='unformatted',status='old')
  open(27,file='in/yy',access='direct', &
    recl=nx*ny*sizereal,form='unformatted',status='old')
  open(28,file='in/zz',access='direct', &
    recl=nx*ny*sizereal,form='unformatted',status='old')
	  
  read(26,rec=ihypo(3)) ((xx(j,l),j=1,nx),l=1,ny)
  read(27,rec=ihypo(3)) ((yy(j,l),j=1,nx),l=1,ny)
  read(28,rec=ihypo(3)) ((zz(j,l),j=1,nx),l=1,ny)
      
  close(26)
  close(27)
  close(28)

!     compute origin coordinate
  origin(1)=xx(ihypo(1),ihypo(2))
  origin(2)=yy(ihypo(1),ihypo(2))
  origin(3)=zz(ihypo(1),ihypo(2))

!  read normal vector( magnitude is area)   
  open(66,file='in/n_x',access='direct', &
	recl=nx*ny*sizereal,form='unformatted',status='old')
  read(66,rec=1) ((normvector(1,j,l),j=1,nx),l=1,ny)
  close(66)
      
  open(67,file='in/n_y',access='direct', &
      recl=nx*ny*sizereal,form='unformatted',status='old')
  read(67,rec=1) ((normvector(2,j,l),j=1,nx),l=1,ny)
  close(67)
      
  open(68,file='in/n_z',access='direct', &
      recl=nx*ny*sizereal,form='unformatted',status='old')
  read(68,rec=1) ((normvector(3,j,l),j=1,nx),l=1,ny)
  close(68)
end if

call MPI_BCAST(origin(1),size(origin),MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(normvector(1,1,1),size(normvector),MPI_REAL,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(xx(1,1),size(xx),mpi_real,0,mpi_comm_world,ierr)
call mpi_bcast(yy(1,1),size(yy),mpi_real,0,mpi_comm_world,ierr)
call mpi_bcast(zz(1,1),size(zz),mpi_real,0,mpi_comm_world,ierr)
!print *,myid,xx(10,10),yy(10,10),zz(10,10),normvector(:,1,1)
call arrays
area=sqrt(sum(normvector*normvector, 1))
unit_norm(1,:,:) = normvector(1,:,:)/area(:,:)
unit_norm(2,:,:) = normvector(2,:,:)/area(:,:)
unit_norm(3,:,:) = normvector(3,:,:)/area(:,:)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine

end module
