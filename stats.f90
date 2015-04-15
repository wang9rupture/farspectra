module m_stats
implicit none
real,parameter :: up=0.42,us=0.59
contains

subroutine stats
use m_globals
use mpi
integer :: ierr
real :: k(2), strdrop(2)
! reduce spherical area
call mpi_allreduce(saar,saar0,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
call mpi_allreduce(safc,safc0,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
call mpi_allreduce(safr,safr0,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
call mpi_allreduce(samo,samo0,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
call mpi_allreduce(energy,energy0,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
! sum divided by sphere area
samo0 = samo0/saar0
safc0 = safc0/saar0
safr0 = safr0/saar0

samo0(1)=samo0(1)*4*pi*rho*dist*vp**3/up
samo0(2)=samo0(2)*4*pi*rho*dist*vs**3/us

strdrop = 7./16.*samo0/radius**3


if(myid == master) then
write(0,*) 'Spherical area = ', saar0, 'm^2'
write(0,*) 'spherical average fc = ', safc0, 'Hz P & S'
write(0,*) 'Spherical average normalized fc*a/beta = ',safc0*rref/vs, 'set radius'
write(0,*) 'Spherical average normalized fc*a/beta = ',safc0*radius/vs, 'measured radius'
write(0,*) 'Ratio of P wave fc to S wave fc =', safc0(1)/safc0(2)
write(0,*) 'spherical average fall off rate = ', safr0, ' P & S'
write(0,*) 'spherical average moment = ',samo0, ' J P & S'
write(0,*) 'Radiated energy  = ', sum(energy0),'J     P & S is',energy0 ,'J'
write(0,*) 'S/P Energy ratio = ', energy0(2)/energy0(1)
write(0,*) 'Er/M0 ratio of radiated energy to moment =', sum(energy0)/samo0(1), &
            sum(energy0)/samo0(2),' from P and S'
write(0,*) 'Stress drop from moment =', strdrop, 'From P & S '
write(0,*) 'Radiation efficiencies =', 2*mu*sum(energy0)/(strdrop*samo0), 'inferred from P &S'
write(0,*)
end if

end subroutine


end module