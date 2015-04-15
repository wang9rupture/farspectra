program main
use m_globals
use m_parameters
use m_timeseries
use m_spectra
use m_util
use m_stats
use mpi
implicit none

integer :: it,ierr,rc

!! mpi initiation
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
!! 
!! master node read parameter and write files 
write(0,*) myid,'of',numprocs
call read_parameters
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
if(soff == 1) then
  if(myid == master) then

  call rupradius

  open(29,file='in/slipr_x',access='direct', &
     recl=nx*ny*sizereal,form='unformatted',status='old')
  open(30,file='in/slipr_y',access='direct', &
     recl=nx*ny*sizereal,form='unformatted',status='old')
  open(31,file='in/slipr_z',access='direct', &
     recl=nx*ny*sizereal,form='unformatted',status='old')
  end if
  do it=1,nt
      if(yorp == 1) then
        call compute_timeseries_3(it)
      elseif(yorp ==2) then
        call compute_timeseries_2(it)
      else
        call compute_timeseries_1(it)
      end if
      if(myid == master .and. mod(it-1,20).eq.0) then
      write(0,"(A,F5.1)") '% = ',real(it-1)/nt*100.0 
      end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end do 
  if(myid == master) then    
    close(29)
    close(30)
    close(31) 
    write(0,*) ''
    write(0,*) 'Calculation of Waveforms Ends'
    write(0,*) ''
  end if
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call compute_displacement
  call calculate_disspectrum
else
!call read_timeseries_displacement
  call read_disspectrum
end if
  call calculate_discornfall
  call stats
  call MPI_FINALIZE(rc)

end program
