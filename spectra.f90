module m_spectra
implicit none
contains

subroutine calculate_disspectrum
use m_globals
use m_FFTW3
use mpi
integer :: i,j,k,l,ind
integer(kind=8) :: plan

!     length of spectrum	  
write(0,*) 'start to calculate spectrum'
!	  ind=1
do k=1,stnum
	do i=1,2
in(:)=dcmplx(displacement(i,:,k),0)
call dfftw_plan_dft_1d_ (plan,ntt,in,out,FFTW_FORWARD, &
       FFTW_ESTIMATE ) 
call dfftw_execute_ ( plan ) 
call dfftw_destroy_plan_ ( plan )     	  			
disspectrum(i,1:nfreq,k)=real(log10(abs(dt*out(1:nfreq))))      
   end do
end do

!write(0,*) 'Calculation of spectra Ends'
end subroutine

subroutine read_disspectrum
use m_globals
use mpi
integer :: k,l,i

end subroutine

!spec(used globally) entered
subroutine getcorn2
!! provided by Peter Shearer, but not well constrained
use m_globals
real,allocatable,dimension(:) :: resid
real, parameter :: dfall=0.1   !falloff interval
real, parameter :: fc1=0.01,fc2=2,dfc=0.01,cfc=0.5!corner freq range
real :: falloff1,falloff2,falloff
integer :: nfall,ifall,i
real :: f1,f2,fc,f
integer :: nf,nfc,kfc
real :: fitbest,pred,sum,wsum,weight,sig,fit
allocate(resid(ifrq1:ifrq2))
write(0,*) 'Enter getcorn'
if (ifrq2.gt.nfreq) then
	print *,'***ERROR in GETCORN: ',ifrq2,nfreq
      stop
end if

      falloff1 = abs(falloff1in)
      falloff2 = abs(falloff2in)
      
      f1=df
      f2=real(nfreq-1)*df
      nf=ifrq2-ifrq1+1
      
      fitbest=9.e30      
      nfall = 1+nint((falloff2-falloff1)/dfall)
      nfc = 1 + nint((fc2-fc1)/dfc)
      
      do ifall = 1, nfall
         falloff = falloff1 + real(ifall-1)*dfall
      
      do kfc = 1, nfc
        fc = fc1 + real(kfc-1)*dfc
         sum=0.
         wsum=0.
         do i=ifrq1,ifrq2
            f=real(i-1)*df
            pred = -alog10(1.+(f/fc)**falloff)
            
            if (falloff1in > 0.) then
               weight = 1.
            else
               weight = 1./max(f,cfc)
            end if
   		
            resid(i)=spec(i)-pred
            sum=sum+resid(i)*weight
            wsum = wsum + weight
         end do
         sig=sum/wsum
         
         fit=0.
         wsum=0.
         do i=ifrq1,ifrq2
            f=real(i-1)*df         
            if (falloff1in > 0.) then
               weight = 1.
            else
               weight = 1./max(f,cfc)
            end if         
         
            resid(i) = resid(i) - sig
            fit = fit + weight*resid(i)**2
            wsum = wsum + weight
         end do
         fit=sqrt(fit/wsum)
         
         if (fit.lt.fitbest) then
            fitbest=fit
            fcorn=fc
            falloffbest = falloff
            sig0=sig
         end if
      end do
end do

      fitrms=fitbest
do i=1,nfreq
      f=real(i-1)*df
      pred = sig0 - alog10(1.+(f/fcorn)**falloffbest)
      specfit(i)=pred
end do
write(0,*) 'Exit getcorn'
end subroutine
!specfit exited

subroutine getcorn
!! Modified based upon getcorn2 given by Peter Shearer
use m_globals
real,allocatable,dimension(:) :: resid
real, parameter :: dfall=0.1   !falloff interval
real, parameter :: fc1=0.01,fc2=2,dfc=0.01,cfc=1. !corner freq range
real :: falloff1,falloff2,falloff
integer :: nfall,ifall,i
real :: f1,f2,fc,f
integer :: nf,nfc,kfc
real :: fitbest,pred,sum,wsum,weight,sig,fit
allocate(resid(ifrq1:ifrq2))
write(0,*) 'Enter getcorn'
if (ifrq2.gt.nfreq) then
	print *,'***ERROR in GETCORN: ',ifrq2,nfreq
      stop
end if

      falloff1 = abs(falloff1in)
      falloff2 = abs(falloff2in)
      
      f1=df
      f2=real(nfreq-1)*df
      nf=ifrq2-ifrq1+1
      
      fitbest=9.e30      
      nfall = 1+nint((falloff2-falloff1)/dfall)
      nfc = 1 + nint((fc2-fc1)/dfc)
      
      sig0=spec(1)
      
      do ifall = 1, nfall
         falloff = falloff1 + real(ifall-1)*dfall
      
      do kfc = 1, nfc
        fc = fc1 + real(kfc-1)*dfc
         sum=0.
         wsum=0.
         do i=ifrq1,ifrq2
            f=real(i-1)*df
            pred = sig0-alog10(1.+(f/fc)**falloff)
            
            if (falloff1in > 0.) then
               weight = 1.
            else
               weight = 1./max(f,cfc)
            end if
   		
            resid(i)=spec(i)-pred
            sum=sum+weight*resid(i)**2
            wsum = wsum + weight
         end do

         fit=sqrt(sum/wsum)
         
         if (fit<=fitbest) then
            fitbest=fit
            fcorn=fc
            falloffbest = falloff
         end if
      end do
end do
      fitrms=fitbest
do i=1,nfreq
      f=real(i-1)*df
      pred = sig0 - alog10(1.+(f/fcorn)**falloffbest)
      specfit(i)=pred
end do
write(0,*) 'Exit getcorn'
end subroutine
!specfit exited


	  
subroutine calculate_discornfall
use m_globals
use mpi
integer :: i,j,k,l,ind
integer :: info,fh,offset
do j=1,stnum
do i=1,2
spec = 0
spec(:)=disspectrum(i,:,j)
!print *,''
!      cornerfall format:
!      sig0,fcorn,falloffbest,fitrms, specfit(1->nfreq)	
call  getcorn
print *,myid,sig0,fcorn,falloffbest,fitrms,specfit(1)
!(optionally call getcorn2)
asig0(i,j)=sig0
afcorn(i,j)=fcorn
afalloffbest(i,j)=falloffbest
afitrms(i,j)=fitrms
fitspectrum(i,:,j)=specfit(:)		  
end do
end do 

!call MPI_FILE_OPEN(MPI_COMM_WORLD,'fcorn',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
!offset = myid*stnum*2
!call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
!call mpi_file_write(fh,afcorn,stnum*2,mpi_real,mpi_status_ignore,ierr)
!call MPI_FILE_CLOSE(fh,ierr)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)   
end subroutine
      
        

end module
