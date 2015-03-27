module m_spectra
implicit none
contains

subroutine calculate_disspectrum
use m_globals
use m_FFTW3
integer :: i,j,k,l,ind
integer(kind=8) :: plan

!     length of spectrum	  
!write(0,*) 'start to calculate spectrum'
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
real, parameter :: fc1=0.01,fc2=5,dfc=0.01!corner freq range
real :: falloff1,falloff2,falloff
integer :: nfall,ifall,i
real :: f1,f2,fc,f
integer :: nf,nfc,kfc
real :: fitbest,pred,sum,wsum,weight,sig,fit
allocate(resid(ifrq1:ifrq2))
!write(0,*) 'Enter getcorn'
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
!inversely dependent on frequency  
               if (cfc < 0) then          
                 if((f-0.)<1e-8) then
                   weight = 0.
                 else	          
                   weight = 1/f
                 end if
! modified flat weight at low frequency
               else               
                 if(f<cfc) then
                   weight = 1
                 else 
                   weight = cfc/f
                 end if
               end if               
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
!inversely dependent on frequency  
               if (cfc < 0) then          
                 if((f-0.)<1e-8) then
                   weight = 0.
                 else	          
                   weight = 1/f
                 end if
! modified flat weight at low frequency
               else               
                 if(f<cfc) then
                   weight = 1
                 else 
                   weight = cfc/f
                 end if
               end if               
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
!write(0,*) 'Exit getcorn'
end subroutine
!specfit exited

subroutine getcorn
!! Modified based upon getcorn2 given by Peter Shearer
use m_globals
real,allocatable,dimension(:) :: resid
real, parameter :: dfall=0.1   !falloff interval
real, parameter :: fc1=0.01,fc2=5,dfc=0.01 !corner freq range
real :: falloff1,falloff2,falloff
integer :: nfall,ifall,i
real :: f1,f2,fc,f
integer :: nf,nfc,kfc
real :: fitbest,pred,sum,wsum,weight,sig,fit
allocate(resid(ifrq1:ifrq2))
!write(0,*) 'Enter getcorn'
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
            pred = sig0-log10(1.+(f/fc)**falloff)
            
            if (falloff1in > 0.) then
               weight = 1.
            else
!inversely dependent on frequency  
               if (cfc < 0) then
                 if((f-0.)<1e-8) then
                   weight = 0.
                 else	          
                   weight = 1/f
                 end if
! modified flat weight at low frequency
               else               
                 if(f<cfc) then
                   weight = 1
                 else 
                   weight = cfc/f
                 end if
               end if
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
      pred = sig0 - log10(1.+(f/fcorn)**falloffbest)
      specfit(i)=pred
end do
!write(0,*) 'Exit getcorn'
end subroutine
!specfit exited


	  
subroutine calculate_discornfall
use m_globals
use mpi
integer :: i,j,k,l,ind,ierr
integer :: info,fh
integer(kind=mpi_offset_kind) :: offset
safc = 0.0
safr = 0.0
samo = 0.0
do j=1,stnum
do i=1,2
spec = 0
spec(:)=disspectrum(i,:,j)
!print *,''
!      cornerfall format:
!      sig0,fcorn,falloffbest,fitrms, specfit(1->nfreq)	
if(getfc == 1) then
call  getcorn
else
call  getcorn2
end if
print *,myid,sig0,fcorn,falloffbest,fitrms,specfit(1)
!(optionally call getcorn)
asig0(i,j)=sig0
afcorn(i,j)=fcorn
afalloffbest(i,j)=falloffbest
afitrms(i,j)=fitrms
fitspectrum(i,1:nfreq,j)=specfit(1:nfreq)	
!weighted fc
!write(0,*) 'fc',myid,subarea(j),fcorn,falloffbest
safc(i)=safc(i)+subarea(j)*fcorn
safr(i)=safr(i)+subarea(j)*falloffbest
samo(i)=samo(i)+subarea(j)*10**sig0	  
end do
end do 

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/fcorn',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,afcorn(1,1),stnum*2,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/sig0',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,asig0(1,1),stnum*2,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)   

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/falloffbest',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,afalloffbest(1,1),stnum*2,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/fitrms',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,afitrms(1,1),stnum*2,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/fitspectrum',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2*nfreq
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,fitspectrum(1,1,1),stnum*nfreq*2,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'out/disspectrum',MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_info_null,fh,ierr)
offset = myid*stnum*2*nfreq
call mpi_file_set_view(fh,offset*sizereal,mpi_real,mpi_real,"native",mpi_info_null,ierr)
call mpi_file_write(fh,disspectrum(1,1,1),stnum*2*nfreq,mpi_real,mpi_status_ignore,ierr)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine
      
        

end module
