#rm timeseries spectrum displacement disspectrum cornerfall specfit
#gfortran   waveform2_revised.f -O3 -I/export/home/wang/local/fftwgnu/include  -L/export/home/wang/local/fftwgnu/lib -lfftw3 -lm -o waveform2_revised #-vec-report2#  -assume byterecl
#gfortran -Wall -fimplicit-none waveform.f -O3 -o waveform -lfftw3 -lm
#ftn waveform2.f -lm -I/opt/fftw/3.3.4.0/x86_64/include -L/opt/fftw/3.3.4.0/x86_64/lib  -lfftw3
#ifort waveform.f -O3 -assume byterecl -lfftw3 -lm -o waveform
#ifort  -I/usr/local/include fftw.f90 globals.f90  arrays.f90 util.f90 \
#parameters.f90 timeseries.f90 spectra.f90 waveform5.f90  -O3 -assume byterecl -lfftw3 -lm \
#-o ../example/waveform
#rm timeseries displacement specfit disspectrum cornerfall 
#gfortran -I/usr/local/include fftw.f90 globals.f90  arrays.f90 util.f90 \
#parameters.f90 timeseries.f90 spectra.f90 waveform5.f90  -O3 -lfftw3 -lm -o ../example/waveform
mpif90 -I/export/home/wang/local/fftwgnu/include  fftw.f90 globals.f90  arrays.f90 util.f90 \
parameters.f90 timeseries.f90 spectra.f90 waveform6.f90  -O3 -L/export/home/wang/local/fftwgnu/lib -lfftw3 -lm -o farspectra
