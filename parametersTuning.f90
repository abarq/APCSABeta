program parametersTuning
use APCSA
implicit none

integer   ::  io, idRun, idCase, idRunCase, Nmax, counter, versionCostFunction
integer   ::  printScreen, markovLengthDynamic
real(dp)  ::  neighborhoodFactor, deltaFactor, runTime, bestCostFn


character (len=*), parameter :: input = 'input.txt'
character (len=*), parameter :: output = 'output.txt'
character (len=*), parameter :: fmt1 = "(ES8.3e1, ES9.3e1, ES9.3e1, ES12.6e1, ES11.5e1, ES9.3e1, ES8.2e1, ES9.3e1, ES9.3e1)"

versionCostFunction = 1
printScreen = 0
markovLengthDynamic = 1

open(4,file=input)
open(5,file=output)

counter =1 
do 
   write(*,*) counter
    read(4,*, IOSTAT=io) idRun, idCase, idRunCase, semilla, Nmax, neighborhoodFactor, deltaFactor
	
	call startOptimization(neighborhoodFactor, Nmax, deltaFactor, runTime, bestCostFn, versionCostFunction, printScreen,&
& markovLengthDynamic)
	
    if (io.eq.0) then
        write(5,fmt1) float(idRun), dfloat(idCase), dfloat(idRunCase), dfloat(semilla), dfloat(Nmax), neighborhoodFactor, &
		& deltaFactor, runTime, bestCostFn  
	else
        exit
    end if
	counter = counter +1 
    IFRST = 0
end do

close(4)
close(5)


end program


