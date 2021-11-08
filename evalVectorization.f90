program evalVectorization
use APCSA
implicit none

integer   ::  idCase, Nmax, counter, nRep, i, versionCostFunction, printScreen, markovLengthDynamic
real(dp)  ::  neighborhoodFactor, deltaFactor, runTime, bestCostFn


character (len=*), parameter :: output = 'outputEvalVec.txt'
character (len=*), parameter :: fmt1 = "(ES9.3e1, ES9.3e1, ES9.3e1, ES9.3e1)"
CHARACTER(100)               :: num1char
CHARACTER(100)               :: num2char

   

call GET_COMMAND_ARGUMENT(1,num1char) 
call GET_COMMAND_ARGUMENT(2,num2char) 

read(num1char,*) versionCostFunction  
read(num2char,*) idCase  


semilla = 1811695
Nmax = 15000
neighborhoodFactor =  1500.00
deltaFactor = 1.00
nRep = 1
markovLengthDynamic = 0
printScreen = 1

open(5,file=output, access='append')


do i=1,nRep   
	
	call startOptimization(neighborhoodFactor, Nmax, deltaFactor, runTime, bestCostFn,versionCostFunction, printScreen,&
	& markovLengthDynamic)

        write(5,fmt1) dfloat(idCase), dfloat(versionCostFunction), runTime, bestCostFn  
    IFRST = 0
end do


close(5)


end program


