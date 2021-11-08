module APCSA
use randomNum
use inputData
use instanceProblem
implicit none


integer, parameter :: maxTempCycles = 500
integer, parameter :: N = 200000
integer, parameter :: nMov = 30
real(dp)           :: tolTemp = 1.0e-4_dp
real(dp)           :: initialAcceptance = 0.90_dp
real(dp)           :: finalAcceptance = 1.0e-6_dp
real(dp)           :: eta = 0.005_dp
real(dp)           :: alfa

contains

subroutine startOptimization(neighborhoodFactor, Nmax, deltaFactor, runTime, bestCostFn, versionCostFunction, printScreen,&
& markovLengthDynamic)
implicit none
integer,  intent (IN)              :: Nmax
real(dp), intent (IN)              :: neighborhoodFactor
real(dp), intent (IN)              :: deltaFactor
real(dp), intent (OUT)             :: runTime
real(dp), intent (OUT)             :: bestCostFn

real(dp), dimension(nParam)        :: L, U, deltaP, stepP, x, xOld,freqCam, xMem
real(dp), dimension(maxTempCycles) :: fnAverage, EminArray, acceptanceRateTemp, bestCostFnTemp
real(dp), dimension(N)             :: energiaMovAcep, historicAcceptance, beta
real(dp), dimension(nSamples)      :: energy, epsUnoE, epsDosE
real(dp)                           :: T, deltaE, E,  promEAcept, prev_Temp
real(dp)                           :: xn, EOld, probAcep, sigma, desviacion
real(dp)                           :: start, finale, EMem, T_APCSA, deltaFnProm, uncertainty
real(dp)                           :: tolEqTer, EMin, Tini
integer                            :: m, i, nMov_Accepted, nMov_Rejected, flagTole, versionCostFunction,printScreen
integer                            :: markovLengthDynamic
logical                            :: isMovAccepted, flag, solidification

call cpu_time(start)


x = 0.0_dp
xOld = 0.0_dp
tolEqTer = 5.0e-2_dp
deltaFnProm = 0.0_dp
fnAverage = 0.0_dp
beta = 0.0_dp
sigma = sqrt(-maxTempCycles**2/(2*log(finalAcceptance/initialAcceptance)))

call readData(energy, epsUnoE, epsDosE)
call getBounds(U,L)

deltaP = U-L
stepP = deltaP/neighborhoodFactor


do i=1,nParam
    call randNum(xn)
    xOld(i) = L(i) + xn*deltaP(i)
end do


call get_Init_Temp(Tini, L, U, energy, epsUnoE, epsDosE,versionCostFunction)
T = Tini
T_APCSA = Tini


call fn(EOld, uncertainty,  xold, energy, epsUnoE, epsDosE,0,versionCostFunction)

EMem = EOld
xMem = xOld

prev_Temp = T_APCSA



solidification = .true.
flagTole = 1
freqCam = 1.00_dp

m = 1
bestCostFnTemp(m) = EOld

do while ((m.le.maxTempCycles).and.solidification)

    nMov_Accepted = 0
    nMov_Rejected = 0
    flag = .true.
    EMin = 1.0e10_dp
    i = 1
    
    !Metropolis Algorithm
    do while (flag.and.(i.le.Nmax))
                
        call makeTransition(xOld, x, L, U, stepP, freqCam)
        call fn(E, uncertainty, x, energy, epsUnoE, epsDosE,0,versionCostFunction)
        deltaE = E-EOld

        if (deltaE.lt.0.0_dp) then
                
            EOld = E
            xOld = x
            call updateRegister(xMem, EMem, Emin,  xOld, EOld)
            nMov_Accepted = nMov_Accepted + 1
            energiaMovAcep(nMov_Accepted) = EOld
            historicAcceptance(nMov_Accepted) = 1.0_dp
            isMovAccepted = .true.
            beta(nMov_Accepted) = EOld/prev_Temp

        else 
            probAcep = exp(-deltaE/T_APCSA)
            call randNum(xn)
                        
            if (xn.lt.probAcep) then
                EOld = E
                xOld = x
                call updateRegister(xMem, EMem, Emin, xOld, EOld)
                nMov_Accepted = nMov_Accepted + 1
                energiaMovAcep(nMov_Accepted) = EOld
                historicAcceptance(nMov_Accepted) = probAcep
                isMovAccepted = .true.
                beta(nMov_Accepted) = EOld/prev_Temp
            else
                nMov_Rejected = nMov_Rejected + 1
                isMovAccepted = .false.                                
            end if
                        
        end if

        if (isMovAccepted.and.(markovLengthDynamic.eq.1)) then
            call evaluate_equilibriumCondition(nMov_Accepted, beta, tolEqTer, flag)
        end if         
        
        i = i+1
 
    end do
        

    
    bestCostFnTemp(m) = EMem        
    if (nMov_Accepted.eq.0) then

        solidification = .false.

    else
        EminArray(m) = Emin
        acceptanceRateTemp(m) = float(nMov_Accepted)/(float(nMov_Accepted+nMov_Rejected))
        call calcularDesviacion(energiaMovAcep,nMov_Accepted, desviacion, promEAcept)

        fnAverage(m) = promEAcept

        prev_Temp = T_APCSA
        T_APCSA = -desviacion/(dlog(initialAcceptance) - ((dfloat(m)**2)/(2*sigma**2)))
            
            
        T = alfa*T    
        call evaluate_solidification_condition(m, EminArray, solidification)
    call fn(EOld, uncertainty, xOld, energy, epsUnoE, epsDosE, 1,versionCostFunction)
    
    if (m.ge.2) then
        deltaFnProm = abs(fnAverage(m-1)-fnAverage(m)) 
        
        if (deltaFnProm.gt.uncertainty) then
            call update_stepSize(x,stepP, deltaFactor)
        else
            if (flagTole.eq.1) then
                tolEqTer = 10_dp*tolEqTer
                flagTole = 0
            end if
                
        end if
    end if
    
    if (printScreen.eq.1) then
        
    	write(*,*) m, Emin, T_APCSA,  deltaFnProm, uncertainty, i
    end if
    !write(*,*) Emin, acceptanceRateTemp(m)
    
    m = m+1
    end if
          
end do


bestCostFn = EMem


call cpu_time(finale)
runTime = finale-start

!write(*,*) xMem

end subroutine




subroutine getBounds(U,L)
use instanceProblem
implicit none
real(dp), intent (OUT) :: L(:)
real(dp), intent (OUT) :: U(:)


L(1) = 1.0_dp
U(1) = 10_dp

L(2) = 0.0_dp
U(2) = 1.0_dp

L(3) = 0.0_dp
U(3) = 0.1_dp

L(4) = 0.0_dp
U(4) = 1.0_dp

L(5) = 0.5_dp
U(5) = 1.5_dp

!L(6) = 0.5_dp
!U(6) = 2.5_dp

!L(7) = 2.0_dp
!U(7) = 3.0_dp

L(6) = 0.0_dp
U(6) = 1.0_dp

L(7) = 0.0_dp
U(7) = 1.0_dp

L(8) = 0.0_dp
U(8) = 5.0_dp

L(9) = 0.0_dp
U(9) = 5.0_dp

L(10) = 0.40_dp
U(10) = 0.42_dp

L(11) = 0.0_dp
U(11) = 5.0_dp

L(12) = 0.0_dp
U(12) = 5.0_dp

L(13) = 1.214_dp
U(13) = 1.226_dp

L(14) = 0.0_dp
U(14) = 5.0_dp

L(15) = 0.0_dp
U(15) = 5.0_dp

L(16) = 2.7_dp
U(16) = 2.9_dp

L(17) = 0.0_dp
U(17) = 5.0_dp

L(18) = 0.0_dp
U(18) = 5.0_dp

L(19) = 4.9_dp
U(19) = 5.0_dp

L(20) = 0.0_dp
U(20) = 5.0_dp

L(21) = 0.0_dp
U(21) = 5.0_dp

L(22) = 11.0_dp
U(22) = 11.1_dp

L(23) = 0.1_dp
U(23) = 0.4_dp



end

subroutine solWilliam(xold)
implicit none
real(dp) :: xold(:)

xold(1) = 5.145_dp
xold(2) = 0.552_dp
xold(3) = 0.009_dp
xold(4) = 0.244_dp
xold(5) = 0.965_dp
xold(6) = 1.287_dp
xold(7) = 2.265_dp
xold(8) = 0.344_dp
xold(9) = 0.666_dp
xold(10) = 0.411_dp
xold(11) = 0.710_dp
xold(12) = 1.969_dp
xold(13) = 1.216_dp
xold(14) = 0.372_dp
xold(15) = 3.423_dp
xold(16) = 2.777_dp
xold(17) = 0.358_dp
xold(18) = 3.093_dp
xold(19) = 4.974_dp
xold(20) = 0.999_dp
xold(21) = 4.352_dp
xold(22) = 11.06_dp
xold(23) = 0.37_dp

end subroutine


subroutine makeTransition(xOld, x, L, U, stepP, freqCam)
use instanceProblem
use randomNum
implicit none
real(dp), intent (IN)  :: xOld(:)
real(dp), intent (IN)  :: L(:)
real(dp), intent (IN)  :: U(:)
real(dp), intent (IN)  :: stepP(:)
real(dp), intent (IN)  :: freqCam(:)
real(dp), intent (OUT) :: x(:)

real(dp)  :: numAlea, xn, dif
logical :: flag
integer :: i, j

flag = .true.
j = 1
do while (flag)
        
    do i=1, nParam
        
        xn = 0.01_dp

        if (xn.le.freqCam(i)) then
        
        
            call randNum(xn)
            numAlea = -1.0_dp + 2.0_dp*xn
            x(i) = xOld(i) +  numAlea*stepP(i)
                        
            if (x(i).gt.U(i)) then
                dif = x(i) - U(i)
                x(i) = U(i)-dif
            end if
                        
            if (x(i).lt.L(i)) then
                dif = abs(x(i)-L(i))
                x(i) = L(i)+dif
            end if
            flag = .false.
        end if
                
                
    end do
    
    
end do


end subroutine



subroutine update_stepSize(x,stepP,deltaFactor)
implicit none

real(dp), intent (INOUT) :: stepP(:)
real(dp), intent (IN)    :: x(:)
real(dp), intent (IN)    :: deltaFactor

integer      :: i
real(dp)     :: stepFactor, dummy1


do i=1,nParam
        
    dummy1 = stepP(i)/x(i)
        
    stepFactor = 1/deltaFactor
        
    if (dummy1.gt.eta) then
        stepP(i) = stepP(i)*stepFactor
    end if
        
end do

end subroutine

subroutine updateRegister(xMem, EMem, Emin,  xOld, EOld)
use instanceProblem
implicit none
real(dp), dimension(nParam) :: xMem, xOld
real(dp) :: EOld, EMem, Emin

    if (EOld.lt.EMem) then
        EMem = EOld
        xMem = xold
    end if

    if (EOld.lt.EMin) then
        EMin = EOld
    end if


end subroutine


subroutine calcularDesviacion(energiaMovAcep,nMov_Accepted, desviacion,promEAcept)
implicit none
real(dp) :: energiaMovAcep(N), suma, prom, desviacion, promEAcept
integer :: nMov_Accepted, i

suma = 0.0_dp
do i=1,nMov_Accepted
        suma = suma + energiaMovAcep(i)
end do

prom = suma/float(nMov_Accepted)
promEAcept = prom
suma = 0.0_dp

do i=1,nMov_Accepted
        suma = suma + abs(energiaMovAcep(i)-prom)
end do

desviacion = suma/float(nMov_Accepted)
end subroutine


subroutine updateMoveFrequencies(xOld, L, U, freqCam, stepP, energy, epsUnoE, epsDosE)
use inputData
implicit none
real(dp), dimension(nParam) :: xOld, freqCam, L, U, stepP, x, z, deltaEk
real(dp), dimension(nSamples) :: energy, epsUnoE, epsDosE
real(dp) :: yavg, dyavg, E, dekMax, xn, rn, dif, sumaFreq, uncertainty
integer :: i, j, versionCostFunction

dekMax = 0.0_dp
z = 0.0_dp

do i = 1,nParam
    yavg = 0.0_dp
    dyavg = 0.0_dp
    x = xOld
    do j= 1,nParam
        
        call randNum(xn)
        rn = -1.0_dp + 2.0_dp*xn
        x(i) = xOld(i) + rn*stepP(i)

        if (x(i).gt.U(i)) then
            dif = x(i) - U(i)
            x(i) = U(i)-dif
        end if
                    
        if (x(i).lt.L(i)) then
            dif = abs(x(i)-L(i))
            x(i) = L(i)+dif
        end if

        
        call fn(E, uncertainty, x, energy, epsUnoE, epsDosE,0,versionCostFunction)
        !write(*,*) E
        z(j) = E 
        
    end do
    
    yavg = sum(z)/dfloat(nParam)
    
    dyavg = sum(abs(z-yavg))/dfloat(nParam)

    deltaEk(i) = dyavg
    
    if (deltaEk(i).gt.dekMax) then
        dekMax = deltaEk(i)
    end if
    
end do

do i=1,nParam
    freqCam(i) = 0.80_dp*deltaEk(i)/dekMax
    
end do

sumaFreq = sum(freqCam)


do i=1,nParam
    freqCam(i) = freqCam(i)/sumaFreq
    write(*,*) freqCam(i)
end do



end subroutine


subroutine get_Init_Temp(Tini, L, U, energy, epsUnoE, epsDosE,versionCostFunction)
use randomNum
implicit none
real(dp), intent (IN)  :: L(:)
real(dp), intent (IN)  :: U(:)
real(dp), intent (IN)  :: energy(:)
real(dp), intent (IN)  :: epsUnoE(:)
real(dp), intent (IN)  :: epsDosE(:)
real(dp), intent (OUT) :: Tini

real(dp), dimension(nParam) ::  x
real(dp), dimension(nMov)   :: energias
real(dp)                    :: xn, E, suma, average, deviation, uncertainty
integer                     :: i, j, versionCostFunction

x=0.0_dp
energias=0.0_dp
suma = 0.0_dp
average = 0.0_dp
deviation = 0.0_dp

do i=1,nMov

    do j=1,nParam
        
        call randNum(xn)
        x(j) = L(j) + xn*(U(j)-L(j))
                
    end do
     
    call fn(E, uncertainty, x, energy, epsUnoE, epsDosE,0,versionCostFunction)
	energias(i) = E
        
end do

average = sum(energias)/dfloat(nMov)

deviation = sum(abs(energias-average))/dfloat(nMov)

Tini=-deviation/dlog(initialAcceptance)

end subroutine

subroutine evaluate_solidification_condition(m, EminArray, solidification)
implicit none

real(dp), intent (IN)   :: EminArray(:)
integer,  intent (IN)   :: m
logical,  intent (OUT)  :: solidification

real(dp)                :: conver, conver2, conver3

if (m.ge.4) then
    conver =  abs((EminArray(m)-EminArray(m-1))/EminArray(m)) 
    conver2 = abs((EminArray(m)-EminArray(m-2))/EminArray(m))
    conver3 = abs((EminArray(m)-EminArray(m-3))/EminArray(m))
    
    if ((conver.lt.tolTemp).and.(conver2.lt.tolTemp).and.(conver3.lt.tolTemp)) then 
        solidification = .false.
    end if
end if  

end subroutine

subroutine evaluate_equilibriumCondition(nMov_Accepted, beta, tolEqTer, flag)
implicit none
real(dp), intent (IN)   :: beta(:)
real(dp), intent (IN)   :: tolEqTer
integer,  intent (IN)   :: nMov_Accepted
logical,  intent (INOUT)  :: flag

real(dp) :: d0embeta, nn, told
integer  :: j

d0embeta = 0.0_dp
do j=1,nMov_Accepted-1
    d0embeta = d0embeta + dexp(beta(nMov_Accepted)-beta(j))
end do
nn = float(nMov_Accepted)    
told = 1.0_dp - ((nn-1.0_dp)/nn)*(1.0_dp+1.0_dp/d0embeta)
told = dabs(told)
if (told.le.tolEqTer) then
    flag = .false.
end if

end subroutine



end module APCSA
