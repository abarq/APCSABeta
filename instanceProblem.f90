module instanceProblem
use randomNum
use inputData
use kinds
implicit none
integer, parameter :: nParam = 24
integer            :: nOs = 5
contains

subroutine fn(epsilonSumaTotal, uncertainty, x, energy, epsUnoE, epsDosE, flag,versionCostFunction)
implicit none
real(dp), intent (IN)  :: x(:)
real(dp), intent (IN)  :: energy(:)
real(dp), intent (IN)  :: epsUnoE(:)
real(dp), intent (IN)  :: epsDosE(:)
integer,  intent (IN)  :: flag
real(dp), intent (OUT) :: epsilonSumaTotal
real(dp), intent (OUT) :: uncertainty

real(dp)                      :: omegaPe, f0, gammaOE, beta, eta,epshf, omegaPe2,OPh2, gammaE, gammaH 
real(dp)                      :: omega, omega2, gammaE2, gammaH2, deps1, deps2
real(dp)                      :: dummy1, dummy2, dummy3, dummy4, xi
real(dp), dimension(nSamples) :: drude1, drude2, lorentz1, lorentz2, epsilon1, epsilonTotal, epsilon2
real(dp), dimension(nSamples) :: difEpsUno, difEpsDos, array1, array2
integer                       :: i, k,j, kk, nDrudeParameters, casoFuncionCosto,versionCostFunction
complex*16                    :: II, eps

deps1 = 0.01_dp
deps2 = 0.01_dp

nDrudeParameters = 8
omegaPe = x(1)
f0 = x(2)
gammaOE = x(3)
beta = x(4)
eta = x(5)
epshf = 1/f0
omegaPe2 = omegaPe*omegaPe
OPh2 = omegaPe2*beta
xi = x(8)

dummy1 = x(6) 
dummy2 = x(7)
dummy3 = 1/x(nParam)


select case(versionCostFunction)

    case(1)

    do i=1,nSamples
        omega = energy(i)
        omega2 = omega*omega
        gammaE = gammaOE*(1+(omega*dummy1)**2)
        gammaE2 = gammaE*gammaE
        gammaH = eta*gammaOE*(1+(omega*dummy2)**2)
        gammaH2 = gammaH*gammaH
        drude1(i) = epshf - omegaPe2/(omega2 + gammaE2) - OPh2/(omega2 + gammaH2)
        drude2(i) = (gammaE*omegaPe2)/(omega*(omega2 + gammaE2)) + (OPh2*gammaH)/(omega*(omega2 + gammaH2))

        lorentz1(i) = 0.0_dp
        lorentz2(i) = 0.0_dp

        k = 1*3
        kk = nDrudeParameters + k
        lorentz1(i) = lorentz1(i) + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
        lorentz2(i) = lorentz2(i) + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)    

        k = 2*3
        kk = nDrudeParameters + k
        lorentz1(i) = lorentz1(i) + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
        lorentz2(i) = lorentz2(i) + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)    
		
        k = 3*3
        kk = nDrudeParameters + k
        lorentz1(i) = lorentz1(i)  + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
        lorentz2(i) = lorentz2(i)  + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)    
		
        k = 4*3
        kk = nDrudeParameters + k
        lorentz1(i) = lorentz1(i) + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
        lorentz2(i) = lorentz2(i) + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)    
		
        k = 5*3
        kk = nDrudeParameters + k
        lorentz1(i) = lorentz1(i) + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
        lorentz2(i) = lorentz2(i) + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)    
		

    end do


    epsilon1 = drude1 + dummy3*lorentz1
    epsilon2 = drude2 + dummy3*lorentz2
	
	case(2)
	
    do i=1,nSamples
        omega = energy(i)
        omega2 = omega*omega
        gammaE = gammaOE*(1+(omega*dummy1)**2)
        gammaE2 = gammaE*gammaE
        gammaH = eta*gammaOE*(1+(omega*dummy2)**2)
        gammaH2 = gammaH*gammaH
        drude1(i) = epshf - omegaPe2/(omega2 + gammaE2) - OPh2/(omega2 + gammaH2)
        drude2(i) = (gammaE*omegaPe2)/(omega*(omega2 + gammaE2)) + (OPh2*gammaH)/(omega*(omega2 + gammaH2))

        lorentz1(i) = 0.0_dp
        lorentz2(i) = 0.0_dp
		
		do j=1,nOs
		    k = j*3
		    kk = nDrudeParameters + k
			lorentz1(i) = lorentz1(i) + ((x(kk-2)*omegaPe2)*(x(kk)**2-omega2))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2)
            lorentz2(i) = lorentz2(i) + ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega2)**2 + (omega*x(kk-1))**2) 			
		end do

    end do


    epsilon1 = drude1 + dummy3*lorentz1
    epsilon2 = drude2 + dummy3*lorentz2	
	
    case(3)
	
    II = dcmplx(0.0D0,1.D0)
    do i=1,nSamples
        omega = energy(i)
        gammaE = gammaOE*(1+(omega*dummy1)**2)
        gammaH = eta*gammaOE*(1+(omega*dummy2)**2)
        eps = epshf - 1/(omega)*(omegaPe2/(omega + II*gammaE) + OPh2/(omega + II*gammaH))
	
        do j=1,nOs
            k = j*3
            kk = nDrudeParameters + k
            eps = eps +  dummy3*(((x(kk-2)*omegaPe2))/(x(kk)**2 - omega**2 - II*(omega*x(kk-1)))); 
        end do
		
        epsilon1(i) = dreal(eps)
        epsilon2(i) = dimag(eps)
    end do

end select

dummy4 = 1/dble(2*nSamples-nParam-1);

difEpsUno = ((epsilon1 - epsUnoE)/epsUnoE)**2
difEpsDos = ((epsilon2 - epsDosE)/epsDosE)**2

epsilonTotal = difEpsUno + difEpsDos

epsilonSumaTotal = sum(epsilonTotal)*dummy4


uncertainty = 0.0_dp
if (flag.eq.1) then

    array1 = dabs(dabs(epsilon1/(epsUnoE))-1.0_dp)
    array2 = dabs(dabs(epsilon2/(epsDosE))-1.0_dp)
    uncertainty = sum(2.0_dp*array1*DABS(epsilon1)*deps1/(epsUnoE**2)) + sum(2.0_dp*array2*DABS(epsilon2)*deps2/(epsDosE**2))
    uncertainty = uncertainty*dummy4

end if


end subroutine


end module instanceProblem
