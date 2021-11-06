module inputData
use kinds
implicit none
integer, parameter           :: nSamples = 200
character (len=*), parameter :: expData = 'PDWB200'
contains

subroutine readData(energy, epsUnoE, epsDosE)
implicit none
real(dp), dimension(nSamples) :: energy, epsUnoE, epsDosE, coefN, coefK
integer                       :: j


	! open(3,file=expData)

	! do j=1,nSamples
		! read(3,*) energy(j), coefN(j), coefK(j)
	! end do

	! close(3)

	! epsUnoE = coefN**2 - coefK**2
	! epsDosE = 2*coefN*coefK
	
    open(3,file=expData)

    do j=1,nSamples
        read(3,*) energy(j), epsUnoE(j), epsDosE(j)
    end do

    close(3)



end subroutine
end module inputData