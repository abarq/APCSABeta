module randomNum
use kinds
implicit none

integer :: semilla 
integer, parameter :: MPLIER=16807,MODLUS=2147483647,MOBYMP=127773
integer, parameter :: MOMDMP=2836
integer :: IFRST = 0

contains

subroutine randNum(xn)
implicit none
integer hvlue,lvlue,testv,nextn
real(dp) xn
save nextn
if (ifrst.eq.0) then
   nextn=semilla
   ifrst=1
endif
hvlue=nextn/mobymp
lvlue=mod(nextn, mobymp)
testv=mplier*lvlue-momdmp*hvlue
if(testv.gt.0) then
 nextn=testv
else
 nextn=testv+modlus
endif
xn=dfloat(nextn)/dfloat(modlus)
return

end subroutine randNum
end module randomNum