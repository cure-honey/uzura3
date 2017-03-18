MODULE mod_crc
IMPLICIT NONE
PRIVATE
PUBLIC :: crc16, crc16dim
INTEGER, PARAMETER :: igenerator = 32773 ! Z'8005' !B'1000000000000101'  ! x^16+x^15+x^2+1
CONTAINS
!--------------------------------------------------------------------------
SUBROUTINE crc16(n, in, icrc)         ! ISO 2.4.3.1, Table A.9, Table B.5 
IMPLICIT NONE
INTEGER, INTENT(IN    ) :: n, in
INTEGER, INTENT(IN OUT) :: icrc
INTEGER :: j, ibit1, ibit2
DO j = n - 1, 0, -1
 ibit1 = IBITS(in  ,  j, 1)           ! jth bit                   bit[31......0]
 ibit2 = IBITS(icrc, 15, 1)           ! sign bit of 16bit crc
 icrc  = ISHFT(IBITS(icrc, 0, 15), 1) ! shift up 1bit 16bit crc
 IF (IEOR(ibit1, ibit2) == 1) icrc = IEOR(icrc, igenerator)
END DO
RETURN
END SUBROUTINE crc16
!--------------------------------------------------------------------------
SUBROUTINE crc16dim(n, in, icrc)
IMPLICIT NONE
INTEGER, INTENT(IN    ) :: n, in(:)
INTEGER, INTENT(IN OUT) :: icrc
INTEGER :: i
DO i = 1, SIZE(in)
 CALL crc16(n, in(i), icrc)
END DO
RETURN
END SUBROUTINE crc16dim
!--------------------------------------------------------------------------
END MODULE mod_crc

