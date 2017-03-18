MODULE mod_inner_loop
USE mod_mpg
USE mod_huffman
IMPLICIT NONE
PRIVATE
PUBLIC :: inner_loop ! subroutine
CONTAINS
!------------------------------------------------------------------------------------------------
SUBROUTINE inner_loop(ibit, iblock_type, wk, i_mdct, side, itot_bits) ! ISO C.1.5.4.4
IMPLICIT NONE
INTEGER             , INTENT(IN    ) :: ibit, iblock_type
REAL (KIND = 8)     , INTENT(IN    ) :: wk(:)
TYPE (side_info_sub), INTENT(IN OUT) :: side
INTEGER             , INTENT(   OUT) :: i_mdct(:), itot_bits
INTEGER :: iq0, iqquant, iquantanf, ibigvalues, icount1, & 
           n0, n1, itab, iregion0, iregion1, ndiv
INTEGER, PARAMETER :: ndiv0 = 32, ndiv1 = 32
LOGICAL :: qfirst 
CALL calc_quantanf(wk, iquantanf) ! not wk 
ndiv    = ndiv0
iqquant = 0
iq0     = 0
qfirst  = .TRUE. 
DO WHILE (ndiv /= 0) ! obtain starting quantization step
 CALL quantization(iqquant, iquantanf, wk, i_mdct)
 IF ( qfirst .AND. iquantanf > -210 .AND. MAXVAL(i_mdct) <= 2**13 - 1) THEN ! 15 + 2**13 - 1
  iquantanf = iquantanf - 63 ! In ISO document 2.4.2.7 "big_values", the maximum absolute value is constrained to 8191 = 2^13 - 1.
  CYCLE                      ! However it is possible to use up to 2^13 - 1 + 15 = 8206.   
 ELSE                        ! ISO sample program dist10 uses 8206.
  qfirst = .FALSE.           ! Because ISO documet section 2 is normative, 2^13 - 1 is chosen in Uzura. 
 END IF                              
 IF (MAXVAL(i_mdct) > 2**13 - 1) THEN  ! bi-section method
  iq0 = iqquant                           
  iquantanf = iquantanf + ndiv 
 ELSE
!  ndiv = ndiv / 2
  ndiv = INT(ndiv * 0.618) ! Fibonacchi
  iqquant = iq0 + ndiv
 END IF 
END DO
ndiv    = ndiv1
iqquant = MAX(0, -iquantanf - 210)
iq0     = iqquant
DO  WHILE (ndiv > 0)
 CALL quantization(iqquant, iquantanf, wk, i_mdct)
 itot_bits = 0
 side%iglobal_gain = iqquant + iquantanf + 210 
 CALL divide(i_mdct, ibigvalues, icount1)
 side%ibig_values = ibigvalues 
 side%icount1     = icount1
 CALL count_one(i_mdct, ibigvalues, icount1, itab, itot_bits)
 side%icount1table_select = itab
 SELECT CASE (iblock_type) ! ISO 2.4.2.7 window_switching_flag
  CASE (0)
    CALL sub_divide(ibigvalues, iregion0, iregion1)  ! ISO 
    iregion0 = MIN(15, iregion0)                     ! ISO 2.4.2.7 region0_count 
    iregion1 = MIN( 7, iregion1, 19 - iregion0) 
    side%iregion0_count = iregion0     !  4 bits
    side%iregion1_count = iregion1     !  3 bits
    n0 = iscalefactorband_l(iregion0, 3) + 1 
    n1 = iscalefactorband_l(iregion0 + iregion1 + 1, 3) + 1
    CALL select_table2(i_mdct,      1, n0, itab, itot_bits)
    side%itable_select(1) = itab
    CALL select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
    side%itable_select(2) = itab
    CALL select_table2(i_mdct, n1 + 1, ibigvalues * 2, itab, itot_bits)
    side%itable_select(3) = itab
  CASE (10, 11, 30, 31) ! switching block
    iregion0 = 7                
    iregion1 = 36 
    side%iregion0_count = iregion0      !  4 bits
    side%iregion1_count = iregion1      !  3 bits
    n0 = iscalefactorband_l(iregion0, 3) + 1  
    n1 = ibigvalues * 2
    CALL select_table2(i_mdct,      1, n0, itab, itot_bits)
    side%itable_select(1) = itab
    CALL select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
    side%itable_select(2) = itab
  CASE (20) ! short 
    iregion0 = 8                
    iregion1 = 36 
    side%iregion0_count = iregion0        
    side%iregion1_count = iregion1        
    n0 = iscalefactorband_s( (iregion0 + 1) / 3, 2) * 3   
    n1 = ibigvalues * 2  
    CALL select_table2(i_mdct,      1, n0, itab, itot_bits)
    side%itable_select(1) = itab
    CALL select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
    side%itable_select(2) = itab
  CASE (21) ! mixed
    iregion0 = 7                
    iregion1 = 36 
    side%iregion0_count = iregion0        
    side%iregion1_count = iregion1        
    n0 = iscalefactorband_l(iregion0, 3) + 1
    n1 = ibigvalues * 2  
    CALL select_table2(i_mdct,      1, n0, itab, itot_bits)
    side%itable_select(1) = itab
    CALL select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
    side%itable_select(2) = itab
  CASE DEFAULT
    STOP ' error : subroutine inner_loop ' 
 END SELECT
 IF (itot_bits > ibit) THEN ! bi-section method: speeds up C.1.5.4.4.2
  iq0 = iqquant
  iqquant = iqquant + ndiv 
 ELSE
!  ndiv = ndiv / 2
  ndiv = INT(ndiv * 0.618) ! Fibonacchi
  iqquant = iq0 + ndiv 
 END IF 
END DO
WHERE (wk < 0.0d0) i_mdct = -i_mdct
RETURN
END SUBROUTINE inner_loop
!---------------------------------------------------------------------------------------------
SUBROUTINE calc_quantanf(wk, iquantanf)     ! ISO C.1.5.4.2.1
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: wk(:)
INTEGER, INTENT(OUT)         :: iquantanf
INTEGER :: i
REAL (KIND = 8) :: sfm, sum1, sum2, tmp
sum1 =  0.0d0
sum2 =  0.01d0 
DO i = 1, SIZE(wk) ! 576
 tmp = wk(i)**2.0d0
 IF (tmp > 0.0d0) sum1 = sum1 + LOG(tmp) 
 sum2 = sum2 + tmp
END DO
sfm = EXP( sum1 / 576.0d0) * 576.0d0 / sum2 
iquantanf = INT( 8.0d0 * LOG(sfm) ) 
RETURN
END SUBROUTINE calc_quantanf
!------------------------------------------------------------------------------------------------
SUBROUTINE quantization(iqquant, iquantanf, r_mdct, i_mdct) ! ISO C.1.5.4.4.1
IMPLICIT NONE
INTEGER        , INTENT(IN ) :: iqquant, iquantanf
REAL (KIND = 8), INTENT(IN ) :: r_mdct(:)
INTEGER        , INTENT(OUT) :: i_mdct(:)
REAL (KIND = 8) :: denom, tmp(576)
denom = 2.0d0 ** ( -REAL(iqquant + iquantanf, KIND = 8) / 4.0d0 )
tmp    = ABS(r_mdct) * denom
i_mdct = NINT( SQRT(tmp * SQRT(tmp)) - 0.0946d0 ) ! NINT(tmp**(3/4) - 0.0946)
RETURN
END SUBROUTINE quantization
!------------------------------------------------------------------------------------------------
SUBROUTINE divide(i_mdct, ibigvalues, icount1) ! ISO C.1.5.4.4.3, C.1.5.4.4.4
IMPLICIT NONE
INTEGER, INTENT(IN ) :: i_mdct(:)
INTEGER, INTENT(OUT) :: ibigvalues, icount1
INTEGER :: i, ibig, izero
DO i = 576, 110, -2
 izero = i
 IF (i_mdct(i) /= 0 .OR. i_mdct(i - 1) /= 0) EXIT 
END DO
DO i = izero, 110, -4
 ibig = i
 IF ( ABS(i_mdct(i    )) > 1 .OR. ABS(i_mdct(i - 1)) > 1 .OR. &        !  0, +1, - 1
      ABS(i_mdct(i - 2)) > 1 .OR. ABS(i_mdct(i - 3)) > 1        ) EXIT  
END DO
ibigvalues = ibig / 2 
icount1 = (izero - ibig) / 4 
RETURN
END SUBROUTINE divide 
!------------------------------------------------------------------------------------------------
SUBROUTINE count_one(i_mdct, ibigvalues, icount1, itab, isum) ! ISO 2.4.2.7, C.1.5.4.4.5
IMPLICIT NONE
INTEGER, INTENT(IN    ) :: i_mdct(:), ibigvalues, icount1
INTEGER, INTENT(   OUT) :: itab
INTEGER, INTENT(IN OUT) :: isum
INTEGER :: isum0, isum1, i, k1, k2, k3, k4
isum0 = 0
isum1 = 0
DO i = ibigvalues * 2 + 1, ibigvalues * 2 + icount1 * 4, 4 
 k1 = ABS(i_mdct(i    ))
 k2 = ABS(i_mdct(i + 1))
 k3 = ABS(i_mdct(i + 2))
 k4 = ABS(i_mdct(i + 3))
 isum0 = isum0 + huff_qa(k1, k2, k3, k4)%leng + k1 + k2 + k3 + k4 
 isum1 = isum1 + huff_qb(k1, k2, k3, k4)%leng + k1 + k2 + k3 + k4
END DO
IF (isum0  <= isum1) THEN
 itab = 0 ! use Table A
 isum = isum + isum0
ELSE
 itab = 1 ! use Table B
 isum = isum + isum1 
END IF
RETURN
END SUBROUTINE count_one
!------------------------------------------------------------------------------------------------
SUBROUTINE sub_divide(ibigvalues, iregion0, iregion1) ! ISO C.1.5.4.4.6
IMPLICIT NONE
INTEGER, INTENT(IN ) :: ibigvalues
INTEGER, INTENT(OUT) :: iregion0, iregion1
INTEGER :: n0, n1, i
n0 = 2 * ibigvalues * r0  
n1 = 2 * ibigvalues * r1  
! division suggested in ISO document C.1.5.4.4.6 is 1/3 : 5/12 : 1/4     
DO i = 0, 20
 IF ( n0 >= iscalefactorband_l(i, 3) ) iregion0 = MIN( 15, MAX( 0, i               ) ) 
 IF ( n1 >= iscalefactorband_l(i, 3) ) iregion1 = MIN(  7, MAX( 0, i - iregion0 - 1) )
END DO
RETURN
END SUBROUTINE sub_divide
!------------------------------------------------------------------------------------------------
SUBROUTINE select_table2(ix, n0, n1, itab, isum) ! C.1.5.4.4.5
IMPLICIT NONE
INTEGER, INTENT(IN    ) :: ix(:), n0, n1
INTEGER, INTENT(   OUT) :: itab
INTEGER, INTENT(IN OUT) :: isum
INTEGER :: imax, isum0, isum1, isum2, isum3, i, itab1, itab2
imax = MAXVAL( ABS(ix(n0:n1)) )
itab = -999
IF (imax <= 15) THEN 
 DO i = 13, 0, -1
  IF (imax <= huff(i)%nmax) itab = i
 END DO
 SELECT CASE (itab)
  CASE (0)
   itab = 0
  CASE (1)
   itab = 1
  CASE (2) 
   CALL count_bits(2, ix, n0, n1, isum1)
   CALL count_bits(3, ix, n0, n1, isum2)
   IF (isum1 < isum2) THEN 
    itab = 2
   ELSE
    itab = 3
   END IF
  CASE (5)
   CALL count_bits(5, ix, n0, n1, isum1)
   CALL count_bits(6, ix, n0, n1, isum2)
   IF (isum1 < isum2) THEN 
    itab = 5
   ELSE
    itab = 6
   END IF
  CASE (7)
   CALL count_bits(7, ix, n0, n1, isum1)
   CALL count_bits(8, ix, n0, n1, isum2)
   CALL count_bits(9, ix, n0, n1, isum3)
   isum0 = MIN(isum1, isum2, isum3)
   IF      (isum1 == isum0) THEN 
    itab = 7
   ELSE IF (isum2 == isum0) THEN
    itab = 8
   ELSE
    itab = 9
   END IF
  CASE (10)
   CALL count_bits(10, ix, n0, n1, isum1)
   CALL count_bits(11, ix, n0, n1, isum2)
   CALL count_bits(12, ix, n0, n1, isum3)
   isum0 = MIN(isum1, isum2, isum3)
   IF      (isum1 == isum0) THEN 
    itab = 10
   ELSE IF (isum2 == isum0) THEN
    itab = 11
   ELSE
    itab = 12
   END IF
  CASE (13)
   CALL count_bits(13, ix, n0, n1, isum1)
   CALL count_bits(15, ix, n0, n1, isum2)
   IF (isum1 < isum2) THEN 
    itab = 13
   ELSE
    itab = 15
   END IF
  CASE DEFAULT
   STOP 'something wrong : select_table2 : <= 15'
 END SELECT    
ELSE
 DO i = 16, 23
  IF (imax <= huff(i)%linmax + 15) THEN 
   itab1 = i
   EXIT
  END IF
 END DO
 DO i = 24, 31
  IF (imax <= huff(i)%linmax + 15) THEN 
   itab2 = i
   EXIT
  END IF
 END DO
 CALL count_bits(itab1, ix, n0, n1, isum1)
 CALL count_bits(itab2, ix, n0, n1, isum2)
 IF (isum1 < isum2) THEN 
  itab = itab1
 ELSE
  itab = itab2
 END IF
END IF
CALL count_bits(itab, ix, n0, n1, isum0)
isum = isum + isum0
RETURN
END SUBROUTINE select_table2
!------------------------------------------------------------------------------------------------
SUBROUTINE count_bits(nhuff, ix, n0, n1, isum) ! C.1.5.4.4.8
IMPLICIT NONE
INTEGER   , INTENT(IN ) :: nhuff, ix(:), n0, n1
INTEGER   , INTENT(OUT) :: isum
INTEGER :: i, k1, k2
isum = 0
IF (nhuff ==  0) RETURN
IF (nhuff <= 15) THEN 
 DO i = n0, n1, 2
   k1 = ABS( ix(i    ) )
   k2 = ABS( ix(i + 1) )
  isum = isum + huff(nhuff)%leng(k1, k2)
  IF (k1 /= 0) isum = isum + 1 ! one bit for sign
  IF (k2 /= 0) isum = isum + 1 ! one bit for sign
 END DO
ELSE
 DO i = n0, n1, 2
  k1 = ABS( ix(i    ) )
  k2 = ABS( ix(i + 1) )
  IF (k1 > 14) THEN 
   k1 = 15
   isum = isum + huff(nhuff)%linbits
  END IF
  IF (k2 > 14) THEN 
   k2 = 15
   isum = isum + huff(nhuff)%linbits
  END IF
  isum = isum + huff(nhuff)%leng(k1, k2)
  IF (ix(i    ) /= 0) isum = isum + 1 ! one bit for sign
  IF (ix(i + 1) /= 0) isum = isum + 1 ! one bit for sign
 END DO
END IF
RETURN
END SUBROUTINE count_bits
!------------------------------------------------------------------------------------------------
END MODULE mod_inner_loop
