MODULE mod_mdct
IMPLICIT NONE
PRIVATE
PUBLIC :: mdct_initialize, sub_mdct
REAL (KIND = 8):: pi
INTEGER :: indx9(9)
COMPLEX (KIND = 8), SAVE:: omega9(0:8), omega9s(9), sqrt3_2, omega3(0:2), omega3s(3)
REAL (KIND = 8), PARAMETER :: c(0:7) = (/ -0.6000d0, -0.5350d0, -0.3300d0, -0.1850d0, -0.0950d0, -0.0410d0, -0.0142d0, -0.0037d0/)
REAL (KIND = 8), SAVE :: ca(0:7), cs(0:7), window_n(36), window_s(12), window_start(36), window_stop(36)
REAL (KIND = 8), SAVE, ALLOCATABLE :: subbuff(:, :, :)
CONTAINS
!---------------------------------------------------------------------------------------------
SUBROUTINE mdct_initialize()
IMPLICIT NONE
INTEGER :: i
pi = 4.0d0 * ATAN(1.0d0)
subbuff = 0.0d0
CALL fft9_initialize()
CALL fft3_initialize()
! ISO Table B.9 coefficients for aliasing reduction:
DO i = 0, 7
 cs(i) =  SQRT( 1.0d0    / ( 1.0d0 + c(i)**2 ) )
 ca(i) = -SQRT( c(i)**2  / ( 1.0d0 + c(i)**2 ) )
END DO
! ISO 2.4.3.4.10.3 Windowing, C.1.5.3.3
! normal window
DO i = 1, 36 
 window_n(i) = SIN( REAL(2 * i - 1, KIND = 8) / 72.0d0 * pi)
END DO
! short window
DO i = 1, 12 
 window_s(i) = SIN( REAL(2 * i - 1, KIND = 8) / 24.0d0 * pi)
END DO
! start / stop window
DO i = 1, 18 
 window_start(i     ) = SIN( REAL(2 *  i       - 1, KIND = 8) / 72.0d0 * pi)
 window_stop (i + 18) = SIN( REAL(2 * (i + 18) - 1, KIND = 8) / 72.0d0 * pi)
END DO
DO i = 1, 6
 window_start(i + 24) = SIN( REAL(2 * (i + 6) - 1, KIND = 8) / 24.0d0 * pi)
 window_stop (i +  6) = SIN( REAL(2 *  i      - 1, KIND = 8) / 24.0d0 * pi)
END DO
window_start(19:24) = 1.0d0
window_start(31:36) = 0.0d0
window_stop ( 1: 6) = 0.0d0
window_stop (13:18) = 1.0d0
RETURN
END SUBROUTINE mdct_initialize
!-------------------------------------------------------------------------------------------------
SUBROUTINE cp_subband(subband, subbuff)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: subband(:, :, :) ! 32 * 36 * nchannel
REAL (KIND = 8), INTENT(OUT) :: subbuff(:, :, :) ! 32 * 54 * nchannel
subbuff = EOSHIFT(subbuff, 36, 0.0d0, 2)
subbuff(:     , 19:54  , :) =  subband(:, :, :)           
subbuff(2:32:2, 20:54:2, :) = -subbuff(2:32:2, 20:54:2, :) ! ISO Figure A.4 Layer III decoder diagram
RETURN                                                     !    ~~~~~~~~~~~~ not written in the text but only in Figure 
END SUBROUTINE cp_subband
!--------------------------------------------------------------------------------------------------
SUBROUTINE sub_mdct(subband, r_mdct, mblock_type, q_alias)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: subband(:, :, :)
REAL (KIND = 8), INTENT(OUT) :: r_mdct (:, :, :, :)
INTEGER        , INTENT(IN ) :: mblock_type(:)
LOGICAL        , INTENT(IN ) :: q_alias 
INTEGER :: igranule, ichannel
LOGICAL, SAVE :: qfirst = .TRUE.
IF (qfirst) THEN
 qfirst = .FALSE.
 ALLOCATE( subbuff(32, 54, SIZE(subband, 3)) )
 subbuff = 0.0d0
 CALL mdct_initialize()
END IF
CALL cp_subband(subband, subbuff)
DO ichannel = 1, SIZE(subband, 3)
 DO igranule = 1, 2
  SELECT CASE (mblock_type(igranule))
   CASE ( 0) ! long block 
    CALL sub_mdct_normal(igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL anti_alias(r_mdct(1:32, :, igranule, ichannel))
   CASE (10) ! start block : long -> short
    CALL sub_mdct_start (igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL anti_alias(r_mdct(1:32, :, igranule, ichannel))
   CASE (11) ! start block : long -> mixed
    CALL sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL sub_mdct_start (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL anti_alias(r_mdct(1:32, :, igranule, ichannel))
   CASE (30) ! stop block : short -> long
    CALL sub_mdct_stop  (igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL anti_alias(r_mdct(1:32, :, igranule, ichannel))
   CASE (31) ! stop block : mixed -> long
    CALL sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL sub_mdct_stop  (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL anti_alias(r_mdct(1:32, :, igranule, ichannel))
   CASE (20) ! short block : short mode
    CALL sub_mdct_short(igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
   CASE (21) ! short block : mixed mode
    CALL sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    CALL sub_mdct_short (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
    IF (q_alias) CALL anti_alias(r_mdct(1:2, :, igranule, ichannel))
   CASE DEFAULT
    STOP 'error : SUBROUTINE sub_mdct : block_type error'
  END SELECT 
 END DO
END DO
RETURN
END SUBROUTINE sub_mdct
!--------------------------------------------------------------------------------------------------
SUBROUTINE sub_mdct_normal(igranule, n0, n1, subbuff, r_mdct)
IMPLICIT NONE
INTEGER        , INTENT(IN ) :: igranule, n0, n1 
REAL (KIND = 8), INTENT(IN ) :: subbuff(:, :)
REAL (KIND = 8), INTENT(OUT) :: r_mdct (:, :)
REAL (KIND = 8)              :: wk(36)
INTEGER :: iband, k
k = 18 * (igranule - 1) + 1
DO iband = n0, n1 
 wk = window_n * subbuff(iband, k:k + 35)
 CALL mdct36(wk, r_mdct(iband, :))
END DO
RETURN
END SUBROUTINE sub_mdct_normal
!--------------------------------------------------------------------------------------------------
SUBROUTINE sub_mdct_start(igranule, n0, n1, subbuff, r_mdct)
IMPLICIT NONE
INTEGER        , INTENT(IN ) :: igranule, n0, n1 
REAL (KIND = 8), INTENT(IN ) :: subbuff(:, :)
REAL (KIND = 8), INTENT(OUT) :: r_mdct (:, :)
REAL (KIND = 8)              :: wk(36)
INTEGER :: iband, k
k = 18 * (igranule - 1) + 1
DO iband = n0, n1 
 wk = window_start * subbuff(iband, k:k + 35)
 CALL mdct36(wk, r_mdct(iband, :))
END DO
RETURN
END SUBROUTINE sub_mdct_start
!--------------------------------------------------------------------------------------------------
SUBROUTINE sub_mdct_stop(igranule, n0, n1, subbuff, r_mdct)
IMPLICIT NONE
INTEGER        , INTENT(IN ) :: igranule, n0, n1 
REAL (KIND = 8), INTENT(IN ) :: subbuff(:, :)
REAL (KIND = 8), INTENT(OUT) :: r_mdct (:, :)
REAL (KIND = 8)              :: wk(36)
INTEGER :: iband, k
k = 18 * (igranule - 1) + 1
DO iband = n0, n1 
 wk = window_stop * subbuff(iband, k:k + 35)
 CALL mdct36(wk, r_mdct(iband, :))
END DO
RETURN
END SUBROUTINE sub_mdct_stop
!--------------------------------------------------------------------------------------------------
SUBROUTINE sub_mdct_short(igranule, n0, n1, subbuff, r_mdct)
IMPLICIT NONE
INTEGER        , INTENT(IN ) :: igranule, n0, n1
REAL (KIND = 8), INTENT(IN ) :: subbuff(:, :)
REAL (KIND = 8), INTENT(OUT) :: r_mdct (:, :)
REAL (KIND = 8)              :: wk(12)
INTEGER :: iband, k, m, m0, n 
m0 = 18 * (igranule - 1) + 1 + 6
DO k = 1, 3
 m = m0 + 6 * (k - 1) 
 n =  1 + 6 * (k - 1)
 DO iband = n0, n1 
  wk = window_s * subbuff(iband, m:m + 11)
  CALL mdct12(wk, r_mdct(iband, n:n + 5))
 END DO
END DO
RETURN
END SUBROUTINE sub_mdct_short
!--------------------------------------------------------------------------------------------------
SUBROUTINE anti_alias(x) ! C.1.5.3.3 Aliasing-Butterfly, 2.4.3.4.10.1 Alias Reduction, Table B.9  
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN OUT) :: x(:, :) ! x(32, 18)
REAL (KIND = 8) :: tmp1, tmp2
INTEGER :: iband, k
DO iband = 1, SIZE(x, 1) - 1
 DO k = 0, 7
  tmp1 = x(iband    , 18 - k) 
  tmp2 = x(iband + 1,  1 + k) 
  x(iband    , 18 - k) =  tmp1 * cs(k) + tmp2 * ca(k)
  x(iband + 1,  1 + k) =  tmp2 * cs(k) - tmp1 * ca(k) 
 END DO
END DO
RETURN
END SUBROUTINE anti_alias
!--------------------------------------------------------------------------------------------------
SUBROUTINE mdct36(x_in, x_out) ! ISO C.1.5.3.3 MDCT:   ! MDCT(4n) <=> FFT(n) proof: http://members.at.infoseek.co.jp/kitaurawa/mdct.pdf   (JP)
IMPLICIT NONE                                          !                            http://members.at.infoseek.co.jp/kitaurawa/mdct_e.pdf (E)
REAL (KIND = 8), INTENT(IN ) :: x_in (:) ! 36  
REAL (KIND = 8), INTENT(OUT) :: x_out(:) ! 18
REAL (KIND = 8):: x_shift(36), x_re, x_im
COMPLEX (KIND = 8):: fft(9)
INTEGER:: i
x_shift( 1: 9) = -x_in(28:36)
x_shift(10:36) =  x_in( 1:27)
DO i = 1, 9
 x_re = x_shift(2 * i -  1) - x_shift(38 - 2 * i) 
 x_im = x_shift(2 * i + 17) - x_shift(20 - 2 * i) 
 fft(i) = CMPLX(x_re, x_im, KIND = 8) * omega9s(i)
END DO
CALL fft9(fft)
fft = fft * omega9s
DO i = 1, 9
 x_out(2 * i - 1) =  REAL(      fft(i)      , KIND = 8) 
 x_out(2 * i    ) =  REAL(AIMAG(fft(10 - i)), KIND = 8)  
END DO
RETURN
END SUBROUTINE mdct36
!---------------------------------------------------------------------------------------------------
SUBROUTINE fft9(fft) ! FFT of 3^n (n=2)
IMPLICIT NONE
COMPLEX (KIND = 8), INTENT(IN OUT):: fft(:)
COMPLEX (KIND = 8) :: c1, c2, c3, c4, tmp1, tmp2, tmp3
INTEGER, PARAMETER :: nn = 9
INTEGER :: i, j, iphase1, iphase2, m1, m2, m3, k3, kn3
fft = fft(indx9) / 9.0d0 ! reorder and normalize
DO k3 = 1, 2 ! 3^n (n=2)
 kn3 = 3**(k3 - 1)
 DO i = 1, nn, 3 * kn3 
  DO j = 1, kn3
   iphase1 = 3**(2 - k3) * (j - 1) 
   iphase2 = 2 * iphase1
   c1 = omega9( MOD(iphase1, nn) )
   c2 = omega9( MOD(iphase2, nn) )
   m1 = i + j - 1
   m2 = m1 + kn3
   m3 = m2 + kn3
   tmp1 =      fft(m1)
   tmp2 = c1 * fft(m2)
   tmp3 = c2 * fft(m3)
   fft(m1) = tmp1 + tmp2 + tmp3
   c3 = tmp1 - 0.5d0 * ( tmp2 + tmp3 )
   c4 =     sqrt3_2  * ( tmp2 - tmp3 ) ! sqrt3_2 = i sqrt(3) / 2
   fft(m2) = c3 + c4
   fft(m3) = c3 - c4
  END DO
 END DO
END DO
RETURN
END SUBROUTINE fft9
!-----------------------------------------------------------------------
SUBROUTINE fft9_initialize()
IMPLICIT NONE
INTEGER :: i, j, k, n
pi = 4.0d0 * ATAN(1.0d0)
sqrt3_2 = CMPLX(0.0d0, SQRT(0.75d0), KIND = 8) ! iSQRT(3) / 2
DO i = 1, 9
 n = 0
 k = i - 1
 DO j = 2 - 1, 0, -1
  n = n + MOD(k, 3) * 3**j
  k = k / 3
 END DO
 indx9(i) = n + 1
 !omega9(i - 1) = EXP( CMPLX(0.0d0, 2.0d0 * pi /  9.0d0 * REAL(i - 1, KIND = 8), KIND = 8) )
 !omega9s(i)    = EXP( CMPLX(0.0d0, 2.0d0 * pi / 36.0d0 * ( 1.0d0 / 8.0d0 + REAL(i - 1, KIND = 8) ), KIND = 8) )
 omega9(i - 1) = CMPLX( COS( REAL( 2 * i - 2, KIND = 8) /   9.0d0 * pi ),   &
                        SIN( REAL( 2 * i - 2, KIND = 8) /   9.0d0 * pi ), KIND = 8 ) 
 omega9s(i)    = CMPLX( COS( REAL( 8 * i - 7, KIND = 8) / 144.0d0 * pi ), & 
                        SIN( REAL( 8 * i - 7, KIND = 8) / 144.0d0 * pi ), KIND = 8 ) 
END DO
RETURN
END SUBROUTINE fft9_initialize
!-----------------------------------------------------------------------
SUBROUTINE mdct12(x_in, x_out)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: x_in (:) ! 12
REAL (KIND = 8), INTENT(OUT) :: x_out(:) !  6
REAL (KIND = 8):: x_shift(12), x_re, x_im  
COMPLEX (KIND = 8):: fft(3)
INTEGER:: i  
x_shift( 1: 3) = -x_in(10:12)
x_shift( 4:12) =  x_in( 1: 9)
DO i = 1, 3
 x_re = x_shift(2 * i - 1) - x_shift(14 - 2 * i) 
 x_im = x_shift(2 * i + 5) - x_shift( 8 - 2 * i) 
 fft(i) = CMPLX(x_re, x_im, KIND = 8) * omega3s(i)
END DO
CALL fft3(fft)
fft = fft * omega3s
DO i = 1, 3
 x_out(2 * i - 1) =  REAL(      fft(i)     , KIND = 8) 
 x_out(2 * i    ) =  REAL(AIMAG(fft(4 - i)), KIND = 8)  
END DO
RETURN
END SUBROUTINE mdct12
!-----------------------------------------------------------------------
SUBROUTINE fft3(fft)
IMPLICIT NONE
COMPLEX (KIND = 8), INTENT(IN OUT):: fft(:)
COMPLEX (KIND = 8) :: c2, c3, c4
c2 = fft(2) + fft(3)
c3 = fft(1) - 0.5d0 * c2
c4 = sqrt3_2 * ( fft(2) - fft(3) ) ! sqrt3_2 = i sqrt(3) / 2
fft(1) = fft(1) + c2
fft(2) = c3 + c4
fft(3) = c3 - c4
fft = fft / 3.0d0 ! normalization
RETURN
END SUBROUTINE fft3
!--------------------------------------------------------------------------------------------------
SUBROUTINE fft3_initialize()
IMPLICIT NONE
INTEGER :: i
pi = 4.0d0 * ATAN(1.0d0)
sqrt3_2 = CMPLX(0.0d0, SQRT(0.75d0), KIND = 8) ! iSQRT(3) / 2
DO i = 1, 3
 !omega3(i - 1) = EXP( CMPLX(0.0d0, 2.0d0 * pi /  3.0d0 * REAL(i - 1, KIND = 8), KIND = 8) )
 !omega3s(i)    = EXP( CMPLX(0.0d0, 2.0d0 * pi / 12.0d0 * ( 1.0d0 / 8.0d0 + REAL(i - 1, KIND = 8) ), KIND = 8) )
 omega3(i - 1) = CMPLX( COS( REAL( 2 * i - 2, KIND = 8) /  3.0d0 * pi ),   &
                        SIN( REAL( 2 * i - 2, KIND = 8) /  3.0d0 * pi ), KIND = 8 ) 
 omega3s(i)    = CMPLX( COS( REAL( 8 * i - 7, KIND = 8) / 48.0d0 * pi ),   & 
                        SIN( REAL( 8 * i - 7, KIND = 8) / 48.0d0 * pi ), KIND = 8 ) 
END DO
RETURN
END SUBROUTINE fft3_initialize
!-----------------------------------------------------------------------
END MODULE mod_mdct
