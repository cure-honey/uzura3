MODULE mod_fft
IMPLICIT NONE
PRIVATE
PUBLIC :: init_fft, fft23, fft23han ! SUBROUTINE
INTEGER, PUBLIC :: indx576(1152), indx192(384)
COMPLEX (KIND = 8), PUBLIC :: omega576(0:1151), omega192(0:383), sqrt3_2
REAL     (KIND = 8), PUBLIC :: han576(1152), han192(384)
REAL     (KIND = 8) :: pi
CONTAINS
!-----------------------------------------------------------------------------------------------------------------
SUBROUTINE init_fft
IMPLICIT NONE
INTEGER :: i, k, n, m
pi = 4.0d0 * ATAN(1.0d0)
sqrt3_2 = CMPLX(0.0d0, SQRT(0.75d0), KIND = 8) ! iSQRT(3) / 2
!
indx576 = 1
DO i = 1, 1152
 DO k = 1, 7
  n = 2**(k - 1)
  m = 1152 / 2**k
  indx576(i) = indx576(i) + ( MOD(i - 1, 2 * m) / m ) * n
 END DO
 DO k = 1, 2
  n = 2**6 * 3**(k - 1) 
  m = 1152 / 2**6 / 3**k
  indx576(i) = indx576(i) + ( MOD(i - 1, 3 * m) / m ) * n
 END DO
 omega576(i - 1) = EXP( CMPLX(0.0d0, 2.0d0 * pi / 1152.0d0 * REAL(i - 1, KIND = 8), KIND = 8) )
! han576(i) = SQRT(8.0d0 / 3.0d0) * 0.5d0 *( 1.0d0 - COS(2.0d0 * pi * REAL(i - 1, KIND = 8) / 1152.0d0 ) )
 han576(i) = 0.5d0 *( 1.0d0 - COS(2.0d0 * pi * REAL(i - 1, KIND = 8) / 1152.0d0 ) )
END DO
!
indx192 = 1
DO i = 1, 192
 DO k = 1, 6
  n = 2**(k - 1)
  m = 384 / 2**k
  indx192(i) = indx192(i) + ( MOD(i - 1, 2 * m) / m ) * n
 END DO
 DO k = 1, 1
  n = 2**6 * 3**(k - 1) 
  m = 384 / 2**6 / 3**k
  indx192(i) = indx192(i) + ( MOD(i - 1, 3 * m) / m ) * n
 END DO
 omega192(i - 1) = EXP( CMPLX(0.0d0, 2.0d0 * pi / 384.0d0 * REAL(i - 1, KIND = 8), KIND = 8) )
! han192(i) = SQRT(8.0d0 / 3.0d0) * 0.5d0 *( 1.0d0 - COS(2.0d0 * pi * REAL(i - 1, KIND = 8) / 384.0d0 ) )
 han192(i) = 0.5d0 *( 1.0d0 - COS(2.0d0 * pi * REAL(i - 1, KIND = 8) / 384.0d0 ) )
END DO
END SUBROUTINE init_fft
!-----------------------------------------------------------------------------------------------------------------
SUBROUTINE fft23(np2, np3, indx, omega, fft)
INTEGER           , INTENT(IN    ) :: np2, np3, indx(:)
COMPLEX (KIND = 8), INTENT(IN    ) :: omega(0:)
COMPLEX (KIND = 8), INTENT(IN OUT) :: fft(:)
COMPLEX (KIND = 8) :: c1, c2, c3, c4, tmp1, tmp2, tmp3
INTEGER :: i, j, nn, iphase1, iphase2, m1, m2, m3, k3, kn3, k2, kn2
nn = 2**np2 * 3**np3
fft = fft(indx) / REAL(nn, KIND = 8) ! reorder and normalize
! 3**np3
DO k3 = 1, np3 ! 3^n (n=2)
 kn3 = 3**(k3 - 1)
 DO i = 1, nn, 3 * kn3 
  DO j = 1, kn3
   iphase1 = 2**np2 * 3**(np3 - k3) * (j - 1) 
   iphase2 = 2 * iphase1
   c1 = omega( MOD(iphase1, nn) )
   c2 = omega( MOD(iphase2, nn) )
   m1 = i + j - 1
   m2 = m1 + kn3
   m3 = m2 + kn3
   tmp1 =      fft(m1)
   tmp2 = c1 * fft(m2)
   tmp3 = c2 * fft(m3)
   fft(m1) = tmp1 + tmp2 + tmp3
   c3 = tmp1 - 0.5d0 * ( tmp2 + tmp3 )
   c4 =      sqrt3_2 * ( tmp2 - tmp3 ) ! sqrt3_2 = i sqrt(3) / 2
   fft(m2) = c3 + c4
   fft(m3) = c3 - c4
  END DO
 END DO
END DO
! 2**np2
DO k2 = 1, np2
 kn2 = 2**(k2 - 1) * 3**np3
 DO i = 1, nn, 2 * kn2
  DO j = 1, kn2  
   iphase1 = 2**(np2 - k2) * (j - 1) 
   c1 = omega( MOD(iphase1, nn) )
   m1 = i + j - 1
   m2 = m1 + kn2
   tmp2 = c1 * fft(m2)
   fft(m2) = fft(m1) - tmp2
   fft(m1) = fft(m1) + tmp2
  END DO
 END DO
END DO
RETURN
END SUBROUTINE fft23
!-----------------------------------------------------------------------------------------------------------------
SUBROUTINE fft23han(np2, np3, indx, omega, fft, han)
INTEGER           , INTENT(IN    ) :: np2, np3, indx(:)
COMPLEX (KIND = 8), INTENT(IN    ) :: omega(0:)
COMPLEX (KIND = 8), INTENT(IN OUT) :: fft(:)
REAL     (KIND = 8), INTENT(IN    ) :: han(:)
fft = fft * han
CALL fft23(np2, np3, indx, omega, fft)
RETURN
END SUBROUTINE fft23han
!-----------------------------------------------------------------------------------------------------------------
END MODULE mod_fft
!==================================================================================================================
MODULE mod_psycho
USE mod_mpg
USE mod_fft
IMPLICIT NONE
PRIVATE
PUBLIC :: psycho, calc_mask  ! subroutine
REAL (KIND = 8) :: pi, pi2
REAL (KIND = 8) :: ath_l(576) , ath_s(192, 3)
REAL (KIND = 8) :: sf_l(576, 576), sf_s(192, 192) 
REAL (KIND = 8) :: afft_l(576, 2, 2), afft_s(192, 3, 2, 2)
REAL (KIND = 8) :: phi_l(576, 2, 2), phi_s(192, 3, 2, 2)
REAL (KIND = 8) :: freq_l(576), freq_s(192), bark_l(576), bark_s(192), bw_l(576), bw_s(192)
REAL (KIND = 8) :: weight_l(576), weight_s(192)
INTEGER        :: ibark_l(576), ibark_s(192), ifb_l(25, 0:2), ifb_s(25, 0:2)
!
CONTAINS
!----------------------------------------------------------------
SUBROUTINE init_absolute_threshold(isample_rate)
IMPLICIT NONE
INTEGER, INTENT(IN ) :: isample_rate
REAL(KIND = 8):: freq, temp !, ath(576)
INTEGER :: i, k
pi  = 4.0d0 * ATAN(1.0d0)
pi2 = 2.0d0 * pi 
DO i = 1, 576
 freq = REAL(isample_rate, KIND = 8) / 2.0d0 / 1000.0d0 * (REAL(i, KIND = 8) - 0.0d0) / 576.0d0
!  temp =  3.64d0  * freq ** (-0.8d0) & 
!       -  6.50d0  * EXP(-0.6d0 * (freq -  3.3d0)**2.0d0) &
!       +  0.001d0 * freq ** 4.0d0 &      
!       + ath_min
 temp = 3.64d0 * freq ** (-0.8d0) &                     ! alternative ATH function 
      - 6.50d0 * EXP(-0.6d0 * (freq - 3.3d0)**2.0d0) &  ! reference: LAME ath-type 3  
      + 5.1d0 
!
 IF      (freq >  5.0d0 .AND. freq <=  5.5d0) THEN
  temp =  5.1d0 + 4.0d0 * (freq -  5.0d0) 
 ELSE IF (freq >  5.5d0 .AND. freq <=  8.0d0) THEN 
  temp =  7.1d0 + 2.76d0 * (freq -  5.5d0) 
 ELSE IF (freq >  8.0d0 .AND. freq <= 10.0d0) THEN
  temp = 14.0d0 + 1.00d0 * (freq -  8.0d0) 
 ELSE IF (freq > 10.0d0 .AND. freq <= 11.5d0) THEN
  temp = 16.0d0 + 0.33d0 * (freq - 10.0d0) 
 ELSE IF (freq > 11.5d0 .AND. freq <= 12.0d0) THEN
  temp = 16.5d0 + 2.0d0 * (freq - 11.5d0) 
 ELSE IF (freq > 12.0d0) THEN
  temp = 0.001d0 * freq ** 3.81d0 + 4.60d0  ! 12.93     
!
!  temp = 0.001d0 * freq ** 3.80d0 + 4.92d0  ! 12.61    
!  temp = 0.001d0 * freq ** 3.81d0 + 4.60d0  ! 12.93     
!  temp = 0.001d0 * freq ** 3.82d0 + 4.27d0  ! 13.26     
!  temp = 0.001d0 * freq ** 3.83d0 + 3.94d0  ! 13.59     
!  temp = 0.001d0 * freq ** 3.84d0 + 3.61d0  ! 13.93     
!  temp = 0.001d0 * freq ** 3.85d0 + 3.27d0  ! 14.28     
!  temp = 0.001d0 * freq ** 4.00d0 - 2.86d0  ! 20.70     
 END IF
 temp = MIN(temp + ath_min, ath_max)  
 ath_l(i) = 10.0d0**(temp / 20.0d0) 
END DO
!
DO i = 1, 192
 k = 3 * (i - 1) + 1
 ath_s(i, :) = MINVAL( ath_l(k:k + 2) ) 
END DO
RETURN
END SUBROUTINE init_absolute_threshold
!------------------------------------------------------------------------------------------------
FUNCTION switch_q(wx) RESULT(ires) ! attack detection (UZURA original)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN) :: wx(:, :)
INTEGER :: ires
REAL (KIND = 8), SAVE :: sum0a, sum1a, sum0b, sum1b
sum0a = sum1a
sum1a = SUM( wx(1:36, :) )
sum0b = sum1b
sum1b = SUM( wx(37: , :) ) 
IF  ( sum1a > switch * sum0a .OR. sum1b > switch * sum0b ) THEN
 ires = mblock_type_param
 IF (q_sm .AND. sum1a < xsm * sum0a .AND. sum0a < xsm * sum1a ) ires = 21 ! mixed 
ELSE 
 ires = 0 
END IF
IF (mblock_type_param == 0) ires = 0 ! force long-only mode
!debug info
if ( sum1a >= switch * sum0a ) nn1 = nn1 + 1
if ( sum1b >= switch * sum0b ) nn2 = nn2 + 1
RETURN
END FUNCTION switch_q
!--------------------------------------------------------------------------------------------------------------------
SUBROUTINE attack(wx, mblock_type) !..... ! ISO Figure C.7 (p.95) Window Switching State Diagram
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: wx(:, :, :)
INTEGER        , INTENT(OUT) :: mblock_type(:, :)
INTEGER :: i, iattack
INTEGER, SAVE :: mblock_prev = 0
DO i = 1, 2
 iattack = switch_q(wx(:, i, :))
 SELECT CASE (iattack)
  CASE (0) ! no-attack
    SELECT CASE (mblock_prev)
      CASE ( 0, 30, 31) ! long
        mblock_type(i, :) =  0
      CASE (10) 
        mblock_type(i, :) = 20
      CASE (11) 
        mblock_type(i, :) = 21
      CASE (20) ! short
        mblock_type(i, :) = 30
      CASE (21) ! mixed
        mblock_type(i, :) = 31
      CASE DEFAULT
        WRITE(*, *) 'error: psycho : unexpected block type: case 0-x', i, mblock_prev, mblock_type(i, 1) 
        STOP
      END SELECT
  CASE (20) ! attack-short
    SELECT CASE (mblock_prev)
      CASE ( 0, 30, 31)  
        mblock_type(i, :) = 10
      CASE (10) 
        mblock_type(i, :) = 20
      CASE (11) 
        mblock_type(i, :) = 21
      CASE (20) 
        mblock_type(i, :) = 20
      CASE (21) 
        mblock_type(i, :) = 21 
      CASE DEFAULT
        WRITE(*, *) 'error: psycho : unexpected block type: case 20-x', i, mblock_prev, mblock_type(i, 1) 
        STOP
      END SELECT
  CASE (21) ! attack-mixed
    SELECT CASE (mblock_prev)
      CASE ( 0, 30, 31)  
        mblock_type(i, :) = 11
      CASE (10) 
        mblock_type(i, :) = 20
      CASE (11) 
        mblock_type(i, :) = 21
      CASE (20) 
        mblock_type(i, :) = 20 
      CASE (21) 
        mblock_type(i, :) = 21
      CASE DEFAULT
        WRITE(*, *) 'error: psycho : unexpected block type: case 21-x', i, mblock_prev, mblock_type(i, 1) 
        STOP
      END SELECT
  CASE DEFAULT
    WRITE(*, *) 'error: psycho : unexpected block type'
    STOP
 END SELECT
 !mblock_type = 20 ! for debug
 mblock_prev = mblock_type(i, 1)
!....... debug info ...............................................
 select case ( mblock_type(i, 1) )
  case ( 0) 
   long = long + 1
  case (20)
   nshort = nshort + 1
  case (21) 
   mix = mix + 1
  case (10, 11)
   m1 = m1 + 1
  case (30, 31)
   m3 = m3 + 1
 end select
!
END DO
RETURN
END SUBROUTINE attack
!-----------------------------------------------------------------------------------------------------------------------
! not working correctly ; MDCT based version is better  
SUBROUTINE mid_side(mpg, wx, mblock_type) ! ISO  C.2.4.3.4.9.2,  G.2 MS_Stereo and intensity stereo coding Layer III
IMPLICIT NONE                    
REAL (KIND = 8)       , INTENT(IN    ) :: wx(:, :, :)
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
INTEGER        , INTENT(IN) :: mblock_type(:, :)
INTEGER         :: igranule, nchannel, n0, n1
REAL (KIND = 8) :: tmp1, tmp2
INTEGER, SAVE :: mode_old = 0, mode_ext_old = 0
LOGICAL         :: qms
nchannel = SIZE(wx, 3)
qms = qms_stereo
SELECT CASE (mpg%isample_rate) ! threshold ~ 7kHz (empirical value)
 CASE (0) ! 44.1kHz
  n0 = 183
 CASE (1) ! 48.0kHz
  n0 = 168
 CASE (2) ! 32.0kHz
  n0 = 252
 CASE DEFAULT
  STOP ' sample_rate error : SUBROUTINE mid_side '
END SELECT
n1 = n0 + 1
tmp1 = 0.0d0
tmp2 = 0.0d0
DO igranule = 1, 2
! pop musics often have different L-R behavior in low and high frequency 
 tmp1 = tmp1 + SUM( ABS( wx(1:n0, igranule, 1) - wx(1:n0, igranule, 2) ) ) 
 tmp2 = tmp2 + SUM(    ( wx(1:n0, igranule, 1) + wx(1:n0, igranule, 2) ) ) 
END DO
IF ( tmp1 > xms * tmp2 ) THEN 
 qms = .FALSE.
 ns1 = ns1 + 1
 ns  = ns  + 1
END IF
n1 = n0 + 1
tmp1 = 0.0d0
tmp2 = 0.0d0
DO igranule = 1, 2
! pop musics often have different L-R behavior in low and high frequency 
 tmp1 = tmp1 + SUM( ABS( wx(n1:, igranule, 1) - wx(n1:, igranule, 2) ) ) 
 tmp2 = tmp2 + SUM(    ( wx(n1:, igranule, 1) + wx(n1:, igranule, 2) ) ) 
END DO
IF ( ABS(tmp1) > xms * ABS(tmp2) ) THEN 
 qms = .FALSE.
 ns1 = ns1 + 1
 ns  = ns  + 1
END IF
!
IF (mblock_type(1, 1) /= 0 .AND. mblock_type(1, 1) /= 20 .AND. mblock_type(1, 1) /= 21) THEN  
 mpg%mode = 0 !mode_old 
 mpg%mode_extension = 0 ! mode_ext_old
ELSE
 IF (qms) THEN 
  mpg%mode            =  1 ! joint stereo
  mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
  ms = ms + 1 
 ELSE
  mpg%mode            =  0 ! normal stereo
  mpg%mode_extension  =  0 ! intensity_stereo off / ms_stereo off 
 END IF
END IF 
mode_old     = mpg%mode 
mode_ext_old = mpg%mode_extension 
RETURN
END SUBROUTINE mid_side
!------------------------------------------------------------------------------------------------
SUBROUTINE fft_long(nchannel, pcm, afft_l, phi_l)
IMPLICIT NONE
INTEGER       , INTENT(IN ) :: nchannel
REAL (KIND = 8), INTENT(IN ) :: pcm(:, :)
REAL (KIND = 8), INTENT(OUT) :: afft_l(:, :, :), phi_l(:, :, :)
COMPLEX (KIND = 8) :: fft576(1152, 2, 2)
INTEGER :: ichannel, igranule, m1,m2
DO igranule = 1, 2
 DO ichannel = 1, nchannel
  m1 = 1 + 480 * (igranule - 1)
  m2 = m1 + 1152 - 1
  fft576(:, igranule, ichannel) = CMPLX(PCM(m1:m2, ichannel), 0.0d0, KIND = 8) ! put 1152 real data -> get 576 complex FFT
!  CALL fft23han(7, 2, indx576, omega576, fft576(:, igranule, ichannel), han576 ) ! 2^7 * 3^2 = 1152
  CALL fft23(7, 2, indx576, omega576, fft576(:, igranule, ichannel) )
  afft_l (:, igranule, ichannel) = ABS(fft576(1:576, igranule, ichannel))
  afft_l (:, igranule, ichannel) = afft_l (:, igranule, ichannel) 
  phi_l(:, igranule, ichannel) = ATAN2(AIMAG(fft576(1:576, igranule, ichannel)), REAL(fft576(1:576, igranule, ichannel)))
 END DO
END DO
phi_l = phi_l + pi ! ATAN -pi~pi -> 0~2pi
RETURN
END SUBROUTINE fft_long
!------------------------------------------------------------------------------------------------
SUBROUTINE fft_short(nchannel, pcm, afft_s, phi_s)
IMPLICIT NONE
INTEGER       , INTENT(IN ) :: nchannel
REAL (KIND = 8), INTENT(IN ) :: pcm(:, :)
REAL (KIND = 8), INTENT(OUT) :: afft_s(:, :, :, :), phi_s(:, :, :, :)
COMPLEX (KIND = 8) :: fft192(384, 3, 2, 2)
INTEGER :: ichannel, igranule, iwin, m1, m2
DO igranule = 1, 2
 DO ichannel = 1, nchannel
  DO iwin = 1, 3
   m1 = 1 + 480 * (igranule - 1) + 384 * (iwin - 1)
   m2 = m1 + 384 - 1
   fft192(:, iwin, igranule, ichannel) = CMPLX(PCM(m1:m2, ichannel), 0.0d0, KIND = 8) ! put 384 real data -> get 192 complex FFT
 !  CALL fft23han(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel), han192 ) ! 2^7 * 3^1 = 384
   CALL fft23(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel) )
   afft_s (:, iwin, igranule, ichannel) = ABS(fft192(1:192, iwin, igranule, ichannel))
   afft_s (:, iwin, igranule, ichannel) = afft_s (:, iwin, igranule, ichannel) 
   phi_s(:, iwin, igranule, ichannel) = ATAN2(AIMAG(fft192(1:192, iwin, igranule, ichannel)), REAL(fft192(1:192, iwin, igranule, ichannel)))
  END DO
 END DO
END DO
phi_s = phi_s + pi ! ATAN2 -pi~pi -> 0~2pi
RETURN
END SUBROUTINE fft_short
!------------------------------------------------------------------------------------------------
SUBROUTINE calc_wx(nchannel, wx)
IMPLICIT NONE
INTEGER      , INTENT(IN ) :: nchannel
REAL(KIND = 8), INTENT(OUT) :: wx(:, :, :)
INTEGER :: igranule, ichannel
DO igranule = 1, 2
 DO ichannel = 1, nchannel
  wx(:, igranule, ichannel) = ( afft_l(:, igranule, ichannel) * weight_l )**2.0d0 
 END DO
END DO
RETURN
END SUBROUTINE calc_wx
!------------------------------------------------------------------------------------------------
SUBROUTINE psycho(pcm, mpg, mblock_type)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
INTEGER               , INTENT(   OUT) :: mblock_type(:, :)
REAL (KIND = 8)       , INTENT(IN    ) :: pcm(:, :)
LOGICAL, SAVE :: qfirst = .true.
INTEGER        :: igranule, ichannel, nchannel
real (kind = 8) :: wx(576, 2, 2), pm, pm0(2, 2)
!..... initialization .........................................................................................
IF (qfirst) THEN
 qfirst = .false.               
 !!!- to avoid VBR bug of portable MP3 player DIAMOND MULTIMEDIA RIO500. first frame bitrate must be less than average bitrate
 IF (q_vbr .and. q_rio500) mpg%ibit_rate = 8 ! force first frame 112kbps for RIO500 (firmware 1.15)  
 CALL init_absolute_threshold( mpeg_sample_rates(mpg%isample_rate) )
 CALL init_mask( mpeg_sample_rates(mpg%isample_rate) )
 CALL init_fft()
END IF
!..... FFT 576/192 ............................................................................................
nchannel = SIZE(mblock_type, 2) 
CALL fft_long (nchannel, pcm, afft_l, phi_l)
CALL fft_short(nchannel, pcm, afft_s, phi_s)
!..... weighted intensity .....................................................................................
CALL calc_wx(nchannel, wx)
!..... attack detection .......................................................................................
CALL attack(wx, mblock_type)
!..... MS/NS selection ........................................................................................
!CALL mid_side(mpg, wx, mblock_type) ! not working correctly : MDCT based model works better
!..... VBR  ...................................................................................................
IF (q_vbr) THEN ! simple implimentation 
 DO igranule = 1, 2
  DO ichannel = 1, nchannel
   pm0(igranule, ichannel) = SUM(wx(:, igranule, ichannel) * bark_l) !Psychoacoustic Moment (UZURA original) 
  END DO 
 END DO
 pm = SUM(pm0) / ( SUM(wx) + 1.0d-9 )
 IF (pm < 0.1d0) THEN
   mpg%ibit_rate =  1 
 ELSE IF (pm             < 1.0d0) THEN
   mpg%ibit_rate =  9  
 ELSE IF (pm * pm_factor < 2.5d0) THEN
   mpg%ibit_rate = 10  
 ELSE IF (pm * pm_factor < 3.5d0) THEN
   mpg%ibit_rate = 11  
 ELSE IF (pm * pm_factor < 5.0d0) THEN
   mpg%ibit_rate = 12 
 ELSE IF (pm * pm_factor < 7.0d0) THEN 
   mpg%ibit_rate = 13 
 ELSE 
   mpg%ibit_rate = 14
 END IF
END IF
nbits(mpg%ibit_rate) = nbits(mpg%ibit_rate) + 1
RETURN
END SUBROUTINE psycho
!------------------------------------------------------------------------------------------------
SUBROUTINE init_mask(nsample_rate)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nsample_rate
INTEGER :: i, j
REAL (KIND  = 8) :: f0, f1
DO i = 1, 576
 freq_l(i) = REAL(nsample_rate, KIND = 8) / 2.0d0 * (REAL(i, KIND = 8) - 0.5d0) / 576.0d0 ! kHz
 bark_l(i) = bark(freq_l(i) / 1000.0d0) 
 ibark_l(i) = INT(bark_l(i) + 0.1d0) + 1 
 f0 = REAL(nsample_rate, KIND = 8) / 2000.0d0 * REAL(i - 1, KIND = 8) / 576.0d0
 f1 = REAL(nsample_rate, KIND = 8) / 2000.0d0 * REAL(i    , KIND = 8) / 576.0d0
 bw_l(i) = bark(f1) - bark(f0)
 weight_l(i) = ( bark(f1) - bark(f0) ) / (f1 - f0) 
END DO
DO i = 1, 192
 freq_s(i) = REAL(nsample_rate, KIND = 8) / 2.0d0 * (REAL(i, KIND = 8) - 0.5d0) / 192.0d0 ! kHz
 bark_s(i) = bark(freq_s(i) / 1000.0d0)  
 ibark_s(i) = INT(bark_s(i)) + 1
 f0 = REAL(nsample_rate, KIND = 8) / 2000.0d0 * REAL(i - 1, KIND = 8) / 192.0d0
 f1 = REAL(nsample_rate, KIND = 8) / 2000.0d0 * REAL(i    , KIND = 8) / 192.0d0
 bw_s(i) = bark(f1) - bark(f0)
 weight_s(i) = ( bark(f1) - bark(f0) ) / (f1 - f0) 
END DO
!
ifb_l(1, 1) = 1
DO i = 1, 25
 DO j = 1, 576
  IF (ibark_l(j) == i - 1) ifb_l(i, 1) = j + 1 
  IF (ibark_l(j) == i    ) ifb_l(i, 2) = j  
 END DO
 ifb_l(i, 0) = ifb_l(i, 2) - ifb_l(i, 1) + 1
END DO
ifb_s(1, 1) = 1
DO i = 1, 25
 DO j = 1, 192
  IF (ibark_s(j) == i - 1) ifb_s(i, 1) = j + 1 
  IF (ibark_s(j) == i    ) ifb_s(i, 2) = j  
 END DO
 ifb_s(i, 0) = ifb_s(i, 2) - ifb_s(i, 1) + 1
END DO
!
DO i = 1, 576
 DO j = 1, 576
  sf_l(i, j) = 10.0d0 ** ( spreading_function( bark_l(i) - bark_l(j) ) / 20.0d0 ) 
 END DO
END DO 
DO i = 1, 192
 DO j = 1, 192
  sf_s(i, j) = 10.0d0 ** ( spreading_function( bark_s(i) - bark_s(j) ) / 20.0d0 ) 
 END DO
END DO 
RETURN
END SUBROUTINE init_mask
!------------------------------------------------------------------------------------------------
SUBROUTINE calc_mask(igranule, ichannel, mblock_type, xmask, xnoise)
IMPLICIT NONE
INTEGER       , INTENT(IN    ) :: igranule, ichannel, mblock_type
REAL (KIND = 8), INTENT(   OUT) :: xmask(:, :), xnoise(:, :)
REAL (KIND = 8) :: x0_l(576), x0_s(192, 3), y0_l(576), y0_s(192, 3)
REAL (KIND = 8) :: d2phi_l(576), tone_l(576), fk_l(576), fl_l(576), tn_l(576)
REAL (KIND = 8) :: d2phi_s(192), tone_s(192), fk_s(192), fl_s(192), tn_s(192)
REAL (KIND = 8) :: yn(25)
REAL (KIND = 8), SAVE :: x1_l(576) = 0.0d0, x1_s(192, 3) = 0.0d0
REAL (KIND = 8), SAVE :: phi1_l(576, 2) = 0.0d0, phi2_l(576, 2) = 0.0d0
REAL (KIND = 8), SAVE :: phi1_s(192, 2) = 0.0d0, phi2_s(192, 2) = 0.0d0
REAL (KIND = 8), SAVE :: p0_l(576, 2) = 0.0d0, p1_l(576, 2) = 0.0d0, p2_l(576, 2) = 0.0d0
REAL (KIND = 8), SAVE :: p0_s(192, 2) = 0.0d0, p1_s(192, 2) = 0.0d0, p2_s(192, 2) = 0.0d0
REAL (KIND = 8), SAVE :: p3_s(192, 2) = 0.0d0, p4_s(192, 2) = 0.0d0, p5_s(192, 2) = 0.0d0
INTEGER :: icritical_band, iwin
!---------------------------------------------------------------------------------------------
! masking / allowed noise : reference Bosse Lincoln "An Experimental High Fidelity Perceptual Audio Coder Project in MUS420 Win 97"
!---------------------------------------------------------------------------------------------
d2phi_l = phi_l(:, igranule, ichannel) + phi2_l(:, ichannel) - 2.0d0 * phi1_l(:, ichannel) 
p0_l(:, ichannel) = MOD( ABS(d2phi_l), pi2 ) / pi2  
tone_l = 1.0d0 - MAX(p0_l(:, ichannel), p1_l(:, ichannel), p2_l(:, ichannel) )  
! mask for long block
fk_l =  0.3d0 * tone_l +  0.5d0 * (1 - tone_l) 
fl_l = 34.0d0 * tone_l + 20.0d0 * (1 - tone_l) 
tn_l = 10.0d0**( - ( fk_l * bark_l + fl_l + offset ) / 20.0d0 )  
x0_l = MATMUL(sf_l, afft_l(:, igranule, ichannel) * tn_l * weight_l) + ath_l
x1_l = tempo * x1_l + (1.0d0 - tempo) * x0_l ! temporal masking
x0_l = MAX(x0_l, x1_l)
! allowed noise for long block : average mask over a critical band
DO icritical_band = 1, 25
 yn(icritical_band) = SUM( x0_l(ifb_l(icritical_band, 1):ifb_l(icritical_band, 2)) ) / REAL(ifb_l(icritical_band, 0), KIND = 8)
END DO
y0_l = MAX( ath_l, yn(ibark_l) ) 
! save old data
p1_l(:, ichannel) = p0_l(:, ichannel)
p2_l(:, ichannel) = p1_l(:, ichannel)
phi1_l(:, ichannel) = phi_l(:, igranule, ichannel)
phi2_l(:, ichannel) = phi1_l(:, ichannel)
!
DO iwin = 1, 3
 d2phi_s = phi_s(:, iwin, igranule, ichannel) + phi2_s(:, ichannel) - 2.0d0 * phi1_s(:, ichannel)
 p0_s(:, ichannel) = MOD( ABS(d2phi_s), pi2 ) / pi2
 tone_s = 1.0d0 - MAX( p0_s(:, ichannel), p1_s(:, ichannel), p2_s(:, ichannel), & 
                       p3_s(:, ichannel), p4_s(:, ichannel), p5_s(:, ichannel)  )
! mask for short block
 fk_s =  0.3d0 * tone_s +  0.5d0 * (1 - tone_s)  
 fl_s = 34.0d0 * tone_s + 20.0d0 * (1 - tone_s) 
 tn_s = 10.0d0**( - ( fk_s * bark_s + fl_s + offset ) / 20.0d0 )
 x0_s(:, iwin) = MATMUL(sf_s, afft_s(:, iwin, igranule, ichannel) * tn_s * weight_s) + ath_s(:, iwin)
! save old data
 p1_s(:, ichannel) = p0_s(:, ichannel)
 p2_s(:, ichannel) = p1_s(:, ichannel)
 p3_s(:, ichannel) = p2_s(:, ichannel)
 p4_s(:, ichannel) = p3_s(:, ichannel)
 p5_s(:, ichannel) = p4_s(:, ichannel)
 phi1_s(:, ichannel) = phi_s(:, iwin, igranule, ichannel)
 phi2_s(:, ichannel) = phi1_s(:, ichannel)
END DO
x1_s(:, 1) = tempo * x1_s(:, 3) + (1.0d0 - tempo) * x0_s(:, 1) ! temporal masking   
x1_s(:, 2) = tempo * x1_s(:, 1) + (1.0d0 - tempo) * x0_s(:, 2) ! 
x1_s(:, 3) = tempo * x1_s(:, 2) + (1.0d0 - tempo) * x0_s(:, 3) !  
x0_s  = MAX(x0_s, x1_s)
! allowed noise for short block  : average mask over a critical band
DO iwin = 1, 3
 DO icritical_band = 1, 25
  yn(icritical_band) = SUM( x0_s(ifb_s(icritical_band, 1):ifb_s(icritical_band, 2), iwin) ) / REAL(ifb_s(icritical_band, 0), KIND = 8) 
 END DO
 y0_s(:, iwin) = MAX( ath_s(:, iwin), yn(ibark_s(:)) ) 
END DO
! order to r_mdct style
SELECT CASE (mblock_type)
 CASE (0, 10, 11, 30, 31)
  CALL deorder_l(x0_l, xmask )
  CALL deorder_l(y0_l, xnoise)
 CASE (20)
  CALL deorder_s(x0_s, xmask )
  CALL deorder_s(y0_s, xnoise)
 CASE (21) 
  CALL deorder_m(x0_l, x0_s, xmask )
  CALL deorder_m(y0_l, y0_s , xnoise)
 CASE DEFAULT
  STOP 'SUBROUTINE: calc_noise'
END SELECT
RETURN
END SUBROUTINE calc_mask
!------------------------------------------------------------------------------------------------
SUBROUTINE deorder_l(zl, xth)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: zl(:)
REAL (KIND = 8), INTENT(OUT) :: xth(:, :)
INTEGER :: i, iband
DO iband = 1, 32
 DO i = 1, 18
  xth(iband, i) = zl(18 * (iband - 1) + i)
 END DO
END DO
RETURN
END SUBROUTINE deorder_l
!------------------------------------------------------------------------------------------------
SUBROUTINE deorder_s(zs, xth)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: zs(:, :)
REAL (KIND = 8), INTENT(OUT) :: xth(:, :)
INTEGER :: i, iwin, iband
DO iband = 1, 32
 DO iwin = 1, 3 
  DO i = 1, 6
   xth(iband, 6 * (iwin - 1) + i) = zs(6 * (iband - 1) + i, iwin) 
  END DO
 END DO  
END DO
RETURN
END SUBROUTINE deorder_s
!------------------------------------------------------------------------------------------------
SUBROUTINE deorder_m(zl, zs, xth)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: zl(:), zs(:, :)
REAL (KIND = 8), INTENT(OUT) :: xth(:, :)
INTEGER :: i, iwin, iband
DO iband = 1, 2
 DO i = 1, 18
  xth(iband, i) = zl( 18 * (iband - 1) + i)
 END DO
END DO
DO iband = 3, 32
 DO iwin = 1, 3 
  DO i = 1, 6
   xth(iband, 6 * (iwin - 1) + i) = zs(6 * (iband - 1) + i, iwin) 
  END DO
 END DO  
END DO
RETURN
END SUBROUTINE deorder_m
!------------------------------------------------------------------------------------------------
FUNCTION bark(f) RESULT(res)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN) :: f
REAL (KIND = 8) :: res
res = 13.0d0 * ATAN(0.76d0 * f) + 3.5d0 * ATAN( (f / 7.5d0)**2.0d0 )
RETURN
END FUNCTION bark
!------------------------------------------------------------------------------------------------
FUNCTION spreading_function(z) RESULT(res)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN) :: z
REAL (KIND = 8) :: res
res = 15.81d0 + 7.5d0 * (z + 0.474d0) - 17.5d0 * SQRT(1.0d0 + (z + 0.474d0)**2.0d0) 
RETURN
END FUNCTION spreading_function
!------------------------------------------------------------------------------------------------
FUNCTION spreading_function0(z) RESULT(res)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN) :: z
REAL (KIND = 8) :: res
IF (z < 0.0d0) THEN
 res =  25.0d0 * z
ELSE
 res = -10.0d0 * z
END IF 
RETURN
END FUNCTION spreading_function0
!------------------------------------------------------------------------------------------------
FUNCTION spreading_function2(z, y) RESULT(res)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN) :: z, y
REAL (KIND = 8) :: res
res = (15.81d0 - y) + 7.5d0 * (z + 0.474d0) - (17.5d0 - y) * SQRT(1.0d0 + (z + 0.474d0)**2.0d0) 
RETURN
END FUNCTION spreading_function2
!------------------------------------------------------------------------------------------------
END MODULE mod_psycho