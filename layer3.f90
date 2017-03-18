!==============================================================================================
MODULE mod_layer3 
USE mod_mpg
USE mod_psycho
USE mod_inner_loop
IMPLICIT NONE
PRIVATE
PUBLIC  :: alloc_bits, init_scalefactor_bands ! subroutine 
PUBLIC  :: side_info, scfct                   ! variable
PUBLIC  :: scale_factor                       ! type 
INTEGER, PARAMETER :: npretab(0:20, 0:1) = &  ! ISO Table B.6 Layer III preemphasis (pretab)
         RESHAPE((/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2/), SHAPE(npretab))
!                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
TYPE :: scale_factor
 INTEGER :: long(0:20)
 INTEGER :: ishort(0:12, 3)
END TYPE scale_factor
!
TYPE (scale_factor), SAVE :: scfct(2, 2)
TYPE (side_info_  ), SAVE :: side_info
INTEGER :: iscfband_l(0:20, 3, 0:2), iscfband_s(0:11, 3, 0:2)
CONTAINS
!------------------------------------------------------------------------------------------------
SUBROUTINE init_scalefactor_bands(n) ! ISO Table B.8 Layer III scalefactor bands
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: i, j
iscfband_l = 0  
iscfband_s = 0
iscfband_l(:, 1, 2) = (/4, 4, 4, 4, 4, 4, 6, 6, 8, 10, 12, 16, 20, 24, 30, 38, 46, 56, 68, 84, 102/) !32.0kHz
iscfband_l(:, 1, 0) = (/4, 4, 4, 4, 4, 4, 6, 6, 8,  8, 10, 12, 16, 20, 24, 28, 34, 42, 50, 54,  76/) !44.1kHz
iscfband_l(:, 1, 1) = (/4, 4, 4, 4, 4, 4, 6, 6, 6,  8, 10, 12, 16, 18, 22, 28, 34, 40, 46, 54,  54/) !48.0kHz
iscfband_s(:, 1, 2) = (/4, 4, 4, 4, 6, 8, 12, 16, 20, 26, 34, 42/) !32.0kHz
iscfband_s(:, 1, 0) = (/4, 4, 4, 4, 6, 8, 10, 12, 14, 18, 22, 30/) !44.1kHz
iscfband_s(:, 1, 1) = (/4, 4, 4, 4, 6, 6, 10, 12, 14, 16, 20, 26/) !48.0kHz
DO j = 0, 2
 DO i = 1, 20
  iscfband_l(i    , 2, j) = iscfband_l(i - 1, 2, j) + iscfband_l(i - 1, 1, j)
  iscfband_l(i - 1, 3, j) = iscfband_l(i    , 2, j) - 1
 END DO 
 iscfband_l(20, 3, j) = iscfband_l(20, 1, j) + iscfband_l(20, 2, j) - 1
 DO i = 1, 11
  iscfband_s(i    , 2, j) = iscfband_s(i - 1, 2, j) + iscfband_s(i - 1, 1, j)
  iscfband_s(i - 1, 3, j) = iscfband_s(i    , 2, j) - 1
 END DO 
 iscfband_s(11, 3, j) = iscfband_s(11, 1, j) + iscfband_s(11, 2, j) - 1
END DO
iscalefactorband_l = iscfband_l(:, :, n)
iscalefactorband_s = iscfband_s(:, :, n)
RETURN
END SUBROUTINE init_scalefactor_bands
!------------------------------------------------------------------------------------------------
SUBROUTINE alloc_bits(mblock_type, r_mdct, i_mdct, mpg, max_bits, ianc) ! ISO C.1.5.4 
IMPLICIT NONE
INTEGER               , INTENT(IN    ) :: mblock_type(:, :)
REAL (KIND = 8)       , INTENT(IN OUT) :: r_mdct(:, :, :, :)
INTEGER               , INTENT(   OUT) :: i_mdct(:, :, :), max_bits, ianc
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
INTEGER :: ibit, ibit2, ichannel, nchannel, igranule, iused_bits, nused_bits, & 
           mbits(SIZE(r_mdct, 3), SIZE(r_mdct, 4))
REAL (KIND = 8) ::  tot_int(SIZE(r_mdct, 3), SIZE(r_mdct, 4))
REAL (KIND = 8) ::  x_mask (SIZE(r_mdct, 1), SIZE(r_mdct, 2), SIZE(r_mdct, 3), SIZE(r_mdct, 4))
REAL (KIND = 8) ::  x_noise(SIZE(r_mdct, 1), SIZE(r_mdct, 2), SIZE(r_mdct, 3), SIZE(r_mdct, 4))
REAL (KIND = 8) ::  z_noise(SIZE(r_mdct, 1), SIZE(r_mdct, 2))
i_mdct = 0
nchannel = SIZE(r_mdct, 4)
CALL init_scale_factor()
CALL init_side_info(mblock_type)
CALL mid_side_mdct(mpg, r_mdct)
!CALL mid_side_fft(mpg, r_mdct)
CALL calc_totint(r_mdct, tot_int)
CALL get_maxbits(mpg, max_bits)
CALL mean_bits  (mpg, max_bits, tot_int, mbits, nused_bits)
ibit2 = 0 ! remain bits
DO igranule = 1, 2
 DO ichannel = 1, nchannel 
  ibit = mbits(igranule, ichannel) + ibit2
  CALL calc_mask(igranule, ichannel, mblock_type(igranule, ichannel), &
                                     x_mask(:, :, igranule, ichannel), x_noise(:, :, igranule, ichannel) )
  IF (q_mask) WHERE (ABS(r_mdct(:, :, igranule, ichannel)) < x_mask(:, :, igranule, ichannel)) &
                          r_mdct(:, :, igranule, ichannel) = 0.0d0 
! begin test
  if ( mpg%mode == 1 .and. ichannel == 1) then
   z_noise(:, :) = ABS( SIGN(1.0d0, r_mdct(:, :, igranule, 1)) * x_noise(:, :, igranule, 1) &
                      + SIGN(1.0d0, r_mdct(:, :, igranule, 2)) * x_noise(:, :, igranule, 2) ) / SQRT(2.0d0) ! Left
  else if ( mpg%mode == 1 .and. ichannel == 2) then
   z_noise(:, :) = ABS( SIGN(1.0d0, r_mdct(:, :, igranule, 1)) * x_noise(:, :, igranule, 1) &
                      - SIGN(1.0d0, r_mdct(:, :, igranule, 2)) * x_noise(:, :, igranule, 2) ) / SQRT(2.0d0) ! Right
  else
   z_noise(:, :) = x_noise(:, :, igranule, ichannel)
  end if
! end test
  IF ( tot_int(igranule, ichannel) /= 0.0d0 ) THEN
   CALL outer_loop(ibit, mblock_type(igranule, ichannel), &
                       r_mdct (:, :, igranule, ichannel), z_noise(:, :), &
                       i_mdct (:   , igranule, ichannel), side_info%sub(igranule, ichannel), &
                               scfct(igranule, ichannel), iused_bits)
  ELSE 
   iused_bits = 0
   i_mdct(:, igranule, ichannel) = 0
  END IF
  nused_bits = nused_bits + iused_bits  
  side_info%sub(igranule, ichannel)%ipart2_3_length = iused_bits  
  ibit2 = ibit - iused_bits
! debug info
  ntable (side_info%sub(igranule, ichannel)%itable_select(1)) = ntable (side_info%sub(igranule, ichannel)%itable_select(1)) + 1
  ntable (side_info%sub(igranule, ichannel)%itable_select(2)) = ntable (side_info%sub(igranule, ichannel)%itable_select(2)) + 1
  ntable (side_info%sub(igranule, ichannel)%itable_select(3)) = ntable (side_info%sub(igranule, ichannel)%itable_select(3)) + 1
  ntab_ab(side_info%sub(igranule, ichannel)%icount1table_select) &
                                                              = ntab_ab(side_info%sub(igranule, ichannel)%icount1table_select) + 1
 END DO
END DO
ianc = max_bits - nused_bits
RETURN
END SUBROUTINE alloc_bits
!------------------------------------------------------------------------------------------------
SUBROUTINE init_scale_factor()
IMPLICIT NONE
INTEGER :: igranule, ichannel
DO ichannel = 1, 2
 DO igranule = 1, 2
  scfct(igranule, ichannel)%long   = 0
  scfct(igranule, ichannel)%ishort = 0
 END DO
END DO
RETURN
END SUBROUTINE init_scale_factor
!------------------------------------------------------------------------------------------------
SUBROUTINE init_side_info(mblock_type)
IMPLICIT NONE
INTEGER, INTENT(IN) :: mblock_type(:, :)
INTEGER :: igranule, ichannel
side_info%main_data_begin = 0                                     !  9 bits
side_info%iprivate_bits   = 0                                     !  5 / 3 bits
side_info%iscfsi(4, 2)    = 0                                     !  1 bit 
DO igranule = 1, 2
 DO ichannel = 1, SIZE(mblock_type, 2)
  side_info%sub(igranule, ichannel)%ipart2_3_length        = 0    ! 12 bits
  side_info%sub(igranule, ichannel)%ibig_values            = 0    !  9 bits
  side_info%sub(igranule, ichannel)%iglobal_gain           = 0    !  8 bits
  side_info%sub(igranule, ichannel)%iscalefac_compress     = 0    !  4 bits   
  SELECT CASE (mblock_type(igranule, ichannel))
   CASE ( 0) ! long-block
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 0    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 0    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
   CASE (10) ! start-block for short
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 1    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
   CASE (11) ! start-block for mixed
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 1    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
   CASE (30) ! stop-block for short
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 3    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
   CASE (31) ! stop-block for short
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 3    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
   CASE (20) ! short-block
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 2    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
   CASE (21) ! mixed-block
    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
    side_info%sub(igranule, ichannel)%iblock_type            = 2    !  2 bits   
    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
   CASE DEFAULT
    STOP ' error : mblock_type '
  END SELECT
  side_info%sub(igranule, ichannel)%itable_select(:)       = 0    !  5 bits   
  side_info%sub(igranule, ichannel)%isubblock_gain(:)      = 0    !  3 bits   
  side_info%sub(igranule, ichannel)%iregion0_count         = 0    !  4 bits
  side_info%sub(igranule, ichannel)%iregion1_count         = 0    !  3 bits
  side_info%sub(igranule, ichannel)%ipreflag               = 0    !  1 bit    
  side_info%sub(igranule, ichannel)%iscalefac_scale        = 0    !  1 bit   
  side_info%sub(igranule, ichannel)%icount1table_select    = 0    !  1 bit
!.....
  side_info%sub(igranule, ichannel)%icount1                = 0    ! local use
 END DO
END DO
RETURN
END SUBROUTINE init_side_info
!------------------------------------------------------------------------------------------------
SUBROUTINE mid_side_mdct(mpg, r_mdct) ! ISO  C.2.4.3.4.9.2,  G.2 MS_Stereo and intensity stereo coding Layer III
IMPLICIT NONE                    
REAL (KIND = 8)       , INTENT(IN OUT) :: r_mdct(:, :, :, :)
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
INTEGER         :: igranule, nchannel, n0, n1
REAL (KIND = 8) :: tmp_ms(32, 18, 2, 2), tmp1, tmp2
LOGICAL         :: qms
nchannel = SIZE(r_mdct, 4)
qms = qms_stereo
SELECT CASE (mpg%isample_rate) ! threshold ~ 7kHz (empirical value)
 CASE (0) ! 44.1kHz
  n0 = 10
 CASE (1) ! 48.0kHz
  n0 =  9
 CASE (2) ! 32.0kHz
  n0 = 14
 CASE DEFAULT
  STOP ' sample_rate error : SUBROUTINE mid_side '
END SELECT
n1 = n0 + 1
DO igranule = 1, 2
! pop music often have different L-R behavior in low and high frequency 
 tmp1 = SUM( ABS( r_mdct( 1:n0, :, igranule, 1)**2 - r_mdct( 1:n0, :, igranule, 2)**2 ) ) 
 tmp2 = SUM(    ( r_mdct( 1:n0, :, igranule, 1)**2 + r_mdct( 1:n0, :, igranule, 2)**2 ) ) 
 IF ( tmp1 > xms * tmp2 ) THEN 
   qms = .FALSE.
   ns1 = ns1 + 1
   ns  = ns  + 1
 ELSE
   tmp1 = SUM( ABS( r_mdct(n1:32, :, igranule, 1)**2 - r_mdct(n1:32, :, igranule, 2)**2 ) ) 
   tmp2 = SUM(    ( r_mdct(n1:32, :, igranule, 1)**2 + r_mdct(n1:32, :, igranule, 2)**2 ) ) 
  IF ( tmp1 > xms * tmp2 * 0.95d0) THEN ! cheat
   qms = .FALSE.
   ns2 = ns2 + 1
   ns  = ns  + 1
  END IF 
 END IF
END DO

!
IF (qms) THEN 
 IF (nchannel == 2) THEN
  tmp_ms(:, :, :, 1) = ( r_mdct(:, :, :, 1) + r_mdct(:, :, :, 2) ) / SQRT(2.0d0)
  tmp_ms(:, :, :, 2) = ( r_mdct(:, :, :, 1) - r_mdct(:, :, :, 2) ) / SQRT(2.0d0)
 ELSE
  tmp_ms(:, :, :, 1) = ( r_mdct(:, :, :, 1) ) / SQRT(2.0d0)
  tmp_ms(:, :, :, 2) = ( r_mdct(:, :, :, 1) ) / SQRT(2.0d0)
 END IF
 r_mdct(:, :, :, 1) = tmp_ms(:, :, :, 1)
 r_mdct(:, :, :, 2) = tmp_ms(:, :, :, 2)
 mpg%mode            =  1 ! joint stereo
 mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
 ms = ms + 1
ELSE
 mpg%mode            =  0 ! normal stereo
 mpg%mode_extension  =  0 ! intensity_stereo off / ms_stereo off
END IF
RETURN
END SUBROUTINE mid_side_mdct
!------------------------------------------------------------------------------------------------
SUBROUTINE mid_side_fft(mpg, r_mdct) ! ISO  C.2.4.3.4.9.2,  G.2 MS_Stereo and intensity stereo coding Layer III
IMPLICIT NONE                    
REAL (KIND = 8)       , INTENT(IN OUT) :: r_mdct(:, :, :, :)
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
REAL (KIND = 8) :: tmp_ms(32, 18, 2, 2)
INTEGER :: nchannel
!
nchannel = SIZE(r_mdct, 4)
IF (mpg%mode ==  1) THEN 
 IF (nchannel == 2) THEN
  tmp_ms(:, :, :, 1) = ( r_mdct(:, :, :, 1) + r_mdct(:, :, :, 2) ) / SQRT(2.0d0)
  tmp_ms(:, :, :, 2) = ( r_mdct(:, :, :, 1) - r_mdct(:, :, :, 2) ) / SQRT(2.0d0)
 ELSE
  tmp_ms(:, :, :, 1) = ( r_mdct(:, :, :, 1) ) / SQRT(2.0d0)
  tmp_ms(:, :, :, 2) = ( r_mdct(:, :, :, 1) ) / SQRT(2.0d0)
 END IF
 r_mdct(:, :, :, 1) = tmp_ms(:, :, :, 1)
 r_mdct(:, :, :, 2) = tmp_ms(:, :, :, 2)
END IF
RETURN
END SUBROUTINE mid_side_fft
!------------------------------------------------------------------------------------------------
SUBROUTINE calc_totint(r_mdct, tot_int)
IMPLICIT NONE
REAL (KIND = 8) , INTENT(IN ) :: r_mdct (:, :, :, :)
REAL (KIND = 8) , INTENT(OUT) :: tot_int(:, :)
INTEGER :: igranule, ichannel
DO igranule = 1, 2
 DO ichannel = 1, SIZE(r_mdct, 4)
!  tot_int(igranule, ichannel) = SUM( r_mdct(:, :, igranule, ichannel) ** 2.0d0  )
  tot_int(igranule, ichannel) = SUM( ABS(r_mdct(:, :, igranule, ichannel)) )
 END DO
END DO
RETURN
END SUBROUTINE calc_totint
!------------------------------------------------------------------------------------------------
SUBROUTINE get_maxbits(mpg, max_bits) ! after ISO 2.4.3.1 & 2.4.2.3 padding (p.22)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN OUT) :: mpg
INTEGER               , INTENT(   OUT) :: max_bits
INTEGER                                :: idiff, islot_size
INTEGER, SAVE                          :: irest = 0
islot_size = 144000 * mpeg_bit_rates(mpg%ibit_rate, mpg%layer) / mpeg_sample_rates(mpg%isample_rate) 
idiff  = MOD(144000 * mpeg_bit_rates(mpg%ibit_rate, mpg%layer),  mpeg_sample_rates(mpg%isample_rate))
irest  = irest - idiff
IF (irest < 0) THEN     
 mpg%ipadding = 1
 irest = irest + mpeg_sample_rates(mpg%isample_rate)
ELSE
 mpg%ipadding = 0
END IF
max_bits = ( islot_size + mpg%ipadding ) * 8
RETURN
END SUBROUTINE get_maxbits 
!-------------------------------------------------------------------------------------------------
SUBROUTINE mean_bits(mpg, max_bits, tot_int, mbits, iused_bits)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN ) :: mpg
INTEGER               , INTENT(IN ) :: max_bits
REAL (KIND = 8)       , INTENT(IN ) :: tot_int(:, :)
INTEGER               , INTENT(OUT) :: mbits(:, :), iused_bits
INTEGER                             :: nbits, nchannel
nchannel = SIZE(mbits, 2)
IF (nchannel == 1) THEN
 iused_bits =  32 + 136 ! monoral
ELSE
 iused_bits =  32 + 256 ! stereo 
END IF
IF (mpg%icrc == 0) iused_bits = iused_bits + 16
nbits = max_bits - iused_bits
!
!!  distribute bits between 2 granules * n channels according to total intensity 
mbits = INT( (max_bits - iused_bits) * &
             (factor * tot_int / SUM(tot_int + EPSILON(0.0)) + (1.0d0 - factor) / 2 / nchannel ) )
!!             variavle part                                     consant part      
mbits(1, 1) = mbits(1, 1) + nbits - SUM(mbits)
RETURN
END SUBROUTINE mean_bits
!------------------------------------------------------------------------------------------------
SUBROUTINE outer_loop(iallowed_bits, iblock_type, r_mdct, x_noise, iwk, side, sc_fac, iused_bits) ! ISO C.1.5.4.3
IMPLICIT NONE
INTEGER             , INTENT(IN    ) :: iallowed_bits, iblock_type
INTEGER             , INTENT(   OUT) :: iused_bits
REAL (KIND = 8)      , INTENT(IN    ) :: r_mdct(:, :), x_noise(:, :)
INTEGER             , INTENT(   OUT) :: iwk(:)
TYPE (side_info_sub), INTENT(IN OUT) :: side
TYPE (scale_factor) , INTENT(IN OUT) :: sc_fac
INTEGER :: iover_l(0:20), iover_s(0:11, 3) 
INTEGER :: iscfac_bits, ihuff_bits, ibits_best, iscale, iwk_best(SIZE(iwk))
REAL (KIND = 8)     :: wk(SIZE(iwk)), th(SIZE(iwk)), sc(SIZE(iwk))
REAL (KIND = 8)     :: distortion, distortion_min
TYPE (side_info_sub):: side_best
TYPE (scale_factor) :: sc_fac_best
LOGICAL             :: qexit, qfirst
distortion_min = 1.0d10
scale0:DO iscale = 0, 1                  ! scalefactor_scale LOOP                          ! ISO C.1.5.4.3
 side%iscalefac_scale = iscale           ! 0: sqrt(2) or 1: 2
 side%ipreflag = 0                       ! start with preemphasis off                      ! ISO C.1.5.4.3.4
 side%isubblock_gain = 0
 DO                                      ! subblock_gain LOOP (only for short/mixed block) ! ISO 2.4.3.4.7.1
  qfirst = .TRUE.
  sc_fac%long = 0                        ! clear scale factor 
  sc_fac%ishort = 0                      ! clear scale factor
  DO                                     ! scalefactor LOOP
   CALL select_compress(iblock_type, sc_fac  , side%iscalefac_compress)
   CALL calc_scfac_bit (iblock_type, iscfac_bits, side%iscalefac_compress) 
   CALL reorder(iblock_type, r_mdct , wk, icut)
   CALL reorder(iblock_type, x_noise, th,   32)
   CALL calc_scale(iblock_type, iscale, sc_fac, npretab(:, side%ipreflag), side%isubblock_gain, sc)
   wk = sc * wk 
   th = sc * th 
   CALL inner_loop(iallowed_bits - iscfac_bits, iblock_type, wk, iwk, side, ihuff_bits) ! ISO C.1.5.4.3.2
   CALL calc_distortion(iblock_type, side%iglobal_gain, &                               ! ISO C.1.5.4.3.3 
                        wk, iwk, th, sc, iover_l, iover_s, distortion)
   IF (distortion < distortion_min) THEN      ! save best parameters                    ! ISO C.1.5.4.3.1
    distortion_min = distortion
    sc_fac_best    = sc_fac                   ! scale factors
    side_best      = side                     ! side informations
    iwk_best       = iwk                      ! 576 quantized data
    ibits_best     = iscfac_bits + ihuff_bits ! required bits for scale_factors & Huffman codes
   ELSE
    IF (distortion > skip * distortion_min) EXIT ! short cut to avoid meaningless search
    IF (SUM(iover_l) == 20) EXIT ! short cut :long-block
    IF (SUM(iover_s) == 36) EXIT ! short cut :short-block
    IF (SUM(iover_l) ==  9 .AND. SUM(iover_s) == 27) EXIT ! short cut :mixed-block
   END IF
   IF ( SUM(iover_l) + SUM(iover_s) == 0 ) EXIT scale0 ! if converged RETURN            ! ISO C.1.5.4.3.6
   IF ( qfirst .AND. SUM(iover_l(17:20)) == 4 ) THEN                                    ! ISO C.1.5.4.3.4
    side%ipreflag = 1                         ! restart inner_loop with preemphasis on  
    qfirst        = .FALSE.
    CYCLE
   ELSE
    qfirst        = .FALSE.
   END IF
   CALL increase_scale_factor(iblock_type, iover_l, iover_s, sc_fac, qexit)             ! ISO C.1.5.4.3.5  
   IF (qexit) EXIT                           ! when scale_factor reached maximum value defined by ISO 
  END DO
  CALL calc_subblock_gain(iblock_type, iover_s, side%isubblock_gain, qexit) 
  IF (qexit) EXIT                           ! if no subblock_gain increased exit       ! ISO C.1.5.4.3.6  
 END DO 
END DO scale0
distortion = distortion_min                  ! retrieve best parameters and return
sc_fac     = sc_fac_best
side       = side_best 
iwk        = iwk_best
iused_bits = ibits_best 
! debug info
distortion_max = MAX(distortion_max, distortion)
if ( side%iscalefac_scale == 1 ) n_scale = n_scale + 1
if ( side%ipreflag        == 1 ) n_emph  = n_emph  + 1
if ( any(sc_fac%long     /= 0) ) n_sc_l  = n_sc_l  + 1
if ( any(sc_fac%ishort   /= 0) ) n_sc_s  = n_sc_s  + 1
if ( sum(side%isubblock_gain) /= 0 ) n_sub_gain = n_sub_gain + 1
tot_sc_l = tot_sc_l + REAL(sc_fac%long(0:20)       , KIND = 8)
tot_sc_s = tot_sc_s + REAL(sc_fac%ishort(0:11, 1:3), KIND = 8)
RETURN
END SUBROUTINE outer_loop
!------------------------------------------------------------------------------------------------
SUBROUTINE select_compress(iblock_type, sc_fac, icompress)  ! ISO 2.4.3.4.5
IMPLICIT NONE
INTEGER            , INTENT(IN ) :: iblock_type
TYPE (scale_factor), INTENT(IN ) :: sc_fac
INTEGER            , INTENT(OUT) :: icompress
icompress = 0
SELECT CASE (iblock_type) 
 CASE (0, 10, 11, 30, 31) ! long block
  CALL select_compress_long (sc_fac, icompress)
 CASE (20) ! short block
  CALL select_compress_short(sc_fac, icompress)
 CASE (21) ! mixed block
  CALL select_compress_mixed(sc_fac, icompress)
 CASE DEFAULT
  STOP ' error : SUBROUTINE select_compress '
END SELECT
RETURN
END SUBROUTINE select_compress
!..................................................................................................
SUBROUTINE select_compress_long(sc_fac, icompress)     ! ISO 2.4.3.4.5, 2.4.2.7 scalefac_compress
IMPLICIT NONE
TYPE (scale_factor), INTENT(IN ) :: sc_fac
INTEGER            , INTENT(OUT) :: icompress
INTEGER, PARAMETER :: len_scale_compress(0:15, 2) = &  ! 2.4.2.7 scalefac_compress (ISO p.26)
         RESHAPE( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                     0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
INTEGER :: i, n1, n2, iscband 
n1 = 0
n2 = 0
DO iscband =  0, 10
 n1 = MAX(n1, iget_len(sc_fac%long(iscband)) )  
END DO
DO iscband = 11, 20
 n2 = MAX(n2, iget_len(sc_fac%long(iscband)) )  
END DO
DO i = 0, 15
 IF ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) THEN
  icompress = i
  EXIT
 END IF
END DO
IF (n2 >= 4) STOP ' error : select_compress_long ' 
RETURN
END SUBROUTINE select_compress_long
!..................................................................................................
SUBROUTINE select_compress_short(sc_fac, icompress)     ! ISO 2.4.3.4.5, 2.4.2.7 scalefac_compress
IMPLICIT NONE
TYPE (scale_factor), INTENT(IN ) :: sc_fac
INTEGER            , INTENT(OUT) :: icompress
INTEGER :: i, iwin, n1, n2, iscband 
INTEGER, PARAMETER :: len_scale_compress(0:15, 2) = &   ! 2.4.2.7 scalefac_compress (ISO p.26)
         RESHAPE( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                     0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
n1 = 0
n2 = 0
DO iscband =  0, 5
 DO iwin = 1, 3
  n1 = MAX(n1, iget_len(sc_fac%ishort(iscband, iwin)) )  
 END DO
END DO
DO iscband = 6, 11
 DO iwin = 1, 3
  n2 = MAX(n2, iget_len(sc_fac%ishort(iscband, iwin)) )  
 END DO
END DO
DO i = 0, 15
 IF ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) THEN
  icompress = i
  EXIT
 END IF
END DO
IF (n2 >= 4) STOP ' error : select_compress_short ' 
RETURN
END SUBROUTINE select_compress_short
!..................................................................................................
SUBROUTINE select_compress_mixed(sc_fac, icompress)      ! ISO 2.4.3.4.5, 2.4.2.7 scalefac_compress
IMPLICIT NONE
TYPE (scale_factor), INTENT(IN ) :: sc_fac
INTEGER            , INTENT(OUT) :: icompress
INTEGER :: i, iwin, n1, n2, iscband 
INTEGER, PARAMETER :: len_scale_compress(0:15, 2) = &    ! 2.4.2.7 scalefac_compress (ISO p.26)   
         RESHAPE( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                     0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
n1 = 0
n2 = 0
DO iscband =  0, 7
 n1 = MAX(n1, iget_len(sc_fac%long(iscband)) )  
END DO
DO iscband =  3, 5
 DO iwin = 1, 3
  n1 = MAX(n1, iget_len(sc_fac%ishort(iscband, iwin)) )  
 END DO
END DO
DO iscband = 6, 11
 DO iwin = 1, 3
  n2 = MAX(n2, iget_len(sc_fac%ishort(iscband, iwin)) )  
 END DO
END DO
icompress = 0
DO i = 0, 15
 IF ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) THEN
  icompress = i
  EXIT
 END IF
END DO
IF (n2 >= 4) STOP ' error : select_compress_mixed ' 
RETURN
END SUBROUTINE select_compress_mixed
!------------------------------------------------------------------------------------------------
FUNCTION iget_len(k) RESULT(ires)
IMPLICIT NONE
INTEGER, INTENT(IN) :: k
INTEGER :: ires
SELECT CASE (k)
 CASE (0) 
  ires = 0
 CASE (1)
  ires = 1
 CASE (2:3)
  ires = 2
 CASE (4:7)
  ires = 3
 CASE (8:15)
  ires = 4
 CASE DEFAULT
  WRITE(*, *) 'input value', k
  STOP ' error : FUNCTION iget_len ' 
END SELECT
RETURN
END FUNCTION iget_len
!------------------------------------------------------------------------------------------------
SUBROUTINE increase_scale_factor(iblock_type, iover_l, iover_s, sc_fac, qexit) ! C.1.5.4.3.6
IMPLICIT NONE
INTEGER            , INTENT(IN    ) :: iblock_type, iover_l(0:20), iover_s(0:11, 1:3)
TYPE (scale_factor), INTENT(IN OUT) :: sc_fac
LOGICAL            , INTENT(   OUT) :: qexit
INTEGER, PARAMETER :: max_4 = 15, max_3 = 7 ! limit scale factor ! max_4 = 2**4 - 1, max_3 = 2**3 - 1
qexit = .FALSE.
SELECT CASE (iblock_type)
 CASE (0, 10, 11, 30, 31) ! long-block
   WHERE     (sc_fac%long( 0:10)        <  max_4) sc_fac%long( 0:10) = sc_fac%long( 0:10) + iover_l( 0:10)
   WHERE     (sc_fac%long(11:20)        <  max_3) sc_fac%long(11:20) = sc_fac%long(11:20) + iover_l(11:20)
   IF (MAXVAL(sc_fac%long( 0:10))       >= max_4) qexit  = .TRUE.
   IF (MAXVAL(sc_fac%long(11:20))       >= max_3) qexit  = .TRUE.
 CASE (20) ! short-block
   WHERE     (sc_fac%ishort(0: 5, 1:3)  <  max_4) sc_fac%ishort(0: 5, 1:3) = sc_fac%ishort(0: 5, 1:3) + iover_s(0: 5, 1:3)
   WHERE     (sc_fac%ishort(6:11, 1:3)  <  max_3) sc_fac%ishort(6:11, 1:3) = sc_fac%ishort(6:11, 1:3) + iover_s(6:11, 1:3)
   IF (MAXVAL(sc_fac%ishort(0: 5, 1:3)) >= max_4) qexit = .TRUE.
   IF (MAXVAL(sc_fac%ishort(6:11, 1:3)) >= max_3) qexit = .TRUE.
 CASE (21) ! mixed-block
   WHERE     (sc_fac%long(0:7)          <  max_4) sc_fac%long(0:7) = sc_fac%long(0:7) + iover_l(0:7)
   WHERE     (sc_fac%ishort(3: 5, 1:3)  <  max_4) sc_fac%ishort(3: 5, 1:3) = sc_fac%ishort(3: 5, 1:3) + iover_s(3: 5, 1:3)
   WHERE     (sc_fac%ishort(6:11, 1:3)  <  max_3) sc_fac%ishort(6:11, 1:3) = sc_fac%ishort(6:11, 1:3) + iover_s(6:11, 1:3)
   IF (MAXVAL(sc_fac%long(0:7))         >= max_4) qexit = .TRUE.
   IF (MAXVAL(sc_fac%ishort(3: 5, 1:3)) >= max_4) qexit = .TRUE.
   IF (MAXVAL(sc_fac%ishort(6:11, 1:3)) >= max_3) qexit = .TRUE.
 CASE DEFAULT
   STOP ' error : increase_scale_factor'
END SELECT 
RETURN
END SUBROUTINE increase_scale_factor
!------------------------------------------------------------------------------------------------
SUBROUTINE calc_scfac_bit(iblock_type, nbits, icompress)  ! ISO 2.4.3.4.5
IMPLICIT NONE
INTEGER            , INTENT(OUT) :: nbits
INTEGER            , INTENT(IN ) :: iblock_type, icompress
SELECT CASE (iblock_type)
 CASE (0, 10, 11, 30, 31)
  nbits = 11 * len_scale_compress(icompress, 1) + 10 * len_scale_compress(icompress, 2)  
 CASE (20)
  nbits = 18 * (len_scale_compress(icompress, 1) + len_scale_compress(icompress, 2))
 CASE (21)
  nbits = 17 * len_scale_compress(icompress, 1) +  18 * len_scale_compress(icompress, 2)
 CASE DEFAULT
  STOP 'error : SUBROUTINE calc_scfac_bit '
END SELECT
RETURN
END SUBROUTINE calc_scfac_bit 
!------------------------------------------------------------------------------------------------
SUBROUTINE reorder(iblock_type, r_mdct, wk, icut)  ! ISO 2.4.3.4.8
IMPLICIT NONE
INTEGER       , INTENT(IN ) :: iblock_type, icut
REAL (KIND = 8), INTENT(IN ) :: r_mdct(:, :)
REAL (KIND = 8), INTENT(OUT) :: wk(:) 
INTEGER :: k
k = 1
wk = 0.0d0
SELECT CASE (iblock_type)
 CASE (0) ! long-block
  CALL reorder_long (k, 1, icut,        r_mdct, wk)
 CASE (10, 11, 30, 31) ! start / stop
  CALL reorder_long (k, 1, icut,        r_mdct, wk)
 CASE (20) ! short-block
  CALL reorder_short(k, 1, icut, 0, 11, r_mdct, wk)
 CASE (21) ! mixed-block
  CALL reorder_long (k, 1,    2,        r_mdct, wk)
  CALL reorder_short(k, 3, icut, 3, 11, r_mdct, wk) 
 CASE DEFAULT
  STOP ' error : subroutine reorder '
END SELECT 
RETURN
END SUBROUTINE reorder
!..................................................................................................
SUBROUTINE reorder_long(k, n0, n1, r_mdct, wk)
IMPLICIT NONE
INTEGER        , INTENT(IN OUT) :: k
INTEGER        , INTENT(IN    ) :: n0, n1
REAL (KIND = 8), INTENT(IN    ) :: r_mdct(:, :)
REAL (KIND = 8), INTENT(   OUT) :: wk(:) 
INTEGER :: iband
DO iband = n0, n1
 wk(k:k + 18 - 1) = r_mdct(iband, 1:18)
 k = k + 18
END DO
wk(iscalefactorband_l(20, 3) + 2:) = wk(iscalefactorband_l(20, 3) + 2:) * cut_factor
RETURN
END SUBROUTINE reorder_long
!..................................................................................................
SUBROUTINE reorder_short(k, n0, n1, iscfb0, iscfb1, r_mdct, wk)
IMPLICIT NONE
INTEGER        , INTENT(IN OUT) :: k
INTEGER        , INTENT(IN    ) :: n0, n1, iscfb0, iscfb1
REAL (KIND = 8), INTENT(IN    ) :: r_mdct(:, :)
REAL (KIND = 8), INTENT(   OUT) :: wk(:) 
REAL (KIND = 8)                  :: wk0(SIZE(wk) / 3, 3) 
INTEGER :: iwin, k0, m, n, iband, iscfb
k0 = 1
wk0 = 0.0d0
DO iband = n0, n1
 DO iwin = 1, 3
  m = 6 * (iwin - 1) + 1
  wk0(k0:k0 + 6 - 1, iwin) = r_mdct(iband, m:m + 6 - 1)
 END DO
 k0 = k0 + 6
END DO
! Reorder for short-block
k0 = 1
DO iscfb = iscfb0, iscfb1
 n = iscalefactorband_s(iscfb, 1)
 DO iwin = 1, 3
  wk(k:k + n - 1) = wk0(k0:k0 + n - 1, iwin) 
  k = k + n
 END DO
 k0 = k0 + n
END DO
n = (576 - k + 1) / 3
DO iwin = 1, 3
 wk(k:k + n - 1) = wk0(k0:k0 + n - 1, iwin) * cut_factor
 k = k + n
END DO
RETURN
END SUBROUTINE reorder_short
!------------------------------------------------------------------------------------------
SUBROUTINE calc_subblock_gain(iblock_type, iover_s, isubblock_gain, qexit) ! ISO 2.4.2.7 subblock_gain, 2.4.3.4.7.1
IMPLICIT NONE                                                                   ! algorithm UZURA original
INTEGER        , INTENT(IN    ) :: iblock_type, iover_s(0:11, 3)
INTEGER        , INTENT(IN OUT) :: isubblock_gain(:) 
LOGICAL        , INTENT(   OUT) :: qexit
SELECT CASE (iblock_type)
 CASE (0, 10, 11, 30, 31)
  qexit = .TRUE.
 CASE (20) ! short-block
  CALL sub_subblock_gain(iover_s, isubblock_gain, qexit)
 CASE (21) ! mixed-block
  CALL sub_subblock_gain(iover_s, isubblock_gain, qexit) 
 CASE DEFAULT
  STOP ' error : subroutine reorder '
END SELECT 
RETURN
END SUBROUTINE calc_subblock_gain
!..................................................................................................
SUBROUTINE sub_subblock_gain(iover_s, isubblock_gain, qexit)
IMPLICIT NONE
INTEGER        , INTENT(IN    ) :: iover_s(0:11, 3)
INTEGER        , INTENT(IN OUT) :: isubblock_gain(:)
LOGICAL        , INTENT(   OUT) :: qexit
INTEGER :: iwin
qexit = .TRUE.
DO iwin = 1, 3
 IF ( SUM(iover_s(:, iwin)) /= 0 .AND. isubblock_gain(iwin) < 7) THEN
  isubblock_gain(iwin) = isubblock_gain(iwin) + 1
  qexit = .FALSE.
 END IF
END DO
RETURN
END SUBROUTINE sub_subblock_gain
!------------------------------------------------------------------------------------------
SUBROUTINE calc_scale(iblock_type, iscalefac_scale, sc_fac, ipretab, isubblock_gain, scale) ! ISO 2.4.3.4.7.1
IMPLICIT NONE
INTEGER             , INTENT(IN ) :: iblock_type, iscalefac_scale, ipretab(0:20), isubblock_gain(:)
TYPE (scale_factor) , INTENT(IN ) :: sc_fac
REAL (KIND = 8)     , INTENT(OUT) :: scale(:)
INTEGER :: k
scale = 1.0d0
SELECT CASE(iblock_type)
 CASE (0, 10, 11, 30, 31)
  k = 1
  CALL calc_scale_long (k, 0, 20, iscalefac_scale, sc_fac, ipretab       , scale) 
 CASE (20)
  k = 1
  CALL calc_scale_short(k, 0, 11, iscalefac_scale, sc_fac, isubblock_gain, scale) 
 CASE (21)
  k = 1
  CALL calc_scale_long (k, 0,  7, iscalefac_scale, sc_fac, ipretab,        scale) 
  CALL calc_scale_short(k, 3, 11, iscalefac_scale, sc_fac, isubblock_gain, scale) 
 CASE DEFAULT
   STOP ' error : SUBROUTINE rescale '
END SELECT
RETURN 
END SUBROUTINE calc_scale
!..................................................................................................
SUBROUTINE calc_scale_long(k, n0, n1, iscalefac_scale, sc_fac, ipretab, wk)
IMPLICIT NONE
INTEGER             , INTENT(IN OUT) :: k 
INTEGER             , INTENT(IN    ) :: n0, n1, iscalefac_scale, ipretab(0:20)
TYPE (scale_factor) , INTENT(IN    ) :: sc_fac
REAL (KIND = 8)     , INTENT(   OUT) :: wk(:)
INTEGER :: iscband, i
DO iscband = n0, n1
 DO i = 1, iscalefactorband_l(iscband, 1)
  wk(k) = SQRT( REAL( 2**( (1 + iscalefac_scale) * &
                           (sc_fac%long(iscband) + ipretab(iscband)) ), KIND = 8 )  )
  k = k + 1
 END DO
END DO
RETURN
END SUBROUTINE calc_scale_long
!..................................................................................................
SUBROUTINE calc_scale_short(k, n0, n1, iscalefac_scale, sc_fac, isubblock_gain, wk)
IMPLICIT NONE
INTEGER             , INTENT(IN OUT) :: k 
INTEGER             , INTENT(IN    ) :: n0, n1, iscalefac_scale, isubblock_gain(:)
TYPE (scale_factor) , INTENT(IN    ) :: sc_fac
REAL (KIND = 8)     , INTENT(   OUT) :: wk(:)
INTEGER :: iscband, i, iwin
DO iscband = n0, n1
 DO iwin = 1, 3
  DO i = 1, iscalefactorband_s(iscband, 1)
   wk(k) = SQRT( REAL(2**((1 + iscalefac_scale) * sc_fac%ishort(iscband, iwin)), KIND = 8 ) ) & 
         * REAL(4**isubblock_gain(iwin), KIND = 8) 
   k = k + 1
  END DO
 END DO
END DO
RETURN
END SUBROUTINE calc_scale_short
!------------------------------------------------------------------------------------------
SUBROUTINE calc_distortion(iblock_type, iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion) ! ISO C.1.5.4.3.3
IMPLICIT NONE
REAL (KIND = 8) , INTENT(IN ) :: x(:), th(:), sc(:)
INTEGER         , INTENT(IN ) :: iblock_type, iglobal_gain, ix(:)
INTEGER         , INTENT(OUT) :: iover_l(0:), iover_s(0:, :) 
REAL (KIND = 8) , INTENT(OUT) :: distortion 
iover_l = 0
iover_s = 0
SELECT CASE (iblock_type)
 CASE (0, 10, 11, 30, 31)
  CALL calc_dist_long (iglobal_gain, x, ix, th, sc, iover_l, distortion)  
 CASE (20)
  CALL calc_dist_short(iglobal_gain, x, ix, th, sc, iover_s, distortion)  
 CASE (21)
  CALL calc_dist_mixed(iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion)  
 CASE DEFAULT
  STOP ' error : calc_distortion '
END SELECT
RETURN
END SUBROUTINE calc_distortion
!..................................................................................................
SUBROUTINE calc_dist_long(iglobal_gain, x, ix, th, sc, iover_l, distortion)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: x(:), th(:), sc(:)
INTEGER        , INTENT(IN ) :: iglobal_gain, ix(:)
INTEGER        , INTENT(OUT) :: iover_l(0:20)
REAL (KIND = 8), INTENT(OUT) :: distortion
REAL (KIND = 8) ::  dx, ds2(0:20), as2(0:20), d2, a2, bw
INTEGER :: i, istart, iend, iscband
iover_l = 0
ds2 = 0.0d0   
as2 = 0.0d0
distortion = 0.0d0
DO iscband = 0, 20
 bw     = REAL(iscalefactorband_l(iscband, 1), KIND = 8)  
 istart =      iscalefactorband_l(iscband, 2) + 1
 iend   =      iscalefactorband_l(iscband, 3) + 1
 DO i = istart, iend
  dx = ABS(x(i)) - REAL(ABS(ix(i)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)   
  ds2(iscband) = ds2(iscband) + dx   **2.0d0 / bw       
  as2(iscband) = as2(iscband) + th(i)**2.0d0 / bw
  distortion = distortion + ( dx / sc(i) )**2.0d0 / bw
 END DO
END DO
WHERE (ds2 > as2) iover_l = 1
! band 21
d2 = 0.0d0   
a2 = 0.0d0
istart = iscalefactorband_l(20, 3) + 2
iend   = 576
bw     = REAL(iend - istart + 1, KIND = 8)  
DO i = istart, iend
 dx = ABS(x(i)) - REAL(ABS(ix(i)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)   
 d2 = d2 + dx   **2.0d0 / bw       
 a2 = a2 + th(i)**2.0d0 / bw
 distortion = distortion + dx**2.0d0 / bw     
END DO
RETURN
END SUBROUTINE calc_dist_long
!.............................................................................................
SUBROUTINE calc_dist_short(iglobal_gain, x, ix, th, sc, iover_s, distortion)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: x(:), th(:), sc(:)
INTEGER        , INTENT(IN ) :: iglobal_gain, ix(:)
INTEGER        , INTENT(OUT) :: iover_s(0:11, 1:3)
REAL (KIND = 8), INTENT(OUT) :: distortion
REAL (KIND = 8) :: dx, ds2(0:11, 1:3), as2(0:11, 1:3), d2, a2, bw
INTEGER :: i, istart, iend, k, iscband, iwin
k = 0
iover_s = 0
ds2 = 0.0d0
as2 = 0.0d0
distortion = 0.0d0
DO iscband = 0, 11
 bw = REAL(iscalefactorband_s(iscband, 1), KIND = 8)
 istart =  iscalefactorband_s(iscband, 2) + 1
 iend   =  iscalefactorband_s(iscband, 3) + 1  
 DO iwin = 1, 3
  DO i = istart, iend
   k = k + 1
   dx = ABS(x(k)) & 
      - REAL(ABS(ix(k)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)  
   ds2(iscband, iwin) = ds2(iscband, iwin) + dx   **2.0d0 / bw   
   as2(iscband, iwin) = as2(iscband, iwin) + th(k)**2.0d0 / bw   
   distortion = distortion + ( dx / sc(i) )**2.0d0 / bw  
  END DO
 END DO
END DO
WHERE (ds2 > as2) iover_s = 1
! band 12
d2 = 0.0d0   
a2 = 0.0d0
istart = iscalefactorband_s(11, 3) + 2
iend   = 576
bw     = REAL(iend - istart + 1, KIND = 8)  
DO i = istart, iend
 dx = ABS(x(i)) - REAL(ABS(ix(i)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)   
 d2 = d2 + dx   **2.0d0 / bw       
 a2 = a2 + th(i)**2.0d0 / bw
 distortion = distortion + dx**2.0d0 / bw     
END DO
RETURN
END SUBROUTINE calc_dist_short
!..................................................................................................
SUBROUTINE calc_dist_mixed(iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion)
IMPLICIT NONE
REAL (KIND = 8), INTENT(IN ) :: x(:), th(:), sc(:)
INTEGER        , INTENT(IN ) :: iglobal_gain, ix(:)
INTEGER        , INTENT(OUT) :: iover_l(0:20), iover_s(0:11, 1:3)
REAL (KIND = 8), INTENT(OUT) :: distortion
REAL (KIND = 8) :: dx, ds2_l(0:7), as2_l(0:7), ds2_s(3:11, 1:3), as2_s(3:11, 1:3), d2, a2, bw
INTEGER :: i, istart, iend, k, iscband, iwin
k = 0
iover_l = 0
iover_s = 0
ds2_l = 0.0d0
as2_l = 0.0d0
ds2_s = 0.0d0
as2_s = 0.0d0
distortion = 0.0d0
DO iscband = 0, 7
 bw = REAL(iscalefactorband_l(iscband, 1), KIND = 8)
 istart =  iscalefactorband_l(iscband, 2) + 1
 iend   =  iscalefactorband_l(iscband, 3) + 1
 DO i = istart, iend
  k = k + 1
  dx = ABS(x(i)) - REAL(ABS(ix(i)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)  
  ds2_l(iscband) = ds2_l(iscband) + dx   **2.0d0 / bw
  as2_l(iscband) = as2_l(iscband) + th(k)**2.0d0 / bw
  distortion = distortion + ( dx / sc(i) )**2.0d0 / bw
 END DO
 IF (ds2_l(iscband) > as2_l(iscband)) iover_l(iscband) = 1
END DO
DO iscband = 3, 11
 DO iwin = 1, 3 
  bw = REAL(iscalefactorband_s(iscband, 1), KIND = 8)
  istart =  iscalefactorband_s(iscband, 2) + 1
  iend   =  iscalefactorband_s(iscband, 3) + 1  
  DO i = istart, iend
   k = k + 1
   dx = ABS(x(k)) & 
      - REAL(ABS(ix(k)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)  
   ds2_s(iscband, iwin) = ds2_s(iscband, iwin) + dx   **2.0d0 / bw   
   as2_s(iscband, iwin) = as2_s(iscband, iwin) + th(k)**2.0d0 / bw   
   distortion = distortion + ( dx / sc(k) )**2.0d0 / bw  
  END DO
  IF (ds2_s(iscband, iwin) > as2_s(iscband, iwin)) iover_s(iscband, iwin) = 1
 END DO
END DO
! band 12
d2 = 0.0d0   
a2 = 0.0d0
istart = iscalefactorband_s(11, 3) + 2
iend   = 576
bw     = REAL(iend - istart + 1, KIND = 8)  
DO i = istart, iend
 dx = ABS(x(i)) - REAL(ABS(ix(i)), KIND = 8)**(4.0d0 / 3.0d0) * 2.0d0**( REAL(iglobal_gain - 210, KIND = 8) / 4.0d0)   
 d2 = d2 + dx   **2.0d0 / bw       
 a2 = a2 + th(i)**2.0d0 / bw
 distortion = distortion + dx**2.0d0 / bw     
END DO
RETURN
END SUBROUTINE calc_dist_mixed
!------------------------------------------------------------------------------------------------
END MODULE mod_layer3
