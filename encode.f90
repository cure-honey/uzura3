MODULE mod_encode
USE bit_io
USE mod_mpg
USE mod_layer3
USE mod_huffman
USE mod_crc
PRIVATE
PUBLIC encode_all
CONTAINS
!-------------------------------------------------------------------
SUBROUTINE encode_all(mpg, ix, nchannel, ianc) ! ISO 2.4.1
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN) :: mpg
INTEGER               , INTENT(IN) :: ix(:, :, :), nchannel, ianc
CALL clear_bit_buff  ()  ! bit strings are written to bit_string (PRIVATE) in MODULE bit_io
CALL encode_header   (mpg) 
CALL encode_crc      (mpg, nchannel)
CALL encode_side_info(     nchannel) ! side_info, scalefactor
CALL encode_part2_3  (ix,  nchannel)
CALL encode_ancillary(ianc)
RETURN
END SUBROUTINE encode_all
!-------------------------------------------------------------------
SUBROUTINE encode_header(mpg)         ! ISO 2.4.1.3 
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN) :: mpg
CALL put_bits_c('11111111111'      )  !sync word
CALL put_bits(2, mpg%mtype         )  !mpeg1
CALL put_bits(2, mpg%layer         )  !layer 
CALL put_bits(1, mpg%icrc          )  !CRC check 
CALL put_bits(4, mpg%ibit_rate     )  !bitrate 
CALL put_bits(2, mpg%isample_rate  )  !sampling frequency 44.1
CALL put_bits(1, mpg%ipadding      )  !ipadding
CALL put_bits(1, mpg%iprivate      )  !private bit : unused
CALL put_bits(2, mpg%mode          )  !stereo
CALL put_bits(2, mpg%mode_extension)  !mode
CALL put_bits(1, mpg%icopyright    )
CALL put_bits(1, mpg%ioriginal     )
CALL put_bits(2, mpg%iemphasis     )
RETURN
END SUBROUTINE encode_header
!-------------------------------------------------------------------
SUBROUTINE encode_crc(mpg, nchannel) ! ISO 2.4.3.1, Table A.9, Table B.5 
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN) :: mpg
INTEGER, INTENT(IN) :: nchannel
INTEGER :: ichannel, igranule, icrc
IF (mpg%icrc == 0) THEN ! if CRC is on
 icrc = 65535           ! Z'0000FFFF' CRC16 initial value
 CALL crc16(4, mpg%ibit_rate     , icrc)
 CALL crc16(2, mpg%isample_rate  , icrc)
 CALL crc16(1, mpg%ipadding      , icrc)
 CALL crc16(1, mpg%iprivate      , icrc)
 CALL crc16(2, mpg%mode          , icrc)
 CALL crc16(2, mpg%mode_extension, icrc)
 CALL crc16(1, mpg%icopyright    , icrc)
 CALL crc16(1, mpg%ioriginal     , icrc)
 CALL crc16(2, mpg%iemphasis     , icrc)
! 
 CALL crc16(9, side_info%main_data_begin, icrc)
 SELECT CASE (nchannel)
  CASE (1) ! mono
   CALL crc16(5, side_info%iprivate_bits, icrc)
  CASE (2) ! stereo
   CALL crc16(3, side_info%iprivate_bits, icrc)
  CASE DEFAULT
   STOP 'illeagal input nchannel: SUBROUTINE encode_crc '
 END SELECT
 DO ichannel = 1, nchannel
  CALL crc16dim(1, side_info%iscfsi(1:4, ichannel), icrc)
 END DO
 DO igranule = 1, 2
  DO ichannel = 1, nchannel
   CALL  crc16(12, side_info%sub(igranule, ichannel)%ipart2_3_length       , icrc)
   CALL  crc16( 9, side_info%sub(igranule, ichannel)%ibig_values           , icrc)
   CALL  crc16( 8, side_info%sub(igranule, ichannel)%iglobal_gain          , icrc)
   CALL  crc16( 4, side_info%sub(igranule, ichannel)%iscalefac_compress    , icrc)
   CALL  crc16( 1, side_info%sub(igranule, ichannel)%iwindow_switching_flag, icrc)
   IF (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1) THEN ! short/mixed-block 
    CALL crc16( 2, side_info%sub(igranule, ichannel)%iblock_type           , icrc)
    CALL crc16( 1, side_info%sub(igranule, ichannel)%mixied_block_flag     , icrc)
    CALL crc16dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:2) , icrc)
    CALL crc16dim( 3, side_info%sub(igranule, ichannel)%isubblock_gain(1:3), icrc)
   ELSE ! long-block  
    CALL crc16dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:3) , icrc)
    CALL crc16( 4, side_info%sub(igranule, ichannel)%iregion0_count        , icrc)
    CALL crc16( 3, side_info%sub(igranule, ichannel)%iregion1_count        , icrc)
   END IF
   CALL  crc16( 1, side_info%sub(igranule, ichannel)%ipreflag              , icrc)
   CALL  crc16( 1, side_info%sub(igranule, ichannel)%iscalefac_scale       , icrc)
   CALL  crc16( 1, side_info%sub(igranule, ichannel)%icount1table_select   , icrc)
  END DO
 END DO
 CALL put_bits(16, icrc) ! write result icrc :16bits
END IF
RETURN
END SUBROUTINE encode_crc
!-------------------------------------------------------------------
SUBROUTINE encode_side_info(nchannel) ! ISO 2.4.1.7
IMPLICIT NONE
INTEGER, INTENT(IN) :: nchannel
INTEGER :: ichannel, igranule
CALL put_bits(9, side_info%main_data_begin)
SELECT CASE (nchannel)
 CASE (1) ! mono
  CALL put_bits(5, side_info%iprivate_bits)
 CASE (2) ! stereo
  CALL put_bits(3, side_info%iprivate_bits)
 CASE DEFAULT
  STOP 'illeagal input nchannel: SUBROUTINE encode_side_info '
END SELECT
DO ichannel = 1, nchannel
 CALL put_bits_dim(1, side_info%iscfsi(1:4, ichannel))
END DO
DO igranule = 1, 2
 DO ichannel = 1, nchannel
  CALL  put_bits(12, side_info%sub(igranule, ichannel)%ipart2_3_length        )
  CALL  put_bits( 9, side_info%sub(igranule, ichannel)%ibig_values            )
  CALL  put_bits( 8, side_info%sub(igranule, ichannel)%iglobal_gain           )
  CALL  put_bits( 4, side_info%sub(igranule, ichannel)%iscalefac_compress     )
  CALL  put_bits( 1, side_info%sub(igranule, ichannel)%iwindow_switching_flag )
  IF (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1) THEN 
   CALL put_bits( 2, side_info%sub(igranule, ichannel)%iblock_type            )
   CALL put_bits( 1, side_info%sub(igranule, ichannel)%mixied_block_flag      )
   CALL put_bits_dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:2) )
   CALL put_bits_dim( 3, side_info%sub(igranule, ichannel)%isubblock_gain(1:3))
  ELSE  
   CALL put_bits_dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:3) )
   CALL put_bits( 4, side_info%sub(igranule, ichannel)%iregion0_count         )
   CALL put_bits( 3, side_info%sub(igranule, ichannel)%iregion1_count         )
  END IF
  CALL  put_bits( 1, side_info%sub(igranule, ichannel)%ipreflag               )
  CALL  put_bits( 1, side_info%sub(igranule, ichannel)%iscalefac_scale        )
  CALL  put_bits( 1, side_info%sub(igranule, ichannel)%icount1table_select    )
 END DO
END DO
RETURN
END SUBROUTINE encode_side_info 
!-------------------------------------------------------------------------------------
SUBROUTINE encode_part2_3(ix, nchannel) ! ISO 2.4.1.7
IMPLICIT NONE
INTEGER, INTENT(IN) :: ix(:, :, :), nchannel
INTEGER :: ichannel, igranule, n
DO igranule = 1, 2
 DO ichannel = 1, nchannel
! scale factors
  IF (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1 .AND. &
      side_info%sub(igranule, ichannel)%iblock_type            == 2         ) THEN
   IF ( side_info%sub(igranule, ichannel)%mixied_block_flag    == 1) THEN ! mixed-block
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
    CALL put_bits_dim(n, scfct(igranule, ichannel)%long(0:7))
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
    CALL put_bits_dim2(n, scfct(igranule, ichannel)%ishort(3:5, 1:3))
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
    CALL put_bits_dim2(n, scfct(igranule, ichannel)%ishort(6:11, 1:3))
   ELSE ! short-block
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
    CALL put_bits_dim2(n, scfct(igranule, ichannel)%ishort(0:5, 1:3))
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
    CALL put_bits_dim2(n, scfct(igranule, ichannel)%ishort(6:11, 1:3))
   END IF
  ELSE ! long block
   IF (side_info%iscfsi(1, ichannel) == 0 .OR. igranule == 1) THEN
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
    CALL put_bits_dim(n, scfct(igranule, ichannel)%long( 0: 5))
   END IF 
   IF (side_info%iscfsi(2, ichannel) == 0 .OR. igranule == 1) THEN 
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
    CALL put_bits_dim(n, scfct(igranule, ichannel)%long(6:10))
   END IF
   IF (side_info%iscfsi(3, ichannel) == 0 .OR. igranule == 1) THEN
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
    CALL put_bits_dim(n, scfct(igranule, ichannel)%long(11:15))
   END IF 
   IF (side_info%iscfsi(4, ichannel) == 0 .OR. igranule == 1) THEN 
    n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
    CALL put_bits_dim(n, scfct(igranule, ichannel)%long(16:20))
   END IF
  END IF
! Huffman codes
  CALL encode_big_region(ix(:, igranule, ichannel), igranule, ichannel)
  CALL encode_quadruples(ix(:, igranule, ichannel), igranule, ichannel)
 END DO
END DO
RETURN
END SUBROUTINE encode_part2_3
!-------------------------------------------------------------------
SUBROUTINE encode_big_region(ix, igranule, ichannel) ! ISO 2.4.1.7
IMPLICIT NONE
INTEGER, INTENT(IN) :: ix(:), igranule, ichannel
INTEGER :: n0, n1, itab, iregion0, iregion1
SELECT CASE (side_info%sub(igranule, ichannel)%iblock_type)  
 CASE (0) ! long-block
   iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
   iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
   n0   = 1
   n1   = iscalefactorband_l(iregion0, 3) + 1 
   itab = side_info%sub(igranule, ichannel)%itable_select(1)
   CALL encode_big(itab, n0, n1, ix)
   n0   = n1 + 1
   n1   = iscalefactorband_l(iregion0 + iregion1 + 1, 3) + 1
   itab = side_info%sub(igranule, ichannel)%itable_select(2)
   CALL encode_big(itab, n0, n1, ix)
   n0   = n1 + 1
   n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
   itab = side_info%sub(igranule, ichannel)%itable_select(3) 
   CALL encode_big(itab, n0, n1, ix)
 CASE (1, 3) ! long-block start/stop
   iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
   iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
   n0   = 1
   n1   = iscalefactorband_l(iregion0 + 1, 2) 
   itab = side_info%sub(igranule, ichannel)%itable_select(1)
   CALL encode_big(itab, n0, n1, ix)
   n0   = n1 + 1
   n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
   itab = side_info%sub(igranule, ichannel)%itable_select(2)
   CALL encode_big(itab, n0, n1, ix)
 CASE (2) ! short-block or mixed-block
  IF (side_info%sub(igranule, ichannel)%mixied_block_flag == 0) THEN ! short block
   iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
   iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
   n0   = 1
   n1   = iscalefactorband_s( (iregion0 + 1) / 3, 2) * 3 
   itab = side_info%sub(igranule, ichannel)%itable_select(1)
   CALL encode_big(itab, n0, n1, ix)
   n0   = n1 + 1
   n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
   itab = side_info%sub(igranule, ichannel)%itable_select(2)
   CALL encode_big(itab, n0, n1, ix)
  ELSE IF (side_info%sub(igranule, ichannel)%mixied_block_flag == 1) THEN ! mixed block
   iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
   iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
   n0   = 1
   n1   = iscalefactorband_l(iregion0, 3) + 1 
   itab = side_info%sub(igranule, ichannel)%itable_select(1)
   CALL encode_big(itab, n0, n1, ix)
   n0   = n1 + 1
   n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
   itab = side_info%sub(igranule, ichannel)%itable_select(2)
   CALL encode_big(itab, n0, n1, ix)
  ELSE
   STOP 'error : encode_big_region : mixed_block_flag '
  END IF
 CASE DEFAULT
  STOP 'error : encode_big_region : iblock_type '
END SELECT
RETURN
END SUBROUTINE encode_big_region
!-------------------------------------------------------------------
SUBROUTINE encode_quadruples(ix, igranule, ichannel) ! ISO 2.4.1.7, 2.4.2.7 huffmancodebits()
IMPLICIT NONE
INTEGER, INTENT(IN) :: ix(:), igranule, ichannel
INTEGER :: i, n0, n1, itab, nbits, k1, k2, k3, k4, is1, is2, is3, is4
n0 =      1 + side_info%sub(igranule, ichannel)%ibig_values * 2 
n1 = n0 - 1 + side_info%sub(igranule, ichannel)%icount1 * 4 
itab = side_info%sub(igranule, ichannel)%icount1table_select
DO i = n0, n1, 4
 k1 = ABS( ix(i    ) )
 k2 = ABS( ix(i + 1) )
 k3 = ABS( ix(i + 2) )
 k4 = ABS( ix(i + 3) )
 is1 = IAND( 1, ISHFTC( ix(i    ), 1 ) ) ! get sign bit 
 is2 = IAND( 1, ISHFTC( ix(i + 1), 1 ) ) ! bit 31 of [31....0]
 is3 = IAND( 1, ISHFTC( ix(i + 2), 1 ) ) ! cyclic shift + AND B'00...001' 
 is4 = IAND( 1, ISHFTC( ix(i + 3), 1 ) ) ! positive 0 / negative 1
 IF      (itab == 0) THEN  ! Table A      ISO Table B.7
  nbits = huff_qa(k1, k2, k3, k4)%leng 
  CALL put_bits(nbits, huff_qa(k1, k2, k3, k4)%icod)
 ELSE IF (itab == 1) THEN  ! Table B      ISO Table B.7
  nbits = huff_qb(k1, k2, k3, k4)%leng 
  CALL put_bits(nbits, huff_qb(k1, k2, k3, k4)%icod)
 ELSE
  STOP 'error '
 END IF 
 IF (k1 /= 0) CALL put_bits(1, is1)
 IF (k2 /= 0) CALL put_bits(1, is2)
 IF (k3 /= 0) CALL put_bits(1, is3)
 IF (k4 /= 0) CALL put_bits(1, is4)
END DO
RETURN
END SUBROUTINE encode_quadruples
!-------------------------------------------------------------------
SUBROUTINE encode_big(itab, n0, n1, ix)  ! ISO 2.4.1.7, 2.4.2.7 huffmancodebits(), C.1.5.3.7
IMPLICIT NONE
INTEGER, INTENT(IN) :: itab, n0, n1, ix(:)
INTEGER :: i, is1, is2, k1, k2, linbitsx, linbitsy
DO i = n0, n1, 2
 IF (itab == 0) EXIT
 k1 = ABS( ix(i    ) )
 k2 = ABS( ix(i + 1) )
 is1 = IAND( 1, ISHFTC(ix(i    ), 1) ) ! get sign bit ! bit 31 of [31..0]  ! positive = 0  
 is2 = IAND( 1, ISHFTC(ix(i + 1), 1) ) ! cyclic shift + AND B'00...001'    ! negative = 1
 IF (itab <=15) THEN
  CALL put_bits( huff(itab)%leng(k1, k2), huff(itab)%icod(k1, k2) )   
  IF (k1 /= 0) CALL put_bits(1, is1)
  IF (k2 /= 0) CALL put_bits(1, is2)
 ELSE
  IF (k1 > 14) THEN 
   linbitsx = k1 - 15
   k1 = 15   
  END IF
  IF (k2 > 14) THEN 
   linbitsy = k2 - 15
   k2 = 15   
  END IF
  CALL put_bits( huff(itab)%leng(k1, k2), huff(itab)%icod(k1, k2) )   
  IF (k1 == 15       ) CALL put_bits(huff(itab)%linbits, linbitsx)  
  IF (ix(i)     /=  0) CALL put_bits(1, is1)
  IF (k2 == 15       ) CALL put_bits(huff(itab)%linbits, linbitsy)  
  IF (ix(i + 1) /=  0) CALL put_bits(1, is2)
 END IF
END DO
RETURN
END SUBROUTINE encode_big
!-------------------------------------------------------------------
SUBROUTINE encode_ancillary(ianc) ! ISO C.1.5.3.6
IMPLICIT NONE
INTEGER, INTENT(IN) :: ianc
INTEGER :: i
DO i = 1, ianc        ! fill remaining bits with 0 
 CALL put_bits(1, 0)
END DO
RETURN
END SUBROUTINE encode_ancillary
!-------------------------------------------------------------------
END MODULE mod_encode