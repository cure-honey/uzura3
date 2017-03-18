PROGRAM uzura3
USE arguments
USE wav_io
USE bit_io
USE mod_mpg
USE mod_polyphase
USE mod_psycho
USE mod_mdct
USE mod_layer3
USE mod_huffman
USE mod_encode
IMPLICIT NONE
TYPE (riff_chunk) :: riff
TYPE (mpeg_parameters) :: mpg
REAL (KIND = 8), ALLOCATABLE :: pcm(:, :), subband(:, :, :), r_mdct(:, :, :, :)
INTEGER        , ALLOCATABLE :: i_mdct(:, :, :), mblock_type(:, :), mblock_prev(:)
INTEGER :: iframe_length, ianc
INTEGER :: nchannel, iframe, itotal_frames, i
CHARACTER (LEN = 40) :: file_name, fn_in, fn_out
!
mpg%mtype           =  3 ! 0:mpeg2.5, 1:---, 2:mpeg2, 3:mpeg1
mpg%layer           =  1 ! layer { 1 = III, 2 = II, 3 = I }
mpg%ibit_rate       =  9 ! 128kbps
mpg%isample_rate    =  0
mpg%ipadding        =  0
mpg%icrc            =  1 ! CRC16  0 enabled / 1 disabled
mpg%iprivate        =  0
mpg%mode            =  1 ! joint stereo
mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
mpg%icopyright      =  0
mpg%ioriginal       =  0
mpg%iemphasis       =  0 ! 
CALL init_huffman()
CALL init_huffman_quadruples()
CALL get_option(mpg, file_name)
fn_in  = TRIM(file_name) // '.wav'
fn_out = TRIM(file_name) // '.mp3'
CALL pr_info(mpg)
CALL open_wav_file(10, fn_in)
CALL check_riff_chunk(riff)
CALL play_time(riff)
CALL sample_rate(riff, mpg)
CALL init_scalefactor_bands(mpg%isample_rate)
nchannel = riff%fmt%ichannels
itotal_frames = riff%dat%ichunk_size / (mpeg_frame_size(mpg%layer) * nchannel * 2)   ! ISO 2.4.2.3 padding bit (p22)
IF ( MOD(riff%dat%ichunk_size, mpeg_frame_size(mpg%layer) * nchannel * 2) > 480 ) itotal_frames = itotal_frames + 1
ALLOCATE( pcm(1632, 2), subband(32, 36, nchannel) )
ALLOCATE( r_mdct(32, 18, 2, nchannel), i_mdct(576, 2, nchannel), mblock_type(2, nchannel), mblock_prev(2) )
CALL open_mpg_file(9, fn_out)
mblock_type = mblock_type_param
pcm = 0.0d0
iframe = 0
CALL read_pcm0(pcm)                ! read 480 pcm data (arbitrary)
DO WHILE (iframe < itotal_frames)  ! main loop
 iframe = iframe + 1
 CALL read_pcm_1frame(pcm)
 CALL psycho(pcm, mpg, mblock_type)
 CALL polyphase_filter36(pcm, subband)
 CALL sub_mdct(subband, r_mdct, mblock_type(:, 1), q_alias)
 CALL alloc_bits(mblock_type, r_mdct, i_mdct, mpg, iframe_length, ianc)
 CALL encode_all(mpg, i_mdct, nchannel, ianc)
 CALL write_bits_1frame(iframe_length)
 IF ( MOD(iframe, 50) == 1 ) CALL update_status(iframe, itotal_frames) 
END DO
CALL update_status(iframe, itotal_frames) 
WRITE(*, *) 'total frames', iframe
CALL close_wav_file()
CALL close_mpg_file()
WRITE(*, '(a)') ' Normal End. '
CALL print_debug_info()
STOP
CONTAINS
!------------------------------------------------------------------
SUBROUTINE print_debug_info()
IMPLICIT NONE
if (q_info) then
print '(a)', ' ======== info ==============================================================='
print '(5(a, i7))', ' Block type:long', long, ':short', nshort, ':mixed', mix, ':type1', m1, ':type3', m3
print '(a, 2i7, a, 2i6, a, 3i7)', ' MS/NS select   ', ms, ns, ' [', ns1, ns2, '] :Long/Short switch', nn1, nn2
print '(a)', '..............................................................................'
print '(a, i7)', ' Average scale factor (long 0-20)', n_sc_l
print '(11F6.2)', tot_sc_l( 0: 9)  / REAL( MAX(n_sc_l, 1), KIND = 8)  
print '(11F6.2)', tot_sc_l(10:20)  / REAL( MAX(n_sc_l, 1), KIND = 8) 
print '(a, i7)', ' Average scale factor (short 0-11)', n_sc_s
print '(12F6.2)', ( tot_sc_s(:, 1) + tot_sc_s(:, 2) + tot_sc_s(:, 3) ) / REAL( MAX(3 * n_sc_s, 1), KIND = 8) 
print '(a)', ' Selected Huffman table (table 0-31) ' 
print '( 4(1x, 10i7/) )', ntable 
print '( (1x, a, 10i7/) )', 'Selected count1 table A, B', ntab_ab 
print '(a)', '..............................................................................'
print '(3(a,  i6))', ' pre-emphasis ', n_emph, ' :scalefactor_scale ', n_scale, ' :subblock_gain ', n_sub_gain 
print '(a)', '..............................................................................' 
print '(a, 14i5  )', ' bit:', (i, i = 1, 14)
print '(a, 14F5.1)', ' (%):', REAL(nbits * 100) / REAL(SUM(nbits))
print '(a, F7.2, a, F10.5)', ' average bit rate (kbps)', SUM(REAL(nbits * mpeg_bit_rates(1:14,mpg%layer))) / REAL(SUM(nbits)) &
                            , '      :maximum distortion ', distortion_max
print '(a)', ' ======== info ==============================================================='
end if
RETURN
END SUBROUTINE print_debug_info
!------------------------------------------------------------------
SUBROUTINE sample_rate(riff, mpg)
IMPLICIT NONE
TYPE (riff_chunk     ), INTENT(IN ) :: riff
TYPE (mpeg_parameters), INTENT(OUT) :: mpg
SELECT CASE (riff%fmt%isamples_per_sec)
 CASE (44100)
  mpg%isample_rate = 0
 CASE (48000)
  mpg%isample_rate = 1
 CASE (32000)
  mpg%isample_rate = 2
 CASE DEFAULT
  WRITE(*, *) 'sampling rate ', riff%fmt%isamples_per_sec, ' is not supported in MPEG-1.'
  STOP 'abnormal end'
END SELECT
RETURN
END SUBROUTINE sample_rate
!------------------------------------------------------------------
SUBROUTINE play_time(riff)
IMPLICIT NONE
TYPE (riff_chunk), INTENT(IN) :: riff
INTEGER :: itot_time, ihour, imin, isec
itot_time = riff%dat%ichunk_size / riff%fmt%ibytes_per_sec
ihour =          itot_time / 3600
imin  =      MOD(itot_time, 3600) / 60
isec  = MOD( MOD(itot_time, 3600) , 60 )
WRITE(*, '(a, i3, a, i2, a, i2)') ' Playtime ', ihour, ':', imin, ':', isec
WRITE(*, *)
RETURN
END SUBROUTINE play_time
!------------------------------------------------------------------
SUBROUTINE pr_info(mpg)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(IN) :: mpg
WRITE(*, *) 'UZURA3 (MP3 Encoder/FORTRAN90) Ver.0.5b (c) H.O. Psychoacoustic model Enoken'
IF (mpg%icrc == 0) WRITE(*, *) 'CRC16 error protection enabled'
RETURN
END SUBROUTINE pr_info
!------------------------------------------------------------------
SUBROUTINE update_status(iframe, itot_frames)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iframe, itot_frames
INTEGER :: it(8), ielapsed, iel_min, iel_sec
INTEGER, SAVE :: istart
LOGICAL :: qfirst = .TRUE.
REAL    :: percent
CHARACTER (LEN = 10) :: time, date, zone
CALL DATE_AND_TIME(date, time, zone, it)
IF (qfirst) THEN
 istart   = it(5) * 3600 + it(6) * 60 + it(7)
 qfirst   = .FALSE.
END IF
ielapsed = it(5) * 3600 + it(6) * 60 + it(7) - istart
iel_min  =     ielapsed / 60
iel_sec  = MOD(ielapsed , 60)
percent = REAL(100 * iframe) / REAL(itot_frames)
WRITE(*, '(a, f6.2, a, i4, 2(a, i2), 3(a, i2.2), a, i4.2, a, i2.2, a)')  &
      '+Processed...', percent, '%  ', &
      it(1), '/', it(2), '/', it(3), ' ', it(5), ':', it(6), ':', it(7), &
      ' time elapsed ', iel_min, 'min ', iel_sec, 'sec'
RETURN
END SUBROUTINE update_status
!----------------------------------------------------------------------------
END PROGRAM uzura3
