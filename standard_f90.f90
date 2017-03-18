MODULE arguments
USE mod_mpg
CONTAINS
!--------------------------------------------------------------------
SUBROUTINE get_option(mpg, fn_in)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(   OUT) :: mpg
CHARACTER (LEN = *)   , INTENT(IN OUT) :: fn_in
mpg%ibit_rate = 9       ! bit_rate 9 = 128kbps
fn_in         = 'fort'  ! input "fn_in.wav"/ output "fn_in.mp3"  
icut          = 25      ! cut subband(icut + 1:)
!qms_stereo    = .TRUE.  ! Mid-Side stereo modo on / off
!q_alias       = .TRUE.  ! take anti-alias in Subband 0-1 in mixed mode
!q_mask        = .TRUE.  ! ATH masking on / off
!q_emphasis    = .TRUE.  ! pre-emphasis
!q_sm          = .TRUE.  ! use both of short & mixed
!q_old_ath     = .FALSE. ! use old ATH function
!q_vbr         = .FALSE. ! VBR mode
!q_rio500      = .FALSE. ! avoid RIO500 VBR bug 
!q_info        = .TRUE.  ! print info at the end 
!mblock_type_param = 20 ! default short block! long = 0, short = 20, mixed = 21
RETURN
END SUBROUTINE get_option
!--------------------------------------------------------------------
END MODULE arguments
!=====================================================================
MODULE raw_io
IMPLICIT NONE
PRIVATE
PUBLIC :: open_wav_file, close_wav_file, open_mpg_file, close_mpg_file, &
          rd_move, rd_word32, rd_word16, rd_int32, rd_int16, wr_cha8 
INTEGER, SAVE :: ir, iw, ir_record = 0, iw_record = 0
LOGICAL, SAVE :: qlittle_endian
CONTAINS
!---------------------------------------------------------------------
SUBROUTINE abort(text)
IMPLICIT NONE
CHARACTER (LEN = *), INTENT(IN) :: text
WRITE(*, *) 'Abort:: ', text
STOP
END SUBROUTINE abort
!--------------------------------------------------------------------
SUBROUTINE check_endian()
IMPLICIT NONE
CHARACTER(LEN = 4) :: ctmp
INTEGER            :: itmp
ctmp = 'ABCD'
itmp = TRANSFER(ctmp, itmp)
IF      ( CHAR( MOD(itmp, 256) ) == 'A' ) THEN
 qlittle_endian = .TRUE.
ELSE IF ( CHAR( MOD(itmp, 256) ) == 'D' ) THEN
 qlittle_endian = .FALSE.
ELSE
 CALL abort('endian indeterminable')
END IF    
RETURN
END SUBROUTINE check_endian
!------------------------------------------------------------------
SUBROUTINE open_wav_file(iread, fname)
IMPLICIT NONE
INTEGER            , INTENT(IN    ) :: iread
CHARACTER (LEN = *), INTENT(IN    ) :: fname
INTEGER :: io
ir = iread
CALL check_endian()
OPEN(ir, FILE = fname, STATUS = 'OLD', IOSTAT = io, ACCESS = 'DIRECT', RECL = 1)  ! RECL = 1byte assumed
IF (io /= 0) THEN
 WRITE(*, '(a, i3, a, i3, 2a)' ) ' I/O error ', io, ' occuerred. file =', iread, ' file name ', fname
 CALL abort('Check input file! Suggestion: Is file name correct?')
END IF
RETURN
END SUBROUTINE open_wav_file
!------------------------------------------------------------------
SUBROUTINE close_wav_file
IMPLICIT NONE
CLOSE(ir)
RETURN
END SUBROUTINE close_wav_file
!------------------------------------------------------------------
SUBROUTINE open_mpg_file(iwrite, fname)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iwrite
CHARACTER (LEN = *) :: fname
INTEGER :: io
iw = iwrite
OPEN(iw, FILE = fname, STATUS = 'REPLACE', IOSTAT = io, ACCESS = 'DIRECT', RECL = 1)
IF (io /= 0) THEN
 WRITE(*, '(a, i3, a, i3, 2a)') ' I/O error ', io, ' occuerred. file =', iw, ' file name ', fname
 CALL abort('Check output file! Suggestion: File may be in use.')
END IF
RETURN
END SUBROUTINE open_mpg_file
!-------------------------------------------------------------------
SUBROUTINE close_mpg_file
IMPLICIT NONE
CLOSE(iw)
RETURN
END SUBROUTINE close_mpg_file
!------------------------------------------------------------------
SUBROUTINE rd_move(n) 
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
ir_record = ir_record + n
RETURN
END SUBROUTINE rd_move
!------------------------------------------------------------------
FUNCTION rd_word32() RESULT(res)
IMPLICIT NONE
CHARACTER (LEN = 4) :: res
INTEGER             :: i, io
DO i = 1, 4
 ir_record = ir_record + 1
 READ(ir, REC = ir_record, IOSTAT = io) res(i:i)
END DO
RETURN
END FUNCTION rd_word32
!---------------------------------------------------------------------
FUNCTION rd_word16() RESULT(res)
IMPLICIT NONE
CHARACTER (LEN = 4) :: res
INTEGER             :: i, io
DO i = 1, 2
 ir_record = ir_record + 1
 READ(ir, REC = ir_record, IOSTAT = io) res(i:i)
END DO
RETURN
END FUNCTION rd_word16
!---------------------------------------------------------------------
FUNCTION rd_int32() RESULT(ires)
IMPLICIT NONE
CHARACTER (LEN = 4) :: tmp
INTEGER :: i, io, ires
IF (qlittle_endian) THEN
 DO i = 1, 4
  ir_record = ir_record + 1
  READ(ir, REC = ir_record, IOSTAT = io) tmp(i:i)
 END DO
ELSE
 DO i = 4, 1, -1
  ir_record = ir_record + 1
  READ(ir, REC = ir_record, IOSTAT = io) tmp(i:i)
 END DO
END IF
ires = TRANSFER(tmp, ires)
RETURN
END FUNCTION rd_int32
!---------------------------------------------------------------------
FUNCTION rd_int16() RESULT(ires)
IMPLICIT NONE
CHARACTER (LEN = 4) :: tmp
INTEGER :: i, io, ires
IF (qlittle_endian) THEN
 DO i = 3, 4
  ir_record = ir_record + 1
  READ(ir, REC = ir_record, IOSTAT = io) tmp(i:i)
 END DO
ELSE
 DO i = 2, 1, -1
  ir_record = ir_record + 1
  READ(ir, REC = ir_record, IOSTAT = io) tmp(i:i)
 END DO
END IF
ires = TRANSFER(tmp, ires) / 2**16
RETURN
END FUNCTION rd_int16
!---------------------------------------------------------------------
SUBROUTINE wr_cha8(n) 
IMPLICIT NONE
INTEGER, INTENT(IN) :: n  
CHARACTER (LEN = 4) :: cha
INTEGER :: io 
cha = TRANSFER(n, cha)
iw_record = iw_record + 1
IF (qlittle_endian) THEN 
 WRITE(iw, REC = iw_record, IOSTAT = io) cha(1:1)
ELSE
 WRITE(iw, REC = iw_record, IOSTAT = io) cha(4:4)
END IF
RETURN
END SUBROUTINE wr_cha8
!---------------------------------------------------------------------
END MODULE raw_io
!==================================================================================
MODULE wav_io 
USE mod_mpg
USE raw_io
IMPLICIT NONE
PRIVATE
PUBLIC :: riff_chunk, fmt_chunk, fact_chunk, data_chunk              ! type
PUBLIC :: check_riff_chunk, open_wav_file, close_wav_file, read_pcm_1frame, read_pcm0 ! subroutine
!
TYPE :: fmt_chunk
 CHARACTER (LEN = 4):: chunk_id  != 'fmt '
 INTEGER :: ichunk_size
 INTEGER :: iformat_type, ichannels
 INTEGER :: isamples_per_sec
 INTEGER :: ibytes_per_sec
 INTEGER :: iblock_size, ibits_per_sample
END TYPE fmt_chunk
!
TYPE :: fact_chunk
 CHARACTER (LEN = 4) :: chunk_id  != 'fact'
 INTEGER :: ichunk_size !( bytes )
!! same unknown data
END TYPE fact_chunk
!
TYPE :: data_chunk
 CHARACTER (LEN = 4) :: chunk_id  != 'data'
 INTEGER :: ichunk_size !( bytes )
!! PCM data follows INTEGER (KIND = 2) :: idata(isize / 2)
END TYPE data_chunk
!
TYPE :: riff_chunk
 CHARACTER (LEN = 4) :: chunk_id  != 'RIFF'
 INTEGER :: ichunk_size
 CHARACTER (LEN = 4) :: format_type ! = 'WAVE'
 TYPE ( fmt_chunk) :: fmt
 TYPE (fact_chunk) :: fct
 TYPE (data_chunk) :: dat
END TYPE riff_chunk
!
CONTAINS
!--------------------------------------------------------------------
SUBROUTINE init_chunk_names(riff)
IMPLICIT NONE
TYPE (riff_chunk), INTENT(OUT) :: riff
riff%chunk_id     = 'RIFF'
riff%format_type  = 'WAVE'
riff%fmt%chunk_id = 'fmt '
riff%fct%chunk_id = 'fact'
riff%dat%chunk_id = 'data'
RETURN
END SUBROUTINE init_chunk_names
!------------------------------------------------------------------
SUBROUTINE check_riff_chunk(riff)
IMPLICIT NONE
TYPE (riff_chunk), INTENT(IN OUT) :: riff
CALL init_chunk_names(riff)
IF ( rd_word32() == riff%chunk_id ) THEN    ! 'RIFF'?
 WRITE(*, '(a)', ADVANCE = 'NO') ' MS RIFF '
ELSE
 CALL abort('This is not MS-RIFF file!')
END IF
riff%ichunk_size = rd_int32()
IF ( rd_word32() == riff%format_type ) THEN ! 'WAVE'?
 WRITE(*, '(a)', ADVANCE = 'NO') 'WAV audio '
ELSE
 WRITE(*, *)
 CALL abort('This is not WAV file!')
END IF
CALL check_fmt_chunk(riff%fmt)
CALL check_fct_chunk(riff%fct)
CALL check_dat_chunk(riff%dat)
RETURN
END SUBROUTINE check_riff_chunk
!------------------------------------------------------------------
SUBROUTINE check_fmt_chunk(fmt)
IMPLICIT NONE
TYPE (fmt_chunk), INTENT(IN OUT) :: fmt
IF ( rd_word32() /= fmt%chunk_id ) CALL abort('Cannot find format chunk!')
fmt%ichunk_size      = rd_int32()
fmt%iformat_type     = rd_int16()
fmt%ichannels        = rd_int16()
fmt%isamples_per_sec = rd_int32()
fmt%ibytes_per_sec   = rd_int32()
fmt%iblock_size      = rd_int16()
fmt%ibits_per_sample = rd_int16()
IF ( fmt%iformat_type     /=  1) CALL abort('Unknown WAVE format!') !linear PCM
IF ( fmt%ibits_per_sample /= 16) CALL abort('Not 16bit data!')
SELECT CASE ( fmt%ichannels )
 CASE (1)
  WRITE(*, '(a, i3, a, i6, a)', ADVANCE = 'NO') &
   'Monoral', fmt%ibits_per_sample, 'bit Sampling rate', fmt%isamples_per_sec, 'Hz '
 CASE (2)
  WRITE(*, '(a, i3, a, i6, a)', ADVANCE = 'NO') &
   'Stereo' , fmt%ibits_per_sample, 'bit Sampling rate', fmt%isamples_per_sec, 'Hz '
 CASE DEFAULT
  WRITE(*, '(a, i1)') ' Number of wave channels is ', fmt%ichannels
  CALL abort('Wave channel must be 1 or 2!')
END SELECT
RETURN
END SUBROUTINE check_fmt_chunk
!------------------------------------------------------------------
SUBROUTINE check_fct_chunk(fct)
IMPLICIT NONE
TYPE (fact_chunk), INTENT(IN OUT) :: fct
IF ( rd_word32() == fct%chunk_id ) THEN 
 fct%ichunk_size = rd_int32()
 CALL rd_move(fct%ichunk_size) ! skip to next chunk 
ELSE
 CALL rd_move(-4) ! back to beginning of this chunk
END IF
RETURN
END SUBROUTINE check_fct_chunk
!------------------------------------------------------------------
SUBROUTINE check_dat_chunk(dat)
IMPLICIT NONE
TYPE (data_chunk), INTENT(IN OUT) :: dat
IF ( rd_word32() /= dat%chunk_id ) CALL abort('Cannot find data chunk!')
dat%ichunk_size = rd_int32()
RETURN
END SUBROUTINE check_dat_chunk
!------------------------------------------------------------------
SUBROUTINE wav_read(pcm)
IMPLICIT NONE
REAL (KIND = 8), INTENT(OUT) :: pcm(:, :)
REAL (KIND = 8), PARAMETER   :: denom = 32768.0d0 !32768 = 2^15
INTEGER        , PARAMETER   :: maxbuff = 1152 * 2
INTEGER :: i, nchannel, ndat
ndat     = SIZE(pcm, 1)
nchannel = SIZE(pcm, 2)
IF (ndat * nchannel > maxbuff) CALL abort('check maxbuff: subroutine wav_get')
SELECT CASE (nchannel)
 CASE (1) !mono
  DO i = 1, ndat
   pcm(i, 1) = REAL( rd_int16(), KIND = 8) / denom 
   pcm(i, 2) = 0.0d0
  END DO
 CASE (2) !stereo
  DO i = 1, ndat
   pcm(i, 1) = REAL( rd_int16(), KIND = 8) / denom  
   pcm(i, 2) = REAL( rd_int16(), KIND = 8) / denom  
  END DO
 CASE DEFAULT
  CALL abort('ichannel must be 1 or 2: subroutine wav_get')
END SELECT
RETURN
END SUBROUTINE wav_read
!------------------------------------------------------------------
SUBROUTINE read_pcm_1frame(pcm)
IMPLICIT NONE
REAL (KIND = 8), INTENT(OUT) :: pcm(:, :)
pcm = EOSHIFT(pcm, 1152, 0.0d0, 1)
CALL wav_read(pcm(481:1632, :))
RETURN
END SUBROUTINE read_pcm_1frame
!------------------------------------------------------------------
SUBROUTINE read_pcm0(pcm)
IMPLICIT NONE
REAL (KIND = 8), INTENT(OUT) :: pcm(:, :)
CALL wav_read(pcm(1153:1632, :))
RETURN
END SUBROUTINE read_pcm0
!------------------------------------------------------------------
END MODULE wav_io
!======================================================================
MODULE bit_io
USE raw_io
USE mod_mpg
USE mod_crc
USE mod_layer3
USE mod_huffman
IMPLICIT NONE
INTEGER,SAVE :: ip
CHARACTER (LEN = 10000) :: bit_string
CONTAINS
!---------------------------------------------------------------------
SUBROUTINE clear_bit_buff()
IMPLICIT NONE
ip = 1
bit_string = REPEAT(' ', LEN(bit_string))
RETURN
END SUBROUTINE clear_bit_buff
!-------------------------------------------------------------------
SUBROUTINE put_bits(n, inp)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, inp
INTEGER :: i, m
IF (n > 32) CALL abort('out of range: n must be < 32: subroutine put_bits')
DO i = 1, n
 m = 2**(n - i)
 IF (MOD(inp / m, 2) == 1) THEN
  bit_string(ip:ip) = '1'
 ELSE
  bit_string(ip:ip) = '0'
 END IF
 ip = ip + 1
END DO
RETURN
END SUBROUTINE put_bits
!-------------------------------------------------------------------
SUBROUTINE put_bits_dim(n, inp)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, inp(:)
INTEGER :: i
DO i = 1, SIZE(inp)
 CALL put_bits(n, inp(i))
END DO
RETURN
END SUBROUTINE put_bits_dim
!-------------------------------------------------------------------
SUBROUTINE put_bits_dim2(n, inp)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, inp(:, :)
INTEGER :: i, j
DO i = 1, SIZE(inp, 1)
 DO j = 1, SIZE(inp, 2)
  CALL put_bits(n, inp(i, j))
 END DO
END DO
RETURN
END SUBROUTINE put_bits_dim2
!-------------------------------------------------------------------
SUBROUTINE put_bits_c(str)
IMPLICIT NONE
CHARACTER (LEN = *) :: str
INTEGER :: i
DO i = 1, LEN_TRIM(str)
 IF (str(i:i) /= '0' .AND. str(i:i) /= '1') CALL abort('invalid string: subroutine put_bit_c')
 bit_string(ip:ip) = str(i:i)
 ip = ip + 1
END DO
RETURN
END SUBROUTINE put_bits_c
!-------------------------------------------------------------------
SUBROUTINE write_bits_1frame(n)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: i, j, ipos, m
IF (MOD(n, 8) /= 0) CALL abort('input error: n must be multiples of 8: SUBROUTINE write_bits_1frame')
ipos = 0
DO i = 1, n, 8
 m = 0
 DO j = 1, 8
  ipos = ipos + 1
  IF (ipos > LEN(bit_string)) EXIT
  IF (bit_string(ipos:ipos) == '1') m = m + 2**(8 - j)
 END DO
 CALL wr_cha8(m)
END DO
RETURN
END SUBROUTINE write_bits_1frame
!-------------------------------------------------------------------
END MODULE bit_io

