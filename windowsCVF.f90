MODULE arguments 
! module for command line option  ! Compaq (DEC) Visual FORTRAN for Intel (Windows) 
USE DFLIB
USE IFPORT ! Intel Visual Fortran
USE mod_mpg
PRIVATE
PUBLIC :: get_option
CONTAINS
!------------------------------------------------------------------
SUBROUTINE get_option(mpg, fn_in)
IMPLICIT NONE
TYPE (mpeg_parameters), INTENT(   OUT) :: mpg
CHARACTER (LEN = *)   , INTENT(IN OUT) :: fn_in
INTEGER   (KIND = 4) :: narg
INTEGER   (KIND = 2) :: iarg, istatus
CHARACTER (LEN = 40) :: buffer
CHARACTER (LEN =  6) :: fmt
iarg = 0
narg = NARGS()
IF (narg == 1) THEN 
 CALL print_option()
 STOP
END IF 
DO
 iarg = iarg + 1
 IF (iarg >= narg) CALL option_error( TRIM(buffer) )
 CALL GETARG(iarg, buffer)
 IF (buffer(1:1) /= '-') EXIT  
 SELECT CASE(TRIM(buffer))
  CASE ('-b') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(I', istatus, ')' 
   READ(buffer, fmt) mpg%ibit_rate
   IF (mpg%ibit_rate < 1 .OR. mpg%ibit_rate > 14) CALL option_error( TRIM(buffer) )  
  CASE ('-switch') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) switch
  CASE ('-ath_min') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) ath_min
  CASE ('-ath_max') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) ath_max
  CASE ('-crc') 
   mpg%icrc       = 0 ! CRC16 on 
  CASE ('-cut') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(I', istatus, ')' 
   READ(buffer, fmt) icut
   IF (icut < 0 .OR. icut > 32) CALL option_error( TRIM(buffer) )
  CASE ('-xms') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) xms
  CASE ('-xsm') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) xsm
  CASE ('-offset') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) offset
  CASE ('-tempo') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) tempo
  CASE ('-factor') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) factor
   IF (factor < 0.0d0 .OR. factor > 1.0d0) THEN 
    WRITE(*, *) 'input out of range: 0 <= factor <= 1'
    CALL option_error( TRIM(buffer) )
   END IF
  CASE ('-r0') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) r0
  CASE ('-r1') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) r1
  CASE ('-pm') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) pm_factor
  CASE ('-skip') 
   iarg = iarg + 1
   IF ( iarg >= narg ) CALL option_error( TRIM(buffer) )
   CALL GETARG(iarg, buffer, istatus)
   WRITE(fmt, '(a, i1, a)') '(F', istatus, '.0)' 
   READ(buffer, fmt) skip
   IF (skip < 1.0d0) THEN
    WRITE(*, *) 'input out of range: skip >= 1.0'
    CALL option_error( TRIM(buffer) )
   END IF
  CASE ('-cuth') 
   cut_factor = 0.0d0 
  CASE ('-l') 
   mblock_type_param = 0
  CASE ('-s') 
   q_sm = .FALSE.
   mblock_type_param = 20
  CASE ('-m') 
   q_sm = .FALSE.
   mblock_type_param = 21
  CASE ('-sm') 
   q_sm = .TRUE.
  CASE ('-v') ! VBR mode
   q_vbr = .TRUE. 
   qms_stereo = .FALSE.
   icut = 32
  CASE ('-rio500') 
   q_rio500 = .TRUE.
  CASE ('-c') 
   mpg%icopyright = 1 ! copyrigt on
  CASE ('-o') 
   mpg%ioriginal  = 1 ! original on
  CASE ('-ms')
   qms_stereo = .TRUE.
  CASE ('-ns')
   qms_stereo = .FALSE.
   pm_factor = 1.1d0 ! for vbr mode
  CASE ('-nomask')
   q_mask = .FALSE.
  CASE ('-noalias')
   q_alias = .FALSE.
  CASE ('-noinfo')
   q_info = .FALSE.
  CASE ('-about')
   CALL print_about()
  CASE DEFAULT
   CALL option_error(TRIM(buffer))
 END SELECT
END DO
fn_in = TRIM(buffer)
RETURN
END SUBROUTINE get_option
!---------------------------------------------------------------
SUBROUTINE option_error(text)
IMPLICIT NONE
CHARACTER (LEN = *), INTENT(IN) :: text
CALL print_option()
WRITE(*, *)
WRITE(*, '(2a)') ' >>>>>> Error near >>>>>>', text
STOP
END SUBROUTINE option_error
!---------------------------------------------------------------
SUBROUTINE print_option()
IMPLICIT NONE
WRITE(*, *) 'Usage : uzura3 -options file_name '
WRITE(*, *) '      : file_name.wav -> file_name.mp3'
WRITE(*, *) 'Option:-b  1..14   bitrate for CBR mode                 (default 9 : 128kbps)'
WRITE(*, *) '       -crc        CRC16 error protection on            (default off)'
WRITE(*, *) '       -c          copyright flag on                    (default off)'
WRITE(*, *) '       -o          original  flag on                    (default off)'
WRITE(*, *) '       -cut 1..32  band cut-off : place after -b option (default 26: 17.9kHz)'
WRITE(*, *) '       -cuth       cut band 21 (l) or 12 (s/m)          (default off) '
WRITE(*, *) '       -v          VBR mode  (ns, icut = 32)            (default off)'
WRITE(*, *) '       -rio500     avoid RIO500 VBR skip bug            (default off)'
WRITE(*, *) '       -l          long-block-only                      (default off)'
WRITE(*, *) '       -s          short-mode for short-block           (default off)'
WRITE(*, *) '       -m          mixed-mode for short-block           (default off)'
WRITE(*, *) '       -sm         short & mixed-mode for short-block   (default on )'
WRITE(*, *) '       -xsm xx     short / mixed switching parameter    (default 1.5)'
WRITE(*, *) '       -switch xx  long/short switching parameter       (default 2.0)'
WRITE(*, *) '       -skip   xx  speeds up outer loop                 (default 1.3)'
WRITE(*, *) '       -ms/-ns     stereo mode (MS/normal)              (default MS)'
WRITE(*, *) '       -xms xx     MS/NS      switching parameter       (default 0.5)'
WRITE(*, *) '       -nomask     masking off                          (default on)'
WRITE(*, *) '       -ath_min xx minimum of ATH  [ dB ]               (default -125)'
WRITE(*, *) '       -ath_max xx ceiling of ATH  [ dB ]               (default  0.0)'
WRITE(*, *) '       -offset  xx offset for mask [ dB ]               (default 40.0)'
WRITE(*, *) '       -tempo   xx temporal masking factor              (default 0.85)'
WRITE(*, *) '       -factor  xx bit distribution among (gr, ch)      (default 0.4)'
WRITE(*, *) '       -noalias    anti-alias for mixed-block off       (default on)'
WRITE(*, *) '       -debug      print debug info                     (default on) '
WRITE(*, *) '       -about      about UZURA3 '
WRITE(*, *) '        '
WRITE(*, *) 'Example CBR 128kbps CRC on : uzura3 -crc   file_name '
WRITE(*, *) 'Example VBR normal stereo  : uzura3 -v -ns file_name '
RETURN
END SUBROUTINE print_option
!---------------------------------------------------------------
SUBROUTINE print_about()
IMPLICIT NONE
WRITE(*, *) ' ******* UZURA3 is an MPEG-I/Layer-III encoder written in FORTRAN90. ********'
WRITE(*, *) ' **** This program comes with absolutely no warranty. (c) H.O. 2000-2004 **** '
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) ' http://members.tripod.co.jp/kitaurawa/index.html   (Japanese page)'
WRITE(*, *) ' http://members.tripod.co.jp/kitaurawa/index_e.html (English  page)'
STOP 
END SUBROUTINE print_about
!---------------------------------------------------------------
END MODULE arguments
!=============================================================================================
MODULE wav_io
USE mod_mpg
IMPLICIT NONE
PRIVATE
PUBLIC :: riff_chunk, fmt_chunk, data_chunk              ! type
PUBLIC :: check_riff_chunk, open_wav_file, close_wav_file, read_pcm_1frame, read_pcm0 ! subroutine
INTEGER, SAVE :: ir
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
!! some unknown data
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
riff%fct%chunk_id = 'fact '
riff%dat%chunk_id = 'data'
RETURN
END SUBROUTINE init_chunk_names
!---------------------------------------------------------------------
FUNCTION word32() RESULT(res)
IMPLICIT NONE
CHARACTER (LEN = 4) :: res
INTEGER :: io
READ(ir, '(a4)', IOSTAT = io) res
RETURN
END FUNCTION word32
!---------------------------------------------------------------------
FUNCTION word16() RESULT(res)
IMPLICIT NONE
CHARACTER (LEN = 2) :: res
INTEGER :: io
READ(ir, '(a2)', IOSTAT = io) res
RETURN
END FUNCTION word16
!---------------------------------------------------------------------
FUNCTION word8() RESULT(res)
IMPLICIT NONE
CHARACTER (LEN = 1) :: res
INTEGER :: io
READ(ir, '(a1)', IOSTAT = io) res
RETURN
END FUNCTION word8
!---------------------------------------------------------------------
FUNCTION int32() RESULT(ires)
IMPLICIT NONE
INTEGER (KIND = 4) :: ires
ires = TRANSFER(word32(), ires) ! little endian assumed
RETURN
END FUNCTION int32
!---------------------------------------------------------------------
FUNCTION int16() RESULT(ires)
IMPLICIT NONE
INTEGER (KIND = 2) :: ires
ires = TRANSFER(word16(), ires) ! little endian assumed
RETURN
END FUNCTION int16
!---------------------------------------------------------------------
SUBROUTINE abort(text)
IMPLICIT NONE
CHARACTER (LEN = *), INTENT(IN) :: text
WRITE(*, *) 'Abort:: ', text
STOP
END SUBROUTINE abort
!------------------------------------------------------------------
SUBROUTINE open_wav_file(iread, fname)
IMPLICIT NONE
INTEGER            , INTENT(IN    ) :: iread
CHARACTER (LEN = *), INTENT(IN    ) :: fname
INTEGER :: io
ir = iread
OPEN(ir, FILE = fname, STATUS = 'old', IOSTAT = io, RECORDTYPE = 'stream') !non-standard  RECORDTYPE "stream"
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
SUBROUTINE check_riff_chunk(riff)
IMPLICIT NONE
TYPE (riff_chunk), INTENT(IN OUT) :: riff
CALL init_chunk_names(riff)
IF ( word32() == riff%chunk_id ) THEN    ! 'RIFF'?
 WRITE(*, '(a)', ADVANCE = 'NO') ' MS RIFF '
ELSE
 CALL abort('This is not MS-RIFF file!')
END IF
riff%ichunk_size = int32()
IF ( word32() == riff%format_type ) THEN ! 'WAVE'?
 WRITE(*, '(a)', ADVANCE = 'NO') 'WAV audio '
ELSE
 WRITE(*, *)
 CALL abort('This is not WAV file!')
END IF
CALL check_fmt_chunk(riff%fmt)
CALL check_dat_chunk(riff%dat, riff%fct)
RETURN
END SUBROUTINE check_riff_chunk
!------------------------------------------------------------------
SUBROUTINE check_fmt_chunk(fmt)
IMPLICIT NONE
TYPE (fmt_chunk), INTENT(IN OUT) :: fmt
IF ( word32() /= fmt%chunk_id ) CALL abort('Cannot find format chunk!')
fmt%ichunk_size      =     int32()
fmt%iformat_type     = INT(int16(), KIND = 4)
fmt%ichannels        = INT(int16(), KIND = 4)
fmt%isamples_per_sec =     int32()
fmt%ibytes_per_sec   =     int32()
fmt%iblock_size      = INT(int16(), KIND = 4)
fmt%ibits_per_sample = INT(int16(), KIND = 4)
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
SUBROUTINE check_dat_chunk(dat, fct)
IMPLICIT NONE
TYPE (data_chunk), INTENT(IN OUT) :: dat
TYPE (fact_chunk), INTENT(IN OUT) :: fct
INTEGER :: i
CHARACTER (LEN = 4) :: chnk_id
CHARACTER (LEN = 1) :: dummy
chnk_id = word32()
IF      ( chnk_id == fct%chunk_id ) THEN 
 fct%ichunk_size = int32()
 DO i = 1, fct%ichunk_size
  dummy = word8()
 END DO
 IF ( word32() == dat%chunk_id ) THEN
  dat%ichunk_size = int32()
 END IF
ELSE  IF ( chnk_id == dat%chunk_id ) THEN
 dat%ichunk_size = int32()
ELSE
 CALL abort('Cannot find fact chunk nor data chunk!')
END IF
RETURN
END SUBROUTINE check_dat_chunk
!------------------------------------------------------------------
SUBROUTINE wav_read(pcm) ! 16bit PCM assumed
IMPLICIT NONE
REAL (KIND = 8), INTENT(OUT) :: pcm(:, :)
REAL (KIND = 8), PARAMETER   :: denom = 32768.0d0 !32768 = 2^15
INTEGER       , PARAMETER   :: maxbuff = 1152 * 2
CHARACTER (LEN = 2) :: cbuff16(maxbuff)
INTEGER :: i, nchannel, ndat
INTEGER  (KIND = 2) :: ibuff16(maxbuff)
EQUIVALENCE (cbuff16, ibuff16)
ndat     = SIZE(pcm, 1)
nchannel = SIZE(pcm, 2)
IF (ndat * nchannel > maxbuff) CALL abort('check maxbuff: subroutine wav_get')
ibuff16 =0
SELECT CASE (nchannel)
 CASE (1) !mono
  CALL wav_read_sub( cbuff16(1:ndat) )
  DO i = 1, ndat
   pcm(i, 1) = REAL( ibuff16(i), KIND = 8) / denom          ! little endian assumed
   pcm(i, 2) = 0.0d0
  END DO
 CASE (2) !stereo
  CALL wav_read_sub( cbuff16(1:2 * ndat) )
  DO i = 1, ndat
   pcm(i, 1) = REAL( ibuff16(2 * i - 1), KIND = 8) / denom  ! little endian assumed
   pcm(i, 2) = REAL( ibuff16(2 * i    ), KIND = 8) / denom  ! little endian assumed
  END DO
 CASE DEFAULT
  CALL abort('ichannel must be 1 or 2: subroutine wav_get')
END SELECT
RETURN
END SUBROUTINE wav_read
!------------------------------------------------------------------
SUBROUTINE wav_read_sub(cha16)
IMPLICIT NONE
CHARACTER (LEN = 2), INTENT(OUT) :: cha16(:)
INTEGER :: io
READ(ir, '(a2)', iostat = io) cha16
SELECT CASE (io)
 CASE (0)
  CONTINUE
 CASE (-1)
  CONTINUE
!  WRITE(*, *) 'End of File!'
 CASE DEFAULT
  WRITE(*, '(a, i3, a, i4)') ' file ', ir, ' iostat ', io
  CALL abort('I/O error occurred while reading wav file.')
END SELECT
RETURN
END SUBROUTINE wav_read_sub
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
!==================================================================================
MODULE bit_io
USE mod_mpg
USE mod_crc
USE mod_layer3
USE mod_huffman
IMPLICIT NONE
PRIVATE
PUBLIC :: open_mpg_file, close_mpg_file, clear_bit_buff, put_bits, put_bits_c &
        , put_bits_dim, put_bits_dim2, write_bits_1frame
INTEGER, SAVE :: iw, ip
CHARACTER (LEN = 10000) :: bit_string
CONTAINS
!---------------------------------------------------------------------
SUBROUTINE abort(text)
IMPLICIT NONE
CHARACTER (LEN = *), INTENT(IN) :: text
WRITE(*, *) 'Abort:: ', text
STOP
END SUBROUTINE abort
!-------------------------------------------------------------------
SUBROUTINE open_mpg_file(iwrite, fname)
IMPLICIT NONE
INTEGER, INTENT(IN) :: iwrite
CHARACTER (LEN = *) :: fname
INTEGER :: io
iw = iwrite
OPEN(iw, FILE = fname, IOSTAT = io, STATUS = 'unknown', RECORDTYPE = 'stream') ! non-standard RECORDTYPE "stream"
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
!-------------------------------------------------------------------
SUBROUTINE clear_bit_buff()
IMPLICIT NONE
ip = 1                      ! initialize buffer ; set position to 1 
bit_string = REPEAT(' ', LEN(bit_string))
RETURN
END SUBROUTINE clear_bit_buff
!-------------------------------------------------------------------
SUBROUTINE put_bits(n, inp) 
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, inp
INTEGER :: i 
IF (n > 32) CALL abort('out of range: n must be < 32: subroutine put_bits')
DO i = 1, n
 IF (IBITS(inp, n - i, 1) == 1) THEN
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
CHARACTER(LEN = 4) ::cm
EQUIVALENCE (m, cm) ! integer*4 assumed for m
IF (MOD(n, 8) /= 0) CALL abort('input error: n must be multiples of 8: SUBROUTINE write_bits_1frame')
ipos = 0
DO i = 1, n, 8
 m = 0
 DO j = 1, 8
  ipos = ipos + 1
  IF (ipos > LEN(bit_string)) EXIT
  IF (bit_string(ipos:ipos) == '1') m = m + 2**(8 - j)
 END DO
 WRITE(iw, '(a1)', ADVANCE = 'no') cm(1:1)   ! little endian assumed
END DO
RETURN
END SUBROUTINE write_bits_1frame
!-------------------------------------------------------------------
END MODULE bit_io
