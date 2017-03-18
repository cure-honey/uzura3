MODULE mod_mpg
IMPLICIT NONE
PUBLIC
TYPE:: mpeg_parameters
 INTEGER :: mtype
 INTEGER :: layer
 INTEGER :: ibit_rate
 INTEGER :: isample_rate
 INTEGER :: ipadding
 INTEGER :: iprivate
 INTEGER :: icrc
 INTEGER :: mode
 INTEGER :: mode_extension
 INTEGER :: icopyright
 INTEGER :: ioriginal
 INTEGER :: iemphasis
END TYPE mpeg_parameters
!MPEG1 / audio
INTEGER, PARAMETER :: mpeg_frame_size(3)      = (/1152, 1152, 384/)         ! ISO 2.4.2.1 frame
INTEGER, PARAMETER :: mpeg_sample_rates(0:3)  = (/44100, 48000, 32000, 0/)  ! ISO 2.4.2.3 sampling_frequency
INTEGER, PARAMETER :: mpeg_bit_rates(0:14, 3) = &                           ! ISO 2.4.2.3 bitrate_index
  RESHAPE( (/ 0, 32, 40, 48,  56,  64,  80,  96, 112, 128, 160, 192, 224, 256, 320,    &
              0, 32, 48, 56,  64,  80,  96, 112, 128, 160, 192, 224, 256, 320, 384,    &
              0, 32, 64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 414, 448 /), &
           (/15, 3/) )
CHARACTER (LEN = 8) :: mpeg_mode_names(4)      = (/'stereo  ', 'j-stereo', 'dual-ch ', 'mono    '/)  ! ISO 2.4.2.3 mode
CHARACTER (LEN = 3) :: mpeg_layer_names(3)     = (/'III', 'II ', 'I  '/)                             ! ISO 2.4.2.3 Layer
CHARACTER (LEN = 7) :: mpeg_version_names(0:3) = (/'MPEG2.5', '       ', 'MPEG-II', 'MPEG-I '/)      ! ISO 2.4.2.3 ID ! MPEG2.5 non-standard 
CHARACTER (LEN = 7) :: mpeg_demp_names(4)      = (/'none   ', '50/15us', '       ', 'CITT   '/)      ! ISO 2.4.2.3 emphasis
!-------------------------------------------------------------------------------------------
!MPEG1 / Layer3:   ISO 2.4.1.7, 2.4.2.7 
TYPE :: side_info_sub
 INTEGER :: ipart2_3_length        ! 12 bits
 INTEGER :: ibig_values            !  9 bits
 INTEGER :: iglobal_gain           !  8 bits
 INTEGER :: iscalefac_compress     !  4 bits
 INTEGER :: iwindow_switching_flag !  1 bit
! if short block 
 INTEGER ::  iblock_type           !  2 bits
 INTEGER ::  mixied_block_flag     !  1 bit
 INTEGER ::  itable_select(3)      !  5 bits
 INTEGER ::  isubblock_gain(3)     !  3 bits
! if long block
!INTEGER ::  itable_select(3, 2, 2) 
 INTEGER ::  iregion0_count        !  4 bits
 INTEGER ::  iregion1_count        !  3 bits
! 
 INTEGER :: ipreflag               !  1 bit
 INTEGER :: iscalefac_scale        !  1 bit
 INTEGER :: icount1table_select    !  1 bit
!........
 INTEGER :: icount1  ! local use (not real side_info)
END TYPE side_info_sub
!
TYPE :: side_info_                 !136 bits for mono / 256 bits for stereo
 INTEGER :: main_data_begin        !  9 bits
 INTEGER :: iprivate_bits          !  5 bits for mono /   3 bits for stereo
 INTEGER :: iscfsi(4, 2)           !  1 bit 
 TYPE (side_info_sub) :: sub(2, 2)
END TYPE side_info_
!
INTEGER, PARAMETER :: len_scale_compress(0:15, 2) = & ! ISO 2.4.2.7 scalefac_compress[gr][ch]
         RESHAPE( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                     0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
INTEGER, SAVE :: iscalefactorband_l(0:20, 3), iscalefactorband_s(0:11, 3) ! PUBLIC
!-------------------------------------------------------------------------------------------
!global variables
!debug variables
real (KIND = 8), save :: tot_sc_l(0:20) = 0.0d0, tot_sc_s(0:11, 3) = 0.0d0
integer, save :: m1 = 0, m3 = 0, mix = 0, long = 0, nshort = 0, ntable(0:31) = 0, ntab_ab(0:1) = 0
integer, save :: ms = 0, ns = 0, ns1 = 0, ns2 = 0, nn1 = 0, nn2 = 0, nbits(14)
integer, save :: n_sc_l = 0, n_sc_s = 0, n_sub_gain = 0, n_scale = 0, n_emph = 0
!system parameters
! switches
LOGICAL, SAVE :: qms_stereo = .TRUE., q_alias = .TRUE., q_mask = .TRUE., q_sm = .TRUE.
LOGICAL, SAVE :: q_vbr = .FALSE., q_rio500 = .FALSE., q_info = .TRUE.
REAL (KIND = 8), SAVE :: cut_factor = 1.0d0, distortion_max = 0.0d0, skip = 1.3d0
! parameters
INTEGER       , SAVE :: icut = 26                ! cut at 24 16.5kHz; 25 17.2kHz; 26 17.9kHz; 27 18.6kHz
INTEGER       , SAVE :: mblock_type_param = 20   ! default short block! long = 0, short = 20, mixed = 21
REAL (KIND = 8), SAVE :: ath_min   = -125.0d0     ! offset  for ATH  90.3dB = 2^-15 16bit wav pcm assumed   
REAL (KIND = 8), SAVE :: ath_max   =   0.0d0      ! ceiling for ATH                    (see init_absolute_threshold inpsycho.f90) 
REAL (KIND = 8), SAVE :: switch    =   2.0d0      ! long/short window switching factor (see switch_q in psycho.f90)
REAL (KIND = 8), SAVE :: xms       =   0.8d0      ! NS/MS switching factor             (see mid_side in layer3.f90)
REAL (KIND = 8), SAVE :: xsm       =   1.5d0      ! Short/Mixed switching factor       (see mid_side in layer3.f90)
REAL (KIND = 8), SAVE :: offset    =  40.0d0 ![dB]! offset for masking                 (see psycho in psycho.f90)
REAL (KIND = 8), SAVE :: tempo     =   0.85d0     ! temporal masking parameter         (see psycho in psycho.f90)
REAL (KIND = 8), SAVE :: pm_factor =   1.0d0      ! factor for psychoacoustic moment   (see psycho in psycho.f90) 
REAL (KIND = 8), SAVE :: factor    =   0.4d0      ! distribute bits between 2 granules * n channels by total intensity (see av_bits in layer3.f90) 
REAL (KIND = 8), SAVE :: r0 = 0.33d0, r1 = 0.75d0 ! ISO suggests r0 = 0.33d0, r1 = 0.75d0 (see layer3.90)
!-------------------------------------------------------------------------------------------
END MODULE mod_mpg
