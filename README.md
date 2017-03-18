# uzura3
mpeg audio layer 3 encoder written in Fortran 90/95

Never touched since around 2003~4. 

#
Mpeg-I audio layer-III encoder

Fortran95 の命令を使っています。

全く無知の状態から作り始めたので、設計されていません。
心理音響解析もフェイクを入れるところから始めたのでヤバめです。
原理理解のための習作なので一応 VBR なども入っています。
恥ずかしいけど出しますｗ


 Usage : uzura3 -options file_name
       : file_name.wav -> file_name.mp3
 Option:-b  1..14   bitrate for CBR mode                 (default 9 : 128kbps)
        -crc        CRC16 error protection on            (default off)
        -c          copyright flag on                    (default off)
        -o          original  flag on                    (default off)
        -cut 1..32  band cut-off : place after -b option (default 26: 17.9kHz)
        -cuth       cut band 21 (l) or 12 (s/m)          (default off)
        -v          VBR mode  (ns, icut = 32)            (default off)
        -rio500     avoid RIO500 VBR skip bug            (default off)
        -l          long-block-only                      (default off)
        -s          short-mode for short-block           (default off)
        -m          mixed-mode for short-block           (default off)
        -sm         short & mixed-mode for short-block   (default on )
        -xsm xx     short / mixed switching parameter    (default 1.5)
        -switch xx  long/short switching parameter       (default 2.0)
        -skip   xx  speeds up outer loop                 (default 1.3)
        -ms/-ns     stereo mode (MS/normal)              (default MS)
        -xms xx     MS/NS      switching parameter       (default 0.5)
        -nomask     masking off                          (default on)
        -ath_min xx minimum of ATH  [ dB ]               (default -125)
        -ath_max xx ceiling of ATH  [ dB ]               (default  0.0)
        -offset  xx offset for mask [ dB ]               (default 40.0)
        -tempo   xx temporal masking factor              (default 0.85)
        -factor  xx bit distribution among (gr, ch)      (default 0.4)
        -noalias    anti-alias for mixed-block off       (default on)
        -debug      print debug info                     (default on)
        -about      about UZURA3

 Example CBR 128kbps CRC on : uzura3 -crc   file_name
 Example VBR normal stereo  : uzura3 -v -ns file_name



 total frames 7104
 Normal End.
 ======== info ===============================================================
 Block type:long  11117:short    945:mixed    222:type1    962:type3    962
 MS/NS select       391  12849 [  9740  3109] :Long/Short switch   1174    952
..............................................................................
 Average scale factor (long 0-20)   4278
  0.17  0.21  0.21  0.29  0.32  0.36  0.39  0.53  0.54  0.63
  0.70  0.73  0.76  0.78  0.89  0.86  0.83  0.81  0.76  0.50  0.00
 Average scale factor (short 0-11)    583
  1.01  1.05  1.05  1.23  1.24  1.25  1.25  1.23  1.22  1.20  1.12  0.01
 Selected Huffman table (table 0-31)
    7321   1302   5065  10304      0   3157   4372   1003   2862   2517
     398   1281   1678   1374      0   5173    199    316    556    748
    1180    305     27      0   4200   4038   5918   7949   7652   3396
     940     17
 Selected count1 table A, B  26285   2131
..............................................................................
 pre-emphasis      0 :scalefactor_scale    523 :subblock_gain    136
..............................................................................
 bit:    1    2    3    4    5    6    7    8    9   10   11   12   13   14
 (%):  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0100.0  0.0  0.0  0.0  0.0  0.0
 average bit rate (kbps) 128.00      :maximum distortion    0.00013
 ======== info ===============================================================