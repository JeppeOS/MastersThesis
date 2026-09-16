@echo off
setlocal enabledelayedexpansion

REM Input and output folders
set INPUT=C:\Users\jeppe\Desktop\Bioinf\MiSeq_8_7\Jeppe-Arch
set OUTPUT=C:\Users\jeppe\Desktop\Bioinf\MiSeq_8_7\cutadapt_arch

if not exist "%OUTPUT%" mkdir "%OUTPUT%"

REM Primer sequences
set FWD=GGGYGCAGCAGKCGMGAA
set REV=GTGCTCCCCCGCCAATTCCT

for %%F in ("%INPUT%\*_R1_001.fastq.gz") do (

    set R1=%%F
    set R2=!R1:_R1_001.fastq.gz=_R2_001.fastq.gz!

    set NAME=%%~nF
    set NAME=!NAME:_R1_001=!

    echo Processing !NAME!

    cutadapt ^
        -j 8 ^
        -g %FWD% ^
        -G %REV% ^
        --discard-untrimmed ^
        -e 0.1 ^
        -O 15 ^
        -o "%OUTPUT%\!NAME!_R1.fastq.gz" ^
        -p "%OUTPUT%\!NAME!_R2.fastq.gz" ^
        "!R1!" "!R2!"
)

echo.
echo Finished!
pause