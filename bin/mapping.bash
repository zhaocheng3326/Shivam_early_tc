#!/bin/bash
DIR=/home/chenzh/My_project/Shivam_early_tc
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin

BDOC=$DIR/big_doc

mkdir -p  $TMPD/X10
conda activate py3.10


#' the first batchcd /home/chenzh/My_project/Shivam_tc/tmp_data/X10
cd $TMPD/X10
/home/chenzh/Cheng_tmpData/Big_tools/cellranger-6.1.1/cellranger  multi --id=SV_ETC --csv=/home/chenzh/My_project/Shivam_early_tc/bin/run.config --localmem 64 --localcores 32


