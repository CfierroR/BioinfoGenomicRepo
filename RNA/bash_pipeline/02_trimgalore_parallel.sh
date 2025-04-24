#!/bin/bash
## Supports N processes at the same time
## For more processes change the MAX_NPROC variable

# $1 --> input_dir
input_dir=$1

NUM=0
QUEUE=""
MAX_NPROC=4

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}

function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

# Main program

find -L ${input_dir} -name "*.fastq.gz" | sed 's/_R1.fastq.gz//g;s/_R2.fastq.gz//g' | sort | uniq | while read sample ;
do
        ########## DEFINE COMMAND #############
	echo "Analyzing sample ${sample}"
        ./02_trimgalore.sh ${input_dir} ${sample} &
        ########### DEFINE COMMAND END ############

        PID=$!
        queue $PID

        while [ $NUM -ge $MAX_NPROC ]; do
                checkqueue
                sleep 0.4
        done
done
wait # wait for all processes to finish before exit

