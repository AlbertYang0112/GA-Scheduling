#!/bin/bash

time=0
pid=0
path_to_flame=""

if [ $# -eq 2 ]; then
    path_to_flame=$1
    time=$2
elif [ $# -eq 3 ]; then
    path_to_flame=$1
    time=$2
    pid=$3
else
    echo "Usage: $0 path-to-flamegraph seconds [pid]"
    exit 1
fi

if [ $pid -gt 0 ]; then
    perf record -a -g -p $pid -o perf.data &
else
    perf record -a -g -o perf.data &
fi

PID=`ps aux| grep "perf record"| grep -v grep| awk '{print $2}'`

if [ -n "$PID" ]; then
    sleep $time
    kill -s INT $PID
fi

sleep 1     # wait until perf exite

perf script -i perf.data &> perf.unfold
perl ${path_to_flame}/stackcollapse-perf.pl perf.unfold &> perf.folded
perl ${path_to_flame}/flamegraph.pl perf.folded >perf.svg

echo "Output : perf.svg"
