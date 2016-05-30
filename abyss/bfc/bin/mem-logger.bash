#!/bin/bash

# log top memory-consuming process in the background
start_mem_log() {
	mem_log=$1
	echo -e 'PID\tRSS\tVIRT\tCMD' > $mem_log
	while true; do
		ps -eo pid,rss,vsz,cmd --width 100 --sort -vsz | awk 'NR==2' >> $mem_log
		sleep 10
	done &
	mem_logger_pid=$!
}

# stop logging memory usage
stop_mem_log() {
	kill $mem_logger_pid
	wait $mem_logger_pid || true
}
