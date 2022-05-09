for i in {1..100}
do
	ps -o pid,ppid,pgid,comm,%cpu,%mem  -u kaist | grep 28835 >> memory.txt
	sleep 1
done
