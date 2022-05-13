# 19-fire-alarm-server

## How to run?
`
 ~/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -r "run('main.m');" 
`

 ## How to subscribe?
`
mosquitto_sub -h 143.248.221.190 -u kaist_fire -P 'KaistFire!123' -t /PRE_S/#
`
 
 ## How to kill
 `
 ps -A | grep MATLAB | awk '{print $1}' | xargs kill -9 $1
 `
