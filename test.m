%  ~/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -r "run('test.m');"
clc; 
clear; 
close all;

addpath(genpath('lib'));
addpath(genpath('sofia'));
addpath(genpath('sofia/mylib'));
addpath(genpath('MQTT'))
javaaddpath('MQTT/jar/org.eclipse.paho.client.mqttv3-1.1.0.jar')
javaaddpath('MQTT/mqttasync.jar')

fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");

fireSub = subscribe(fireMQTT, '/CFD/#', 'Callback', @cb);
function cb(topic, msg)
    disp(topic);
end