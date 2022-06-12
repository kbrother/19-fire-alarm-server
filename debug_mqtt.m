%  ~/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -r "run('debug_mqtt.m');"
% ps -A | grep MATLAB | awk '{print $1}' | xargs kill -9 $1
% mosquitto_sub -h 143.248.221.190 -u kaist_fire -P 'KaistFire!123' -t /PRE_T/SFAS81 > sfas81_orig.txt
addpath(genpath('lib'));
addpath(genpath('sofia'));
addpath(genpath('sofia/mylib'));
addpath(genpath('MQTT'))
javaaddpath('MQTT/jar/org.eclipse.paho.client.mqttv3-1.1.0.jar')
javaaddpath('MQTT/mqttasync.jar')

fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");
fireSub = subscribe(fireMQTT, '/CFD/#', 'Callback', @cb);


function cb(topic, msg)
	fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");    
    jsondata = jsondecode(msg);
    num_channel = size(jsondata.cfd, 1);
    output_json.pre = cell(num_channel, 1);

    for i=1:num_channel
        output_json.pre{i} = jsondata.cfd(i);
    end
    
    output_json = jsonencode(output_json);
    publish(fireMQTT, strcat('/PRE_T/', topic), output_json);
    %disp("finish");
end