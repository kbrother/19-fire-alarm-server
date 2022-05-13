%  ~/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -r "run('test.m');"
clc; 
clear; 
close all;

addpath(genpath('lib'));
addpath(genpath('sofia'));
addpath(genpath('sofia/mylib'));
addpath(genpath('../MATLAB Add-Ons'))

fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 22, 'Password', "KaistFire!123");
