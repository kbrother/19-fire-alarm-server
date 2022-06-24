%  ~/MATLAB/R2022a/bin/matlab -nodisplay -nodesktop -r "run('test.m');"
addpath(genpath('sofia/mylib'));
W_init = ones(10 , 10);
options = optimoptions('fmincon');
hw_add_add_fit(W_init);		