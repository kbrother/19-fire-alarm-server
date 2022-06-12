function [U,X_hat,O,info] = sofia_init(Y,Omega,R,lambda1,lambda3)


%% Default parameters
maxEpoch = 300;
tol = 1e-3;
verbose = 1;
maxIters = 100;
fitchangetol = 1e-3;
als_printitn = 0;

%fprintf('\nsofia_init: start!\n');


%% Parameter setting
Y       = double(Y);
Omega   = logical(double(Omega));
Ysz     = size(Y);
N       = ndims(Y);

U_init = cell(N,1);
for n = 1:N
    U_init{n} = rand(Ysz(n),R);
end


%% Set up for iteration
U = U_init;
O = zeros(Ysz);
X = zeros(Ysz);
normX = 0;
lambda3_init = lambda3;


%% Main Loop
t_begin = tic();
for epoch = 1:maxEpoch
    Xpre = X;
    normXpre = normX;
    
    Yo = Y-O;
    [U,X] = sofia_als(Yo, Omega, R, lambda1, U, maxIters, fitchangetol, als_printitn);

    O = thres_soft(Y-X,lambda3);
    lambda3 = lambda3 * 0.85;
    if lambda3 < lambda3_init/100
        lambda3 = lambda3_init/100;
    end
    
    if epoch > 1
        normresidual = fnorm(Xpre-X);
        relative_change = normresidual/normXpre;
    end
    normX = fnorm(X);
    
    if epoch > 1 && relative_change < tol
        converge = true;
    else
        converge = false;
    end
    
    if verbose >= 1
        if epoch == 1
            %fprintf('  sofia_init: epoch %2d: norm(Xpre) = %e, norm(X) = %e\n', ...
            %    epoch, normXpre, normX);
        else
            %fprintf('  sofia_init: epoch %2d: norm(Xpre) = %e, norm(X) = %e, relative_change = %g\n', ...
            %    epoch, normXpre, normX, relative_change);
        end
    end
    
    if converge
        %fprintf(' sofia_init: epoch %2d: norm(Xpre) = %e, norm(X) = %e, relative_change = %g\n', ...
        %    epoch, normXpre, normX, relative_change);
        %fprintf(' sofia_init: converge!\n');
        break;
    end
end
t_elapsed = toc(t_begin);

info.epoch = epoch;
info.elapsed = t_elapsed;
X_hat = X;

end