function info = method_wrapper(Y, O, Omega, R, forecast_steps, opts)


%% ========================= Initialization ============================ %%
Ysz         = size(Y);
N           = ndims(Y);
ntimes      = Ysz(N);
colons      = repmat({':'}, 1, N-1);

T           = opts.T;
mu          = opts.mu;
phi         = opts.phi;
lambda1     = opts.lambda1;
lambda3     = opts.lambda3;

if isfield(opts, 'maxEpoch'), maxEpoch = opts.maxEpoch;
else, maxEpoch = 500; end

if isfield(opts, 'tol'), tol = opts.tol;
else, tol = 1e-4; end

if isfield(opts, 'init'), init = opts.init;
else, init = 'rand'; end

if isfield(opts, 'Omega_orig'), Omega_orig = tensor(opts.Omega_orig);
else, Omega_orig = tenones(Ysz); end


t_end = Ysz(N)-forecast_steps;
assert(t_end >= T);

Y_corrupt = Y+O;
Y_corrupt = Y_corrupt(colons{:}, 1:t_end);

info.Y = Y(colons{:}, 1:t_end);
info.Y_corrupt = double(Y_corrupt);


%% ============================== Method =============================== %%
[U, W, X_hat, Oout, info_method] = sofia(Y_corrupt, Omega(colons{:}, 1:t_end), ...
    R, T, lambda1, lambda3, mu, phi, false);

info.X_hat = X_hat;

%% ===================================================================== %%

info.nre = zeros(1, t_end);
info.rmse = zeros(1, t_end);

for t=1:t_end
    nre = compute_nre(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*X_hat(colons{:}, t));
    rmse = compute_rmse(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*X_hat(colons{:}, t));
%     nre = compute_nre(T(colons{:}, t), X_hat(colons{:}, t));
%     rmse = compute_rmse(T(colons{:}, t), X_hat(colons{:}, t));
    info.nre(t) = nre;
    info.rmse(t) = rmse;
end

info.time = 1:t_end;
info.elapsed = info_method.elapsed_total;
% info.elapsed_update = info_method.elapsed_online;
info.art = info.elapsed/ntimes;
% info.art_update = info.elapsed_update/ntimes;
info.rae = sum(info.nre(~isinf(info.nre)))/ntimes;
info.U = U;
info.W = W;

% info.scale_time = info_method.scale_time;
% info.scale_elapsed = info_method.scale_elapsed;
% info.scale_iter_elapsed = info_method.scale_iter_elapsed;


if forecast_steps > 0
    info.forecast_time = t_end+1:ntimes;
    L = info_method.L;
    B = info_method.B;
    U{N} = hw_add_add_forecast(L, B, forecast_steps);
    Xpred = full(ktensor(U));
    
    info.forecast_nre = zeros(1, forecast_steps);
    info.forecast_rmse = zeros(1, forecast_steps);
    for t=t_end+1:ntimes
        nre = compute_nre(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*Xpred(colons{:}, t-t_end));
        rmse = compute_rmse(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*Xpred(colons{:}, t-t_end));
        info.forecast_nre(t-t_end) = nre;
        info.forecast_rmse(t-t_end) = rmse;
    end
    info.forecast_rae = sum(info.forecast_nre(~isinf(info.forecast_nre)))/forecast_steps;
end


info.Oout = Oout;
info.internal = info_method;

end