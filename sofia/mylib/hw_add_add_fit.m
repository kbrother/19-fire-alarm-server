% vectorized
function [y_hat, l, b, f] = hw_add_add_fit(y)

y_hat = zeros(size(y));
l = zeros(size(y));
b = zeros(size(y));
k = size(y, 2);
f = zeros(2, k);

for i=1:k
    [y_hat(:,i), l(:,i), b(:,i), f(:,i)] = hw_add_add_fit_internal(y(:,i));    
end

end

function [y_hat, l, b, f] = hw_add_add_fit_internal(y)

assert(isvector(y));

is_y_row = true;
if ~isrow(y)
    is_y_row = false;
    y = y';
end

%% Init Parameters
[l0, b0] = hw_add_add_init_values(y);
init_alpha = 0.5;
init_beta = 0.1 * init_alpha;

x_init = [init_alpha, init_beta, l0, b0];
lower = [0, 0, -Inf, -Inf];
upper = [1, 1, Inf, Inf];
max_fval = Inf;
%disp("init params");

%% Function Handle
funchandle = @(x) hw_add_add_sse_fun(x, y, max_fval);


%% 
eps = 1e-4;
bound = 100;
loc = x_init <= lower;
ub = upper;
ub(ub == Inf) = bound;
x_init(loc) = lower(loc) + eps .* (ub(loc) - lower(loc));

loc = x_init >= upper;
lb = lower;
lb(lb == -Inf) = -bound;
x_init(loc) = upper(loc) - eps .*(upper(loc) - lb(loc));
%disp("function handle");

%% Optimization
%TODO: 
options = optimoptions('fmincon', 'HessianApproximation', 'bfgs', ...
                                  'Display', 'off', ...
                                  'MaxFunctionEvaluations', 500000, ...
                                  'MaxIterations', 1000, ...
                                  'SubproblemAlgorithm', 'cg', ...
                                  'OptimalityTolerance', 1e-10, ...
                                  'FiniteDifferenceStepSize', 1e-9, ...
                                  'StepTolerance', 1e-10, ...
                                  'ConstraintTolerance', 1e-10, ...
                                  'FunctionTolerance', 1e-10);

[x,fval,exitflag,output] = fmincon(funchandle, x_init, [], [], [], [], ...
                                   lower, upper, [], options);
fval;
%disp("optimization");
% exitflag
% output
                               
%% Output
f = [x(1), x(2)];
[y_hat, l, b] = hw_add_add_predict(x, y);

if ~is_y_row
    y_hat = y_hat';
    l = l';
    b = b';
    f = f';
end

end


function [l0, b0] = hw_add_add_init_values(w)
l0 = mean(w);
b0 = mean((w(2:end) - w(1:end-1))/(length(w)-1));
end


function [x, minfval] = hw_add_add_grid_search(x, y, m, N)

factors = linspace(0, 1, N);
% factors(1) = 1e-4;
% factors(end) = 1-1e-4;
minfval = Inf;
for a = 1:N
    for b = 1:N
        for c = 1:N
            x(1) = factors(a);
            x(2) = factors(b);
            x(3) = factors(c);
            fval = hw_add_add_sse_fun(x, y, m, Inf);
            if (fval < minfval)
                minfval = fval;
                alpha = factors(a);
                beta = factors(b);
                gamma = factors(c);
            end
        end
    end
end

x(1) = alpha;
x(2) = beta;
x(3) = gamma;

end


function f = hw_add_add_sse_fun(x, y, max_fval)
len = length(y);
alpha = x(1);
beta = x(2);
alphac = 1 - alpha;
betac = 1 - beta;
y_alpha = alpha .* y;

if alpha * beta == 0
    f = max_fval;
    return;
end
if beta > alpha
    f = max_fval;
    return;
end

l = zeros(1, len);
b = zeros(1, len);

l(1) = x(3);
b(1) = x(4);

for i=2:len
    l(i) = y_alpha(i-1) + alphac * (l(i-1) + b(i-1));
    b(i) = beta * (l(i) - l(i-1)) + betac * b(i-1);
end

f = norm((l + b) - y)^2;
end



function [y_hat, l, b] = hw_add_add_predict(x, y)

len = length(y);
alpha = x(1);
beta = x(2);
alphac = 1 - alpha;
betac = 1 - beta;
y_alpha = alpha .* y;

l = zeros(1, len + 1);
b = zeros(1, len + 1);

l(1) = x(3);
b(1) = x(4);

for i=2:len + 1
    l(i) = y_alpha(i-1) + alphac * (l(i-1) + b(i-1));
    b(i) = beta * (l(i) - l(i-1)) + betac * b(i-1);
end

y_hat = l(1:len) + b(1:len);
l = l(2:end);
b = b(2:end);

end
