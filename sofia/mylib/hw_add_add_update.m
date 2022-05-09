% vectorized.
function [l, b] = hw_add_add_update(y, l, b, f)

alpha = f(1,:);
beta = f(2,:);
alphac = 1 - alpha;
betac = 1 - beta;
y_alpha = alpha .* y;

len = size(l, 1);
n = size(y, 1);

for t=len+1:len+n
    l(t,:) = y_alpha(t-len,:) + alphac .* (l(t-1,:) + b(t-1,:));
    b(t,:) = beta .* (l(t,:) - l(t-1,:)) + betac .* b(t-1,:);
end
end
