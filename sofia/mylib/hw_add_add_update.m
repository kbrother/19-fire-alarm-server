% vectorized.
function [l, b] = hw_add_add_update(y, l, b, f)

alpha = f(1,:);
beta = f(2,:);
alphac = 1 - alpha;
betac = 1 - beta;
y_alpha = alpha .* y;

len = size(l, 1);
n = size(y, 1);

prev_l = l(:,:);
prev_b = b(:,:);

l(1,:) = y_alpha(1,:) + alphac .* (prev_l(len,:) + prev_b(len,:));
b(1,:) = beta .* (l(1,:) - prev_l(len,:)) + betac .* prev_b(len,:);
for t=2:n
    l(t,:) = y_alpha(t,:) + alphac .* (l(t-1,:) + b(t-1,:));
    b(t,:) = beta .* (l(t,:) - l(t-1,:)) + betac .* b(t-1,:);
end
end
