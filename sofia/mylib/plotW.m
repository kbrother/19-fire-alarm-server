function plotW(W)
figure;
k = size(W,2);
for i = 1:k
    subplot(k, 1, i);
    plot(W(:, i));
end

