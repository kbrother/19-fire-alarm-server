function y_forecast = hw_add_add_forecast(l, b, h)

k = size(l, 2);
y_forecast = zeros(h, k);
for t=1:h
    y_forecast(t,:) = l(end,:) + t .* b(end,:);
end

end

