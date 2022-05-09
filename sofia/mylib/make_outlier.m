function O = make_outlier(Tsz, ratio, magnitude)

O = tensor(rand(Tsz) < ratio);
% A = (tenrand(Tsz)-0.5)*2*magnitude;
% O = A .* O;
O = magnitude.*O;
A = (tenrand(Tsz)>0.5)*2-1;
O = A .* O;

% tempSize = [Tsz(1), 1, Tsz(3)];
% 
% O = rand(tempSize) < ratio;
% 
% O = magnitude.*O;
% A = (rand(tempSize)>0.5)*2-1;
% O = A .* O;
% O = repmat(O,1,Tsz(2),1);



end

