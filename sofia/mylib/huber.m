function out = huber(x)
%k=1;
k=3;
out = x;
p = abs(x)>=k;
out(p)=k.*sign(x(p));
end

