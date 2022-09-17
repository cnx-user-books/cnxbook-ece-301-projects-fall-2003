function y = a2d(x,b)
% a2d(x) - converts discrete time, continuous amplitude signal to
%          bitstream (b bits per sample)
%  x - original signal
%  b - number bits per sample 

N = length(x);
qinterval = 2*max(abs(x))/((2^b)-1);
if b<=1;
    for i=1:N
        if x(i)<=0
            bits(i)=0;
        else
            bits(i)=1;
        end
    end
else
    q = quant(x-min(x),qinterval);
    bits=zeros(N*b,1);
    
    for k=0:((2^b)-1)
        u=dec2bin(k);
        for i=1:length(u)
            w(i)=str2num(u(i));
        end
        v=[zeros((b-length(w)),1)',w];
        p(k+1,:)=v;
    end
    for i=0:N-1
        for k=0:((2^b)-1)
            if q(i+1)==k*qinterval
                bits(b*i+1:b*i+b)=p(k+1,:);
                break
            else
            end
        end
    end
end
y=bits;
