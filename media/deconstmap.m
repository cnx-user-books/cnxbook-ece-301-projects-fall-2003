function y = deconstmap(x,M)
% deconstmap(x) - maps constellation into bitstream (2 bits/symbol)
%   x - bitstream

%M=[0 0 1;0 1 j;1 0 -1;1 1 -j];
m=length(M(1,:))-1;   % # of bits/symbol
N = length(x);
bits=zeros(N*m,1);
for i=0:N-1
    for k=1:length(M(:,1))
        if (x(i+1)==M(k,m+1))
        bits(1+(m*i):m*(i+1))=M(k,1:m);
        break
    else
    end
end
end

y=bits;