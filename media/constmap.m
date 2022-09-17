function y = constmap(x,M)
% constmap(x) - maps constellation onto bitstream (2 bits/symbol)
%   x - bitstream

m=length(M(1,:))-1;
N = length(x);
newN=m*round(N/m);
x=x(1:newN);
for i=1:m
    newvector(:,i)=x(i:m:newN-(m-i));
end

alpha=zeros(newN/m,1);
for i=1:newN/m
    for k=1:length(M(:,1))
        if (newvector(i,:)==M(k,1:m))
        alpha(i)=M(k,m+1);
        break
    else
    end
end
end
y=alpha;