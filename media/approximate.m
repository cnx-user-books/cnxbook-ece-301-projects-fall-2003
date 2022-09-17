function y = approximate(x,M)
% approximate(x,M) - approximates noisy symbol vector x with valid symbols
% in M
%   x - noisy symbol vector
%   M - matrix representing constellation mapping

%M=[0 0 1;0 1 j;1 0 -1;1 1 -j];
W=length(M(1,:));
L=length(M(:,1));
N = length(x);
a = M(:,W);
alpha=zeros(N,1);

for i=1:N
    mindist=100;
    symnum=0;
    for k=1:L
        if abs(a(k)-x(i))<=abs(mindist)
            mindist=abs(a(k)-x(i));
            alpha(i)=a(k);
        else
        end
    end
end        
y=alpha;