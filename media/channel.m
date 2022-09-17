function y=channel(x,h,aN,sN)

a=conv(x,h);
N=length(a);
w=a+aN;
m=max(abs(x));
for k=1:sN
    w(round(abs(randn(1,1))*(N-1))+1)=m; 
end
y=w;