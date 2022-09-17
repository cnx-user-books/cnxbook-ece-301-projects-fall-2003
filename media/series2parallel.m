function y = series2parallel(x,N)

L=length(x);
q=floor(L/N);
newvec=zeros(N,q);
for i=1:q
    newvec(1:N,i)=x((1+(i-1)*N):i*N);
end
y=newvec;