function y=filterchannel(X,h)
N=length(X(:,1));
W=length(X(1,:));

H=fft(h,N);

invH=1./H;
for k=1:W
    y(:,k)=invH.*X(:,k);
end
y;
