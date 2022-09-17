function y = matlaberr(X)
R=length(X(1,:));
L=length(X(:,1));
for m=1:R
    for n=1:L
        if abs(real(X(n,m)))<.000001
            a=0;
        else
            a=real(X(n,m));
        end
        if abs(imag(X(n,m)))<.000001
            b=0;
        else 
            b=imag(X(n,m));
        end
        X(n,m)=a+j*b;
    end
end
y=X;
