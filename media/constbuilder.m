function y = constbuilder(N)
b=ceil(log2(N));
p=zeros(N,b);
for k=0:(N-1)
        u=dec2bin(k);
        w=zeros(1,length(u));
        for i=1:length(u)
            w(i)=str2num(u(i));
        end
        v=[zeros((b-length(w)),1)',w];
        p(k+1,:)=v;
end
n=0:N-1;
    if N<=16
        a=exp(j*2*pi/N.*n');
        M=[p,a];
    else
        a=exp(j*2*pi/16.*(0:15)');
        if N<=24
            b=0.5*exp(j*2*pi/(N-17).*(0:N-17)');
            M=[p,[a;b]];
        else
            b=0.66*exp(j*2*pi/8.*(0:8)'+j*pi/16);
            if N<=32
                c=0.33*exp(j*2*pi/(N-26).*(0:N-26)');
                M=[p,[a;b;c]];
            else
                a=exp(j*2*pi/N.*n');
                M=[p,a];
            end
        end
    end
if N==16
    a=[.354+.354*j;sqrt(2)/2;sqrt(2)/2*j;-.354+.354*j;-sqrt(2)/2*j;.354-.354*j;-.354-.354*j;-sqrt(2)/2];
    b=exp(j*2*pi/8.*(0:7)');
    M=[p,[a;b]];
else
end    
y=matlaberr(M);
           
            