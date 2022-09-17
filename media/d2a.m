function y = d2a(x,b)
% d2a(x) - converts bitstream (b bits per sample) to discrete time,
%           continuous amplitude signal 
%  x - bitstream
%  b - number bits per sample 
x=x';
N = length(x);
L=floor(N/b);
qinterval = 2/(2^b-1);

if b<=1;
    for i=1:N
        if x(i)==0
            sig(i)=-1;
        else
            sig(i)=1;
        end
    end
else
    sig=zeros(L,1);  
    
    for k=0:((2^b)-1)
        u=dec2bin(k);
        w=zeros(1,length(u)); % conversion from interval # to binary
        for i=1:length(u)
            w(i)=str2num(u(i));
        end
        v=[zeros((b-length(w)),1)',w];
        p(k+1,:)=v;     % p is a vector of binary labels for interval #'s
    end
    
    if b~=4    
    
      for i=0:L-1                     % matches value of q with label #
        for k=0:((2^b)-1)       % inserts corresponding binary value into bitstream   
            if (x((1+(i*b)):(b*i+b))== p(k+1,:))
                sig(i+1)=qinterval*k-1;
                break
            else
            end
        end
      end
    
  
    else
        for i=0:L-1
           %%for k=1:b
             if x(1+i*b)==0
                if x(2+i*b)==0
                   if x(3+i*b)==0
                      if x(4+i*b)==0
                          sig(i+1) = -1;
                      else
                          sig(i+1) = qinterval*1-1;end
                   else
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*2-1;
                      else
                          sig(i+1) = qinterval*3-1;end; end
             else
                  if x(3+i*b)==0
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*4-1;
                      else
                          sig(i+1) = qinterval*5-1;end
                   else
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*6-1;
                      else
                          sig(i+1) = qinterval*7-1;end; end; end
                   
             else
                if x(2+i*b)==0
                   if x(3+i*b)==0
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*8-1;
                      else
                          sig(i+1) = qinterval*9-1;end
                   else
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*10-1;
                      else
                          sig(i+1) = qinterval*11-1;end; end
             else
                  if x(3+i*b)==0
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*12-1;
                      else
                          sig(i+1) = qinterval*13-1;end
                   else
                      if x(4+i*b)==0
                          sig(i+1) = qinterval*14-1;
                      else
                          sig(i+1) = qinterval*15-1;end; end; end; end
      end  
      end              
end
y=sig;
