B= 4; % # of bits
C= 16; % # of points in constellation
S= 256; % # samples of handel
L=16; % length of channel
noise=.0001; % noise amplitude
sN=0; % # TD spikes
N=L; % length of cyclic pad
F=128; % length of blocks (1/2 # of subchannels)
M=constbuilder(C); % constellation matrix

s=abs(randn(L,1));  % random channel
s=s./sum(s);

load handel
y=y(1:S);
x=a2d(y,B);
z=constmap(x,M);
a=series2parallel(z,F);
b=mirror(a);
c=ifft(b); 
d=cyclicpad(c,N);
e=matlaberr(parallel2series(d)); 

aN = randn(length(conv(e,s)),1)*noise;

f=channel(e,s,aN,sN);
f=f(1:length(e));
g=series2parallel(f,2*F+N);
h=decyclicpad(g,N);%
k=fft(h); 
m=filterchannel(k,s);
n=demirror(m);
n(1,:)=a(1,:);
o=parallel2series(n); 
p=approximate(o,M);
q=deconstmap(p,M);
r=d2a(q,B);

% figure(1)
% stem(y/max(abs(y)));title('input speech');
% axis([1 S -1 1]);
% grid on
%     
% figure(2)
% stem(d2a(x,B));ylabel('quantized input');
% axis([1 S -1 1]);
% grid on
% 
% figure(3)
% plot(real(M(:,length(M(1,:)))),imag(M(:,length(M(1,:)))),'.b');
% grid on
% xlabel('REAL');
% ylabel('IMAG');
% axis([-1.5 1.5 -1.5 1.5]);
% 
% figure(4)
% Q=length(d(:,1));
% stem(1:Q-N,matlaberr(d(1:Q-N,1)),'b');ylabel('cylic padded TD input');
% hold
% stem(1:N,matlaberr(d(1:N,1)),'r');
% stem(Q-N+1:Q,matlaberr(d((Q-N+1):Q,1)),'r');
% hold
% 
% 
% figure(5)
% %stem(e);ylabel('series TD input');
% Q=2*(2*F+N);
% stem(1:Q/2,e(1:Q/2));
% hold
% stem(Q/2+1:Q,e(Q/2+1:Q),'g');
% stem(Q/2-L/2+1/2:Q/2+L/2,s,'r');
% hold
% 
% figure(7)
% plot(real(o),imag(o),'*r');
% hold;
% plot(real(M(:,length(M(1,:)))),imag(M(:,length(M(1,:)))),'.b');
% ylabel('Received Constellation Points');
% hold;
% grid;
% axis([-1.5 1.5 -1.5 1.5]);

% figure(6)
% stem(fftshift(abs(fft(s,256))));ylabel('Transfer Function of Channel');
% axis([0 2*F 0 1]);
% 
% figure(8)
% stem(f);ylabel('TD output of channel');
% 
% figure(9)
% stem(h);ylabel('cyclic pad removed');
% 
% figure(10)
% stem(fftshift(abs(1./fft(s,256))));ylabel('transfer function of equalizer');
%   
% figure(11)
% stem(r);ylabel('received/quantized speech'); 
% 
% figure(12)
% stem(1:128,abs(b(1:length(b)/2,1)),'m')
% hold
% stem(129:256,abs(b(length(b)/2+1:length(b),1)),'k')
% hold

ps=(norm(e))^2;
pn=(norm(aN))^2;
snr=10*log10(ps/pn);
biterr=sum(abs(x-q));
