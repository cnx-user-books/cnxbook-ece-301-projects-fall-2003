function  [c,s,h,T,W] = convchirp(TW,p)
% If p is not given
%       alternatively can be calcluated from W = fs/p
%                              W = 1/(Ts *p)
%                              Switch variables around
%                              p = 1/(Ts*W)
%
Ts = 1/ (20*(10 ^6))           % 1 / (p*W)


T = TW * p *Ts;

W = TW / T;

N = round( T/Ts);   % round(T*W*p)
n = [0:N-1]';

s = exp(j*pi* (W/T) *(n*Ts - T/2).^2 );
h = exp(-j*pi* (W/T) *(n*Ts - T/2).^2 );

c=conv(s,h);

% figure
% plot(real(s))
figure
plot(real(c))
% figure
% plot(abs(c))