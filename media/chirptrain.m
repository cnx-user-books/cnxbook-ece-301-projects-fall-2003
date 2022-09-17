function  [s,h,T,W] = chirptrain(TW,p,L,Ts)
% Inputs--
% TW: Time-bandwidth product
% p: oversampling factor, p greater than or equal to 1
% l: length of chirp needed

% If p is not given
%       alternatively can be calcluated from W = fs/p
%                              W = 1/(Ts *p)
%                              Switch variables around
%                              p = 1/(Ts*W)
%
% Ts = 1/ (20*(10 ^6))           % 1 / (p*W)


T = TW * p *Ts

W = TW / T

N = L; %round( T/Ts);   % round(T*W*p)
n = [0:N-1]';

s = exp(j*pi* (W/T) *(n*Ts - T/2).^2 );
h = exp(-j*pi* (W/T) *(n*Ts - T/2).^2 );

% y = conv(s,h);

% figure
% plot(real(s))
% 
% figure
% plot(real(y))
% 
% figure
% plot(abs(y))

%Exercise 2.2 part a)
% >> [sig2,h,y,T,W]=dchirp(64,3);
% 
% Ts =
% 
%   5.0000e-008
% 
% 
% T =
% 
%   9.6000e-006
% 
% 
% W =
% 
%   6.6667e+006 
%
% for n = 1:length(sig2)        %create 700 point signal with chirp
% column = 253+n;               %   chirp in the middle
% blank2(1,column)=sig2(1,n);
% end
% blank2(1,446:700) = 0;


% h=h';                         %create appropiate impulse response
%  for n = 1:length(h)          %   to be used with match filter.
%  column = 253+n;              %Basically, similar to how blank2
%  h2(1,column)=h(1,n);         %   was created, a 700 point signal
% end                           %   with the corresponding h(t) in 
%  h2(1,446:700) = 0;           %   the middle was made

% n = (.2) .* crandn(1,700) ;
% q=std(n);                 
% q/(sqrt(10))                  %
% mean(real(sig2))              %Because the noise desired
                                %   needs to have standard
                                %   deviation equal to sqrt(10)
                                %   of the chirps' amplitude 
                                %   (the average amplitude),
                                %   dividing the standard deviation
                                %   of the noise signal by sqrt(10)
                                %   should make it equal to the average
                                %   amplitude of a chirp. The proof
                                %   is done by visually comparing the
                                %   two numbers and making sure they are
                                %   roughly the same
% clear blank2noise
% blank2noise=blank2+n;
% figure, plot(real(blank2noise))
% y = conv(blank2noise,h);
% figure, plot(real(y))

