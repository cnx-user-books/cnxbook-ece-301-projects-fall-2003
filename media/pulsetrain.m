function  [s,p,fs,v] = pulsetrain(L,N,M,wd)
% creates a pulse train with
% inputs--
% L: number of pulses
% M: interpulse period
% N: pulse length
% wd: doppler shift frequency
% outputs--
% s: burst waveform of pulses
% p: one pulse
% fs: spectrum of burst waveform (fftshifted)
% v: velocity of object that sensed

% speed of light constant
c=3.0*10^8;


p=ones(N,1);

% TN: total number of points
TN=M*L;
s=zeros(TN,1);

for n=1:N,
    for l=0:(L-1),
        if (n+l*(M))<TN,
            s(n+l*(M))=1;
        end
    end
end

for n=1:TN,
    s(n)=s(n)*exp(j*wd*n);
end

fs=fftshift(fft(s));
[maximum,ind]=max(abs(fs))

% plot the signal and it's spectrum
figure,
stem(cos(wd.*(1:TN)),'r');
hold on;
stem(real(s));

figure,stem(abs(fs));
% HOW DO YOU CALCULATE VELOCITY ???

fc=1/M
v=wd*c/(4*pi*fc)