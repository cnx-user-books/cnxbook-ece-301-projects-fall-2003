function  [s,p,y,fy,v,orig,forig,fyshift] = ptsamp(L,N,M,wd,n)

% TEST: [s,p,y,fy,v,orig,forig,fyshift] = ptsamp(10,128,256,2*pi*(.005),0);

% samples a burst wave so that one sample/pulse
% inputs--
% L: number of pulses
% M: interpulse period
% N: pulse length
% wd: doppler shift frequency
% n: noise factor to add for simulation; (ex: n=.2 for sqrt(10) of chirp amp))
% outputs--
% s: burst waveform of pulses
% p: one pulse
% fs: spectrum of burst waveform (fftshifted
% v: velocity of object that sensed
% close(1,2,3,4)
[s,p,v,orig]=imptrainconv(L,N,M,wd,n);
for n=1:L,
    if ((n-1)*M)<(M*L),        
        x(n)=y(N+(n-1)*M);
    end
end

fx=fft(x);

% [maximum,ind]=max(abs(fx));

% figure,stem(real(x));
% title('sampled peaks from matched filter output of sent signal')
% 
% figure,stem(abs(fx));
% title('fft of sampled peaks');
fyshift=fftshift(fy);
% figure,plot(abs(fyshift));
% title('fft of matched filter output of received signal')
% forig=fftshift(fft(orig));

% [orig_maximum,origmax_ind]=max(abs(forig))
% figure,plot(abs(forig));
% title('fft of original signal');