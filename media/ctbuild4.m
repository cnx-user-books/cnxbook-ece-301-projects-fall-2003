function [s,ssent,h,y,T,W,Ts] = ctbuild4(TW,p,sampfreq,M,L,n)
% creates chirp train and noisy chirp train
% inputs--
% TW: Time bandwidth product
% p: oversampling factor
% sampfreq: sampling frequency to be used
% M: Interpulse period
% L: # of pulses
% n: noise factor
%
% outputs--
% s = single chirp 
% ssent = noise-free burst waveform of L lfm chirps of the same TW, p
% h = match filter impulse response
% y = (used for testing of original code) passing s through match filter
% T = time duration of LFM chirp (units: sec)
% W= swept bandwidth of LFM chirp (units: Hz)
% Ts= sampling period (1/fs)


[s,h,y,T,W,Ts] = dchirp2(TW,p,sampfreq);                %64,2,N,(1/sampfreq));

% TN: total length of signal (total number of points)
TN=M*L;

N=round(p*T*W);

imp=zeros(TN-(N-1),1);

% builds up impulse train (imp) - for convolving with signal s to create 
% final signal sent (s)

for i=1:L,
    imp(1+(i-1)*M)=1;
end

% convolves impulse train with signal p to give signal sent

ssent=conv(imp,s);