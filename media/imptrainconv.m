function  [sr,s,vcomp,original] = imptrainconv(L,N,M,fd,n)
% creates a impulse train convoluted with a signal p
% inputs--
% L: number of pulses
% M: interpulse period
% N: pulse length
% fd: doppler shift frequency
% outputs--
% sr: signal received- burst waveform
% p: one pulse
% fs: spectrum of burst waveform (fftshifted)
% v: velocity of object that sensed

% speed of light constant
c=3e8; %units:m/s

% sampling frequency:
sampfreq=(19*7)*1e9 %units: Hz

% radar frequency, units: Hz
radar_freq=7e9;
% radar_freq=sampfreq/M
fc=radar_freq;


s=ones(N,1);
% [s,h,T,W] = chirptrain(64,2,N,(1/sampfreq));
h=s;

% TN: total length of signal (total number of points)
TN=M*L;

imp=zeros(TN-(N-1),1);

% builds up impulse train (imp) - for convolving with signal s to create 
% final signal sent (s)
for i=1:L,
    imp(1+(i-1)*M)=1;
end

% convolves impulse train with signal p to give return signal

sr=conv(imp,s);
sr=sr-mean(sr);
% saves the orignal sulse train signal
original=sr;

% ADD NOISE
noise =(n.*crandn(1,length(sr)))'; 
                                   %create a noise signal using 
                                   %    crandn function for 
                                   %    complex white
                                   %    gaussian noise
                                   
 q=std(noise);                     %NOTE: To specify how much
 q/(sqrt(10))                      %    noise, need to set parameter
 mean(real(s))                     %    "n". 
                                   %If want noise to have standard
                                   %    deviation equal to sqrt(10) 
                                   %    multiplied against the chirp's
                                   %    amplitude, need to perfrom 
                                   %    visual proof of equality of . 
                                   %    the two numbers
                                   %Dividing the standard deviation
                                   %   of the noise signal by sqrt(10)
                                   %   should make it equal to the average
                                   %   amplitude of a chirp. The proof
                                   %   is done by visually comparing the
                                   %   two numbers and making sure they are
                                   %   roughly the same
                                         

sr =(noise)+ sr ;                  % Add complex white Gaussian
                                   %    noise to the signal
                                   

% fc=1/M;
% % input delay shift                                  
                                   
% for n=1:TN,
%     sr(n)=sr(n)*exp(j*2*pi*fd*n);
% end
n=1:length(sr);
n=n';
sr=sr.*exp(j*fd*2*pi/sampfreq.*n);
% *2*pi

[Horig,Forig]=freqz(original,1,[],sampfreq);
figure,stem(Forig,abs(Horig));
[Orig_max,Orig_maxind]=max(abs(Horig))
Orig_freq=Forig(Orig_maxind)


[Hrec1,Frec1]=freqz(sr,1,[],sampfreq);
% % % fs=fftshift(fft(sr));
figure,stem(Frec1,abs(Hrec1));
[H_maximum1,Hmax_ind1]=max(abs(Hrec1))
highest_freq1=Frec1(Hmax_ind1)

deltafreq=diff(Orig_freq,highest_freq1)

% [Hrec,Frec]=freqz(sr,1,[(highest_freq1-deltafreq):(2*deltafreq/128):(highest_freq1+deltafreq)],sampfreq);
% figure,stem(Frec,abs(Hrec));
% [H_maximum,Hmax_ind]=max(abs(Hrec))
% highest_freq=Frec(Hmax_ind)



fdcomp=highest_freq-Orig_freq

% plot the signal and it's spectrum
% figure,
% stem(cos(wd.*(1:TN)),'r');
% hold on;
% plot(real(sr));
% title('plot of received signal');

% figure,plot(abs(fs));
% title('plot of fft of received signal');
% % HOW DO YOU CALCULATE VELOCITY ???
vgoal=fd*c/2/fc
vcomp=fdcomp*c/(2*fc)
length=length(sr)


% convolving full received signal with filter matched to one 
% of individual pulses:
y=conv(h,sr);

% figure, plot(real(y));
% title('Matched filter output to Train of pulses (y)');
% % % fy=fft(y);


% figure, plot(abs(fy));
% title('abs of FFT of y');

% z=conv(h,original);
% figure, plot(real(z));
% title('convolving 2 trains (old and new)');