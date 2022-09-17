function [noisytestecho,noisyshifttestecho,rsigmatchlocs,timedelay,range,h] = burst4(L,TW,p,M,n,sampfreq,TD)
% creates chirp train, corresponding match filter impulse response
% inputs--
% L: number of pulses
% TW: time bandwidth product
% p: oversampling factor
% M: interpulse period
% n: noise factor
% sampfreq: sampling frequency to be used 
% TD: The amount the transmitted signal is to be shifted 

% outputs--
% noisytestecho: noisy burst waveform of L lfm chirps of the same TW,p 
% noisyshifttestecho: shifted (i.e. time-delayed) of noisy burst waveform
% rsigmatchlocs: vector containing locations of peaks from output of match filter of recieved signal
% timedelay: difference in location of 1st peaks from outputs of both match filters 
% range:  real world target range (units: meters)
% h: impulse response for match filter



[s,chirptrain,h,y,T,W,Ts] = ctbuild4(TW,p,sampfreq,M,L,n);


N=round(p*T*W);

MaxTimeDelay = (M*L) - (M - N)
MaxRange = (1/2)*(3*(10^8))*(MaxTimeDelay)*(1/(sampfreq))

% assuming a row vector

for f = 1:length(chirptrain),
    
    z = f - TD ;
    if z <= 0,
        shifttestecho(f) = 0;
    else 
        shifttestecho(f) = chirptrain(f - TD);
    
    end

end





% ADD NOISE to simulate received signal

noise =(n.*crandn(1,length(shifttestecho)));

                                   %create a noise signal using 
                                   %    crandn function for 
                                   %    complex white
                                   %    gaussian noise
                                   % SNR -10 dB right now
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
                                         


noisesize=size(noise)              %checks to see how many elements
                                   %    are in noise vector
chirptrainsize=size(chirptrain)    %checks to see how many elements
                                   %    are in the pulsetrain 

figure, plot(real(chirptrain))     %Take a look at the plot of 
                                   %    the noise free
                                   %    chirp pulse train
title('Transmitted chirptrain NOISE free')



echosize=size(shifttestecho)

noisyshifttestecho= noise + shifttestecho;

                                   %Add complex white Gaussian
                                   %    noise to time delayed test 
                                   %    echo
  
figure, plot(real(noisyshifttestecho))  %Take a look at the plot of
                                   %    chirp pulse train with 
                                   %    noise added to it
title('Shifted Time-Delayed Chirptrain with NOISE')


noisytestecho=(noise)'+ chirptrain; %Add complex white Gaussian
                                   %    noise to transmitted signal
                                   
figure, plot(real(noisytestecho))  %Take a look at the plot of
                                   %    chirp pulse train with 
                                   %    noise added to it
                                   
title('transmitted Chirptrain with NOISE')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[timedelay,range,rsigmatchlocs] = rangecalc(noisyshifttestecho,chirptrain,h,L,50,Ts)
