function [timedelay,range,rsigmatchlocs] = rangecalc(rsig,tsig,h,number,thresh,Ts)
% calculates range, discrete time version of "time delay" of 1st peaks from match filters' outputs  
%
% inputs--
% rsig:  simulated recieved signal
% tsig:  originally transmitted chirp train
% h:    match filter impulse resposne
% number = # of pulses you want to be picked out
% thresh = threshold value used to set where to start looking for peaks 
% Ts = sampling period

% outputs--
%timedelay = difference in location of 1st peaks from outputs of both match filters 
%range =  corresponding real world target range (units: meters)
%rsigmatchlocs = vector containing locations of peaks from output of match filter of recieved signal



% TRY THIS EXAMPLE BELOW FIRST ALLEN...this is for noise having 
%   standard deviation of sqrt(10) multiplied against a 
%   chirp's amplitude
%   NOTE:   
%[pnoise,match] = burst3(5,64,3,16,.2,0);
%[pnoise,match] =burst3(2,64,3,200,.2,0);

rsigmatch=conv(rsig,h);
origsigmatch=conv(tsig,h);

figure, plot(real(rsigmatch))              %Take a look at match filter output 
                                           %    for recieved signal
                                           
title('Matched filter output for recieved signal')                                           

figure, plot(real(origsigmatch))           %Take a look at match filter output 
                                           %    for transmitted signal
                                           
title('Matched filter output for transmitted signal')  
   

[peaksrsigmatch, rsigmatchlocs] = pkpicker( real(rsigmatch), thresh, number, 0 );
[peaksorigsigmatch, origsigmatchlocs] = pkpicker( real(origsigmatch), thresh, number, 0 );



figure, stem(real(rsigmatchlocs))        %Take a look at locations of picked out 
                                         %    peaks for recieved signal
                                           
title('Locations of picked out peaks in recieved signal')                                           





figure, stem(real(origsigmatchlocs))      %Take a look at locations of picked out 
                                          %    peaks for transmitted signal
                                           
title('Location of picked out peaks in transmitted signal')      


                                     
                                           
%compare the location for the first peaks of the match filter outputs
timedelay = rsigmatchlocs(1,1) - origsigmatchlocs(1,1) ;

%Plug resulting time delaty into equation to find corresponding Range 
%   (units = m)
range = (1/2)*(3.00e+008)*(timedelay)*(Ts) ;

