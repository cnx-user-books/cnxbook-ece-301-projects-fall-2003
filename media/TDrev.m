function y = TD(x, pdetect, ptarget, window, hop, fs)

%y = TD(x, pdetect, ptarget, window, hop, fs)
%
%Input arguments: (all vectors are assumed to be column vectors)
%x - signal to be pitch shifted
%pdetect - vector of detected fund. frequencies in Hz (one for each window)
%ptarget - vector of target fund. frequencies in Hz (one for each window)
%window - number of samples per Hanning window
%hop - hop size between windows
%fs - sampling frequency in Hz

numwinds = floor((size(x, 1) - window) / hop) + 1;
y = zeros(size(x,1),1);
windstart = 1;

for(windnum = 1:numwinds)
    %Fetch the samples for this window.
    windx = x(windstart: windstart + window - 1);
    windy = zeros(window, 1);
    
    if(pdetect(windnum) ~= 0)
        %Compute the detected and target period in seconds.
        Tdetect = 1 / pdetect(windnum);
        Ttarget = 1 / ptarget(windnum);
        
        %Compute the detected and target period in samples.
        Tsampdetect = floor(Tdetect * fs);
        Tsamptarget = floor(Ttarget * fs);
        
        %Put the first two periods of the window in a vector and apply a
        %Hanning analysis window.
        twoperiods = windx(1:2*Tsampdetect) .* hanning(2*Tsampdetect);
        
        %Calculate the ratio between the desired and actual periods.
        Tratio = pdetect(windnum) / ptarget(windnum);
        
        %Calculate the number of periods that will appear in this window.
        numperiods = floor(window / Tsamptarget);
        
        if(Tratio < 1)
            Tstart = 2;
        else
            Tstart = 1;
        end
        
        %Copy the first period of the signal as is.
        windy(1:Tsampdetect) = windx(1:Tsampdetect) .* hanning(Tsampdetect);
        
        %Overlap the first two periods around each value of new desired
        %period.
        for(Tnum = Tstart:numperiods - 1)
            marker = Tnum*Tsamptarget;
            if(marker + Tsampdetect>length(windy))
                Tsampdetect=length(windy)-marker;
                twoperiods=twoperiods(1:length(windy(marker - Tsampdetect + 1:marker + Tsampdetect)));
            end
            windy(marker - Tsampdetect + 1:marker + Tsampdetect) = windy(marker - Tsampdetect + 1:marker + Tsampdetect) + twoperiods;
        end
        
        %Place the modified window into the output vector.
        windy = windy .* hanning(window) .* Tratio;
        y(windstart: windstart + window - 1) = y(windstart: windstart + window - 1) + windy;
    else
        %This is an unvoiced segment, so we will copy it as is.
        length(y(windstart: windstart + window - 1))
        y(windstart: windstart + window - 1) = y(windstart: windstart + window - 1) + hanning(window) .* windx;
    end
    
    %Increment the window start pointer
    windstart = windstart + hop;
end