function y = modpvshift(x, pdetect, ptarget, window, hop, fs, sens)

%y = modpvshift(x, pdetect, ptarget, window, hop, fs)
%
%Input arguments: (all vectors are assumed to be column vectors)
%x - signal to be pitch shifted
%pdetect - vector of detected fund. frequencies in Hz (one for each window)
%ptarget - vector of target fund. frequencies in Hz (one for each window)
%window - number of samples per Hanning window **must be an even number**
%hop - hop size between windows
%fs - sampling frequency in Hz
%
%Note that for the algorithm to perform properly, the windows must have at
%least 75% overlap (i.e. hop <= 25% of window size).  Also note that the
%number of samples per window must be an even number and preferably a power
%of two, though this is not a necessity.

%Compute the amount of shift necessary for each window and convert to
%digital frequency.
pdelta = ptarget - pdetect;
pdelta = pdelta / fs;

%Compute the shift in terms of fft bins using the digital frequency.
binbeta = ptarget ./ pdetect;

%Compute the number of windows required and set up a pointer to the start
%of the first windw.
numwinds = ceil((size(x,1) - window)/hop) + 1;
windstart = 1;

%Initialize variables.
y = zeros(size(x));
lastframeinfl = ones(window / 2 + 1);
lastframepeakphase = 1;


for(windnum = 1:numwinds)
    
    %First fetch the samples to be used in the current window, zeropadding
    %if necessary for the last window.
    if(windnum ~= numwinds)
        windx = x(windstart:windstart + window - 1);
    else
        windx = x(windstart:size(x,1));
        windx(size(windx,1) + 1:window) = 0;
        y(size(x, 1) + 1:windstart + window - 1) = 0;
    end
    
    %Apply the Hanning window function to the samples.
    windx = windx .* hanning(window);
    
    %If the window is unvoiced, leave it alone.
    if(pdetect(windnum) == 0)
        y(windstart:windstart + window - 1) = y(windstart:windstart + window - 1) + windx;
    
    %Increment the window start pointer by the hop size.
    windstart = windstart + hop;
    
    %Otherwise, we need to pitch shift.
    else
        
        
    %Compute the fft of the windowed signal.
    windX = fft(windx);
    
    %Since our signal is purely real, working with all of the fft
    %coefficients is redundant.  The negative frequencies are conjugate
    %symmetric with the positive frequencies, so we need only work with the
    %positive frequencies to do our pitch scaling.  We will reconstruct the
    %negative frequencies as a last step before performing the inverse fft.
    halfwind = (window / 2) + 1;
    windX = windX(1:halfwind);
    
    %Find and record the peaks in the spectrum by bin number.  A peak is
    %defined as a bin whose magnitude is greater than its two nearest
    %neighbors on either side.
    peakindex = 1;
    peak = 0;
    for(bin = 3:(halfwind - 2))
        if((abs(windX(bin)) > abs(windX(bin - 1))) & (abs(windX(bin)) > abs(windX(bin - 2))))
            if((abs(windX(bin)) > abs(windX(bin + 1))) & (abs(windX(bin)) > abs(windX(bin + 2))))
                peak(peakindex, 1) = bin;
                peakindex = peakindex + 1;
            end
        end
    end
    
    %Divide the window into regions of influence from each peak.  Assume
    %the influence of a given peak extends to the midpoint of the distance
    %to the next peak (or to the end of the fft).  (An alternate method
    %would be to assume the turning point is the bin of smallest magnitude
    %between the peaks.)  Define a vector which has a value for
    %every bin and points to the element of peak that exerts influence.
    %For instance, if influence(4) = 3, then bin 4 is under the influence
    %of the peak whose bin number can be found at peak(3).
    %First we deal with the first and last peak, whose influences end at
    %the ends of the fft.
    influence = 0;
    influence(1:peak(1),1) = 1;
    influence(peak(size(peak,1)):halfwind,1) = size(peak,1);
    
    %Now we compute the rest of the influences using the midpoint algorithm
    %described above.
    for(index = 2:size(peak,1))
        diff = peak(index) - peak(index - 1);
        len = ceil((diff + 1) / 2);
        influence(peak(index - 1):peak(index - 1) + len) = index - 1;
        influence(peak(index - 1) + len + 1:peak(index)) = index;
    end
    
    %In order to calculate the amount of shift needed for each peak, we
    %must know the frequency that the peak corresponds to.  We find this by
    %fitting a parabola to the peak bin and its two neighbors.  The vertex
    %will be the actual frequency.  The amount to shift by is the frequency
    %times a constant multiple (the ratio of the target pitch to the
    %detected pitch in this window).
    shift = zeros(size(peak,1),1);
    for(peaknum = 1:size(peak,1))
        peakbin = peak(peaknum);
        bins = [peakbin - 1; peakbin; peakbin + 1]';
        vals = 10*log(abs(windX(bins)));
        a = [bins(1)^2 bins(1) 1; bins(2)^2 bins(2) 1; bins(3)^2 bins(3) 1];
        p = inv(a)*vals;
        vertex = (-p(2) / (2*p(1)));
        shift(peaknum) = vertex * binbeta(windnum) - vertex;
        ndetect(peakbin, windnum) = shift(peaknum);
    end
    
    %Now we must shift the peaks (and their corresponding areas of
    %influence).  
    
    %Since we do not constrain the amount of the shift to be an integer
    %number of bins, we will have to use interpolation to calculate what
    %the new values of the integer bins will be.  For the sake of
    %minimizing algorithmic complexity, we will simply use linear
    %interpolation to implement these fractional shifts.  For a peak to be
    %shifted, we will extract all the values in each region of influence
    %into a two column matrix where each bin number in the first column is
    %paired with its complex value in the second column.  Then we will use
    %linear interpolation to insert it into a new fft vector which includes
    %only the complex values we are computing here (all other values are 0) at
    %the appropriate new indices.  Then we will subtract the values from
    %the old peak regions since they have shifted away.  Additionally, the
    %vector will be paired with a new influence showing how these regions
    %change based on the peak shifting.  Instead of referring to indices in the
    %peak vector, the new influence vector will contain the peak number
    %(n) of the peak it is influenced by.  Once this is done for every shifted
    %peak, we will have two vectors which give information about how to
    %manipulate the coefficients and eventually the phase of the original fft
    %to produce the desired pitch scaling.
    
    %If any of the new bin numbers spill over (go beyond the half window
    %size) into the negative frequencies, we will reflect it back into the
    %positive frequencies and perform a complex conjugation, since the
    %signal is real.
    fftdelta = zeros(halfwind,1);
    newinfl = zeros(halfwind,1);
    for(peaknum = 1:size(peak, 1))
        if(true)
            areaindex = 1;
            areadata = 0;
            
            for(binnum = 1:halfwind)
                if(influence(binnum) == peaknum)
                    areadata(areaindex, 1) = binnum;
                    areadata(areaindex, 2) = windX(binnum);
                    areaindex = areaindex + 1;
                end
            end
            
            %Locate this peak in the last frame and accumulate the phasor
            %from the last frame to this one.
            peakphase(peaknum) = lastframepeakphase(lastframeinfl(peak(peaknum)))*exp(j*shift(peaknum)*2*pi/window*hop);
            
            for(areaindex = 1:size(areadata,1))
                fftdelta(areadata(areaindex,1)) = fftdelta(areadata(areaindex,1)) - areadata(areaindex, 2);
            end
            
            areadata(:,1) = areadata(:,1) + shift(peaknum);
            binoffset = areadata(1,1) - floor(areadata(1,1));
            for(areaindex = 2:size(areadata,1))
                newbin = floor(areadata(areaindex,1));
                newbinval = binoffset * areadata(areaindex - 1, 2);
                newbinval = newbinval + (1 - binoffset) * areadata(areaindex, 2);
                if(newbin > halfwind)
                    newbin = 2 * halfwind - newbin;
                    newbinval = conj(newbinval);
                elseif(newbin <= 0)
                    newbin = abs(newbin) + 2;
                    newbinval = conj(newbinval);    
                end
                fftdelta(newbin) = fftdelta(newbin) + newbinval;
                phasor(newbin) = peakphase(peaknum);
            end
            
            newbin = floor(areadata(1,1));
            if(newbin > halfwind)
                newbin = 2 * halfwind - newbin;
                fftdelta(newbin) = conj(fftdelta(newbin) + (1 - binoffset) * areadata(1,2));
            elseif(newbin <= 0)
                newbin = abs(newbin) + 2
                fftdelta(newbin) = conj(fftdelta(newbin) + (1 - binoffset) * areadata(1,2));
            else
                fftdelta(newbin) = fftdelta(newbin) + (1 - binoffset) * areadata(1,2);
                phasor(newbin) = peakphase(peaknum);
            end
            
            newbin = ceil(areadata(size(areadata,1),1));
            if(newbin > halfwind)
                    newbin = 2 * halfwind - newbin;
                    fftdelta(newbin) = conj(fftdelta(newbin) + binoffset * areadata(size(areadata,1),2));
            elseif(newbin <= 0)
                newbin = abs(newbin) + 2;
                fftdelta(newbin) = conj(fftdelta(newbin) + binoffset * areadata(size(areadata,1),2));
            else
                fftdelta(newbin) = fftdelta(newbin) + binoffset * areadata(size(areadata,1),2);
                phasor(newbin) = peakphase(peaknum);
            end
        end
    end
   
    %Now that we have a vector specifying the amount of change needed
    %for each bin, we can go through and add the value from fftdelta.
    windX = windX + fftdelta;
    
    %Now we must adjust the phases of the shifted peak regions specified in
    %the newinfl vector to ensure inter-frame phase continuity.
    for(binnum = 1:halfwind)
        windX(binnum) = windX(binnum)*phasor(binnum);
    end
    %Finally, we must re-create the truncated fft coefficients by using the
    %knowledge of the conjugate symmetries that must exist because the
    %signal is real and then take the inverse fft to take our modified
    %signal back to the time domain.  Then we will reinsert this window
    %into the output signal.
    for(index = (halfwind + 1):window)
        windX(index) = conj(windX(2 * halfwind - index));
    end
    windstart + window - 1
    y(windstart:windstart + window - 1) = y(windstart:windstart + window - 1) + real(ifft(windX));
    
    %Increment the window start pointer by the hop size.
    windstart = windstart + hop;
    
    %Save the influence and phase information for use in the next frame.
    lastframeinfl = influence;
    lastframepeakphase = peakphase;
    end
end