%Frequency Domain Pitch Detection
%f_y = pitch_detec(x, window, hop, xformlength)

function f_y = pitch_detec(x, window, hop, xformlength)

%Windowing input signal
numwinds = ceil((size(x,1) - window)/hop) + 1;
windstart = 1;

h=1;
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
    
    %STFT (Convert from Time Domain to Freq Domain using fft of length 4096)
     f_x = fft(windx, xformlength);
     

    %HPS
    %function f_y = hps(f_x)
    f_x = f_x(1 : size(f_x,1) / 2);
    
    f_x = abs(f_x);
    
    %HPS, PartI: downsampling
    for i = 1:length(f_x)
        f_x2(i,1) = 1;
        f_x3(i,1) = 1;
        f_x4(i,1) = 1;
    %    f_x5(i,1) = 1;
    end

    for i = 1:floor((length(f_x)-1)/2)
        f_x2(i,1) = (f_x(2*i,1) + f_x((2*i)+1,1))/2;
    end

    for i = 1:floor((length(f_x)-2)/3)
        f_x3(i,1) = (f_x(3*i,1) + f_x((3*i)+1,1) + f_x((3*i)+2,1))/3;    
    end
    
    for i = 1:floor((length(f_x)-3)/4)
        f_x4(i,1) = (f_x(4*i,1) + f_x((4*i)+1,1) + f_x((4*i)+2,1) + f_x((4*i)+3,1))/4;
    end
    
    %   for i = 1:floor((length(f_x)-4)/5)
    %      f_x5(i,1) = (f_x(5*i,1) + f_x((5*i)+1,1) + f_x((5*i)+2,1) + f_x((5*i)+3,1) + f_x((5*i)+4,1))/5;
    %   end
    

    %HPS, PartII: calculate product
    f_ym = (1*f_x) .* (1.0*f_x2);%  .* (1*f_x3)  .* f_x4; %.* f_x5;
    
   
    %HPS, PartIII: find max
    f_y1 = max(f_ym);
    
    for c = 1 : size(f_ym)
        if(f_ym(c, 1) == f_y1)
            index = c;
        end
    end
    
    % Convert that to a frequency
    f_y(h) = (index / xformlength) * 44100;
    
    % Do a post-processing LPF
    if(f_y(h) > 600)
        f_y(h) = 0;
    end
    
    %Don't forget to increment the windstart pointer.
    windstart = windstart + hop;
   
%    f_y(h) = f_y1;
    h=h+1;
    f_y = abs(f_y)';

end
   