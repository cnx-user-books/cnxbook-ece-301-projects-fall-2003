function y = PSOLA(x, pdetect, ptarget, window, hop, fs)

%y = PSOLA(x, pdetect, ptarget, window, hop, fs)
%
%Input arguments: (all vectors are assumed to be column vectors)
%x - signal to be pitch shifted
%pdetect - vector of detected fund. frequencies in Hz (one for each window)
%ptarget - vector of target fund. frequencies in Hz (one for each window)
%window - number of samples per Hanning window **must be an even number**
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
        
        %Calculate the number of periods that will appear in this window.
        numperiods = floor(window / Tsampdetect) - 1;
        
        %Calculate old marker vector, assuming markers are placed every
        %Tsampdetect starting at 1.
        for(i=1:numperiods)
            oldmarker(i,1) = 1 + (i - 1) * Tsampdetect;
        end
        
        %Generate matrix B with each column corresponding to a period in
        %the input window and the containing the Hanning windowed samples
        %around the start of the period.
        B = zeros(Tsampdetect * 2, numperiods);
        newB = zeros(Tsampdetect*2, numperiods);
        C = zeros(Tsampdetect * 2, numperiods);
        
        %Do the first period separately as it is a special case.
        B(1 + Tsampdetect: 2*Tsampdetect, 1) = windx(1:Tsampdetect);
        B(:,1) = B(:,1) .* hanning(Tsampdetect * 2);
        
        %Flip the column.
        for(j = 1:2*Tsampdetect)
            newB(j,1) = B(2*Tsampdetect - j + 1, 1);
        end
        
        %Now do the rest of the periods
        for(i=2:numperiods)
            B(:,i) = windx(oldmarker(i) - Tsampdetect:oldmarker(i) + Tsampdetect - 1) .* hanning(Tsampdetect * 2);
            %Flip the column
            for(j=1:2*Tsampdetect)
                newB(j,i) = B(2*Tsampdetect - j + 1, i);
            end
        end
        
        %Find the optimal path through the matrix B.
        path = pathfinder(newB, 4);
        
        oldmarker = 0;
        newmarker = 0;
        
        %Turn the path vector into pitch markers.
        oldmarker = zeros(numperiods,1);
        for(i = 1:numperiods)
            oldmarker(i) = Tsampdetect - path(i) + (i - 1)*Tsampdetect + 1;
        end
        
        %Copy the first period of the signal as is.
        windy(1:Tsampdetect) = windx(1:Tsampdetect) .* hanning(Tsampdetect);
        
        %Calculate new pitch markers. 
        numnewperiods = floor((window - Tsampdetect - oldmarker(1))/ Tsamptarget);
        for(i=1:numnewperiods)
            newmarker(i,1) = oldmarker(1) + (i - 1) * Tsamptarget;
        end
        
        for(i = 2:numnewperiods)
            %Find the old marker to use for this new marker.
            diff = inf;
            for(num = 2:numperiods)
                if(abs(newmarker(i) - oldmarker(num)) < diff)
                    markertouse = num;
                    diff = abs(newmarker(i) - oldmarker(num));
                end
            end
            thisperiod = zeros(2*Tsampdetect, 1);
            for(t = 1:2*Tsampdetect)  
                if((t - Tsampdetect) + oldmarker(markertouse) >= 1)
                    thisperiod(t, 1) = windx((t - Tsampdetect) + oldmarker(markertouse));
                end
            end
            if(windnum>150 & windnum <= 155)
                C(:,i) = thisperiod;
            end
            thisperiod = thisperiod .* hanning(2*Tsampdetect);
            
            if(newmarker(i) - Tsampdetect >= 1)
                windy(newmarker(i) - Tsampdetect: newmarker(i) + Tsampdetect - 1) = windy(newmarker(i) - Tsampdetect: newmarker(i) + Tsampdetect - 1) + thisperiod;
            end
        end
        %Place the modified window into the output vector.
        windy = windy .* hanning(window);
        y(windstart: windstart + window - 1) = y(windstart: windstart + window - 1) + windy;
    else
        %This is an unvoiced segment, so we will copy it as is.
        y(windstart: windstart + window - 1) = y(windstart: windstart + window - 1) + hanning(window) .* windx;
    end
    
    %Increment the window start pointer
    windstart = windstart + hop;
end

function path = pathfinder(MAT,N)
%
% y = pathfinder(MAT)
%
% This function traces a path from the first to the last columns
% of MAT, one that does not exceed slope == N (N integer >0) when
% assuming that successive rows are separated by one unit, and that
% successive columns are separated by one unit) and has the maximum
% possible cumulative MAT values along the path.  The output
% path y adheres to the sample points of MAT.


% calculate best-path cumulative errors:
[mrows,mcols] = size(MAT);
sf = mean(mean(MAT));
MAT = [-Inf*ones(N,mcols); MAT; -Inf*ones(N,mcols)];
best_paths = zeros(size(MAT));
range = N + (1:mrows);
T = zeros(1+2*N,mrows);
B = zeros(1+2*N,mrows);
R = zeros(1,(1+2*N)*mrows);
for i = -N:N
        B(i+N+1,:) = ones(1,mrows) * sf/sqrt(1+i*i);
        R(mrows*(i+N)+[1:mrows]) = range + i;
end
for col = 2:mcols
        T = reshape(MAT(R,col-1),mrows,1+2*N)';
        [temp1,temp2] = max(T+B);
        MAT(range,col) = MAT(range,col) + temp1';
        best_paths(range,col) = temp2';
end


% trace the optimal path backwards through the cum. error matrix:
best_paths = best_paths - N - 1;
path = zeros(1,mcols);
[total_error,row] = max(MAT(:,mcols));
path(mcols) = row;
for col = mcols:-1:2
        row = row + best_paths(row,col);
        path(col-1) = row;
end

path = path - N;

return