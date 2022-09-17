function [target, change] = correctfreq(w, allowed, avglength)

% [target, change] = correctfreq(w, allowed, avglength)
%
% Computes the target frequencies and relative change for each detected
% pitch given in w.  w is the vector of pitches calculated for each window.
%
% allowed is a vector of allowed notes in the FIRST octave.
%
% avglength is the length of the smoothing window, set to zero if no
% smoothing is required.
%
% Example scale:
%    55.0000
%    58.2700
%    61.7400
%    65.4100
%    69.3000
%    73.4200
%    77.7800
%    82.4100
%    87.3100
%    92.5000
%    98.0000
%   103.8300

% Generate the possible frequencies
possible = zeros(5 * size(allowed,1),1);
for i = 1 : 5
    
    for j = 1 : size(allowed, 1)
        
        possible((i-1) * size(allowed, 1) + j, 1) = allowed(j,1) * (2^i);
    end

end

% Take the logarithm

index = dsearchn([0; log(possible)], log(w));

possible = [0 ; possible];

for i = 1 : size(w, 1)
    
    target(i, 1) = possible(index(i, 1), 1);
    
end

if(avglength ~= 0)
    
    change = target - w;
    
    change = filter(ones(1,avglength)/avglength,1, change);
    
    target = w + change;
else
    change = 0;
end