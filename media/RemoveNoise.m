% This program filters out noise from a sound clip after applying distortion
function [NoiseFree] = RemoveNoise(soundclip)

N = length(soundclip);

% Make the filter that will remove low frequency noise
filter = ones(1,N);

filter = zeros(1,floor(N/100));
filter = zeros(N-N/100,N);
soundclip = real(ifft(filter.*fft(soundclip)));

% Remove all components below a certain value
soundfreq = fft(soundclip);
for i:N
    if(abs(soundfreq(i))<1.5)
        soundfreq(i) = 0;
    end
end

    