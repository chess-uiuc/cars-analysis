function f2 = conv_spk(f,g);

%check the length of the inputs to determine the output vector size
%this should be the length of the largest input vector
if (length(f) < length(g)); M = length(g);else; M = length(f);end

%check length of input vector to treat N = even case by adding a zero to
%the end of the vector
if (mod(length(f),2) == 0); f = [f 0];end
if (mod(length(g),2) == 0); g = [g 0];end

N = abs(length(f)-length(g));
if length(f) < length(g)
    f = [zeros(1,N/2) f zeros(1,N/2)];
else
    g = [zeros(1,N/2) g zeros(1,N/2)];
end

F = fft(f);G = fft(g);
f2 = F.*G;f2 = ifftshift(ifft(f2));

f2 = f2(1:M);

end