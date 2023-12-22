f("speech8khz.wav");
f("music16khz.wav");
%%
function sound(audio, Fs)
% Sound function, that plays all sound files one-by-one instead of simultaneously.
y = audioplayer(audio,Fs);
playblocking(y);         
end

function f(filename)
% displaying file name
disp(filename)

% Reading the audio file
[audio, SamplingRate] = audioread(filename);

% Parameters
N = 199; % for this code to work, N should be odd
disp(["the value of N is", num2str(N)])
Fpass = 0.45;
Fstop = 0.55;
Wpass = 1;
Wstop = 1e-3; 

disp("Magnitude response of H0(z)")

% Designing the filter
array = firpm (N, [0 Fpass Fstop 1], [1 1 0 0], [Wpass Wstop]);
Hsuper0 = dfilt.dffir(array);


% Plotting the filter, commented out now for clarity
%{
figure;
freqz(Hsuper0.Numerator, 1, SamplingRate)
%}

% Finding R, defined to be H0(z) * G0(z) = H^0(z)^2
R = conv(Hsuper0.Numerator, Hsuper0.Numerator);

% Setting T (z) = (R (z) - R(-z))/2 and plotting it
for i=1:2:length(R)
    R(i) = 0;
end
T = R;

% Plotting the zp response, commented out now for clarity
%{
figure;
zerophase(T, 1)
%}

% Finding the polyphase components
Hsuper0_0 = Hsuper0.Numerator(1:2:end);
Hsuper0_1 = Hsuper0.Numerator(2:2:end);

% Getting the two outputs of analysis FB
[vsuper0_d, vsuper1_d] = analysis(audio, Hsuper0_0, Hsuper0_1);

% output of synthesis FB for the good case
disp("Good case")
y_good = synthesis(vsuper0_d, vsuper1_d, Hsuper0_1, Hsuper0_0);
figure;
freqz(y_good, 1, SamplingRate);

% output of synthesis FB for the bad case
disp("Bad case")
y_bad = synthesis(vsuper0_d, vsuper1_d, Hsuper0_0, Hsuper0_1);
figure;
freqz(y_bad, 1, SamplingRate);

% Playing the audio files
% sound(audio, SamplingRate);
% sound(y_good, SamplingRate);
% sound(y_bad, SamplingRate);

end

function [vsuper0_d, vsuper1_d] = analysis(x, Hsuper0_0, Hsuper0_1)

% Finding polyphase components of x
x_0 = x (2:2: end);
x_1 = x (1:2: end) ;

tsuper0 = conv(x_0, Hsuper0_0);
tsuper1 = conv(x_1, Hsuper0_1);

% Carrying out the matrix multiplication by DFT matrix
vsuper0_d = tsuper0 + tsuper1;
vsuper1_d = tsuper0 - tsuper1;

end

function [y] = synthesis(vsuper0_d, vsuper1_d, K_0, K_1)

% Carrying out the matrix multiplication by DFT matrix
k_0_input = vsuper0_d + vsuper1_d;
k_1_input = vsuper0_d - vsuper1_d;

y_0 = filter(K_0, 1, k_0_input);
y_1 = filter(K_1, 1, k_1_input);

% Creating the commutator
y = zeros (2 * max(length (y_0), length(y_1)), 1);
for i = 1:length (y_0)
    y (2*i) = y_0 (i) ;
end
for i = 1:length (y_1)
    y (2*i - 1) = y_1(i) ;
end
end