%% 
% Speech File

analyze("speech8khz.wav");
%% 
% Audio File

analyze("music16khz.wav");
%%
function sound(audio, Fs)
% Sound function, that plays all sound files one-by-one.
y = audioplayer(audio,Fs);
playblocking(y);         
end
%% 
% 

function analyze(filename)
% Analyze function analyzes the sound, as required in the experiment
disp(filename)

% Reading the input audio .wav file
[audio, SamplingRate] = audioread(filename);
disp (["Sampling Rate: ", num2str(SamplingRate)])

% Playing the Audio file
% sound(audio, SamplingRate);

%% Part 1
disp("Part 1");

% Downsampling the audio file
audioDS = downsample(audio, 2);
% sound(audioDS, SamplingRate/2);

% Plotting figure
figure;
[X, w] = freqz(audio, 1, SamplingRate);
subplot(2, 1, 1);
plot(w, abs(X));
title("Original Signal in Frequency Domain");
ylabel("$|X(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");

[XDS, w] = freqz(audioDS, 1, SamplingRate);
subplot(2, 1, 2);
plot(w, abs(XDS));
title("Downsampled Signal in Frequency Domain");
ylabel("$|X_{d}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");

% Part 2
disp("Part 2");

% Designing the low-pass filter
h1 = fdesign.lowpass('N,Fp,Fst', 50, 0.45, 0.55);
filter1 = design(h1, 'equiripple');

% Plotting the filter response. The plots have been commented out for
% clarity

%{
figure
disp("Plotting the filter response:")
freqz(filter1.Numerator, 1, SamplingRate)
%}


% Antialiasing, and then downsampling
out = filter(filter1, audio);
outDS = downsample(out, 2);
[YDS, w] = freqz(outDS, 1, SamplingRate);

% Plotting the magnitude spectrum
figure;
hold on
plot(w, abs(YDS));
title("Output after AA filter and downsampling")
ylabel("$|Y_{d}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");
    
% sound(outDS, SamplingRate/2);
hold off

% Part 3
disp("Part 3");
h2 = fdesign.lowpass('N,Fp,Fst', 50, 0.22, 0.28);
filter2 = design(h2, 'equiripple');
% Plotting filter response
%{
figure
disp("Plotting the filter response:")
freqz(filter2.Numerator, 1, SamplingRate)
%}

% Part 4
disp("Part 4");

% Upsampling and passing through filter and downsampling the audio signal
xup = upsample(audio, 3);
yup = filter(filter2, xup);
ydown = downsample(yup, 4);

% Plotting the figure
figure;
subplot(2, 1, 1);
[XUP, w] = freqz(xup, 1, SamplingRate);
plot(w, abs(XUP));
title("Upsampled Signal in Frequency Domain");
ylabel("$|X_{u}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");

subplot(2, 1, 2);
[YD, w] = freqz(ydown, 1, SamplingRate);
plot(w, abs(YD));
title("Final Signal in Frequency Domain");
ylabel("$|Y_{d}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");

% sound(ydown, SamplingRate * 3/4);

% Part 5
disp("Part 5");

% Plotting interpolation filter
figure
interpolationFilter = intfilt(3, 8, 1);
% freqz(interpolationFilter, 1, SamplingRate)

% Passing through filter, downsampling, and then upsampling
y = filter(filter2, audio);
yd = downsample(y, 4);
yu = upsample(yd, 3);

% Plotting the above signal
[YU, w] = freqz(yu, 1, SamplingRate);
figure
plot(w, abs(YU));
title("Output after Upsampling by 3 in Frequency Domain");
ylabel("$|Y_{d}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");
% sound(yu, SamplingRate * 3/4);

yy = filter(interpolationFilter, 1, yu);
figure
[YY, w] = freqz(yy, 1, SamplingRate);
plot(w, abs(YY));
title("Output after interpolation in Frequency Domain, with image rejection");
ylabel("$|Y_{interpol}(s)|$", 'Interpreter','latex');
xlabel("Frequency (rad)");

% sound(yy, SamplingRate * 3/4);

end