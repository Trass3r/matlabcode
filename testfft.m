%function [ output_args ] = testfft( input_args )
%TESTFFT Summary of this function goes here
%   Detailed explanation goes here

close all;

% number of subplots
numRows = 3;
numCols = 1;

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 2000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
%y = x + 2*randn(size(t));     % Sinusoids plus noise
y = x + 0.1*randn(size(t));     % Sinusoids plus noise

subplot(numRows, numCols, 1);
%plot(Fs*t(1:N), y(1:N), 'o', Fs*t(1:N), 0.7*sin(2*pi*50*t)-5, Fs*t(1:N), x-5);
%plot(Fs*t, x, Fs*t(1:N), 0.7*sin(2*pi*50*t), Fs*t(1:N), sin(2*pi*120*t));
plot(Fs*t(1:50), y(1:50))
%legend x c1 c2
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')

pause(0.001); % jframe nastiness
jframe = get(handle(gcf), 'JavaFrame');
jframe.setMaximized(true); % maximize figure


% converts a finite list of equally spaced samples of a function into the list of coefficients of a
% finite combination of complex sinusoids, ordered by their frequencies
% but: it's the periodic extension of the original data!


NFFT = 2 * 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y, NFFT) / L;

% frequencies are integer multiples of the fundamental frequency
% its period is the length of the sampling interval
% the combination of sinusoids obtained through the DFT is therefore periodic with that same period
f = Fs / 2 * linspace(0, 1, NFFT/2); % 0, ... 512 entries ..., 1
Y2 = 2 * abs(Y(1:NFFT/2));

% plot single-sided amplitude spectrum
subplot(numRows, numCols, 2);
plot(f, Y2);
title('Single-Sided Amplitude Spectrum of y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');

line([50 120; 50 120], [0 0; 1.1 1.1], 'Color', 'r' ,'LineStyle', '-');

dBm = 10 * log10(1000 * Y2.^2 / 50);
subplot(numRows, numCols, 3);
plot(f, dBm);
xlabel('Frequency (Hz)');
ylabel('dBm');

%plot(f, angle(Y(1:NFFT/2))/pi*180);
%xlabel('Frequency (Hz)');
%ylabel('Phase angle (degrees)');
