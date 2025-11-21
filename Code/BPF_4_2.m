%% 4.2 Part (a) Design L=41 filter and measure response
clear; close all; clc;
% Define parameters
L  = 41;         % filter length
wc = 0.25*pi;   % desired frequency
n  = 0:L-1;       % time index

% Computing the impulse response h[n]
w_hamming = 0.54 - 0.46*cos(2*pi*n/(L-1));    %Hamming window
shift     = n - (L-1)/2;                        % time shift
h         = w_hamming .* cos(wc*shift);        %bandpass impulse response



% Normalize the filter coefficients
h = h ./ max(abs(freqz(h,1,1024)));

% Frequency response ( Magnitude & phase) with freqz
[H,w] = freqz(h, 1, 1024);  % H(e^jw) at 1024 frequency axis (0 to pi)
magH  = abs(H);
phH   = angle(H);

% Plot magnitude and phase

subplot(2, 1, 1);
plot(w/pi, magH);
grid on;
title('Magnitude Response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('|H(e^{j\omega})|');
subplot(2, 1, 2);
plot(w/pi, phH);
grid on;
title('Phase Response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

% Find magnitude and phase at specific frequencies:
omega_list = [0, 0.1*pi, 0.25*pi, 0.4*pi, 0.5*pi, 0.75*pi];

% Evaluate magnitude and phase at specified frequencies
H_specific = freqz(h, 1, omega_list);
magAtOmega = abs(H_specific);
phAtOmega  = angle(H_specific);

% Store the results in a table 
resultsTable = table((omega_list'/pi), magAtOmega', phAtOmega', ...
    'VariableNames', {'omega_over_pi', 'Magnitude', 'Phase_rad'});

disp(resultsTable);

%{Explanation-
For L = 41, the magnitude response shows a clear band centered at ω = 0.25π, indicating that the Hamming BPF successfully passes the desired frequency component. The filter exhibits strong attenuation at the other tested frequencies (0, 0.1π, 0.4π, 0.5π, and 0.75π), with magnitude values near zero, confirming that these lie in the stopband. The Hamming window produces very small stopband ripples (<0.01), consistent with the lab description. The phase response is approximately linear in the passband, as expected for a symmetric FIR filter. Overall, the measured values match the theoretical behavior of a Hamming-windowed bandpass filter.
%}

%% 4.2 Part (b)  Passband width vs filter length (50% level)

wc = 0.25*pi;               % desired center frequency
L_list = [21, 41, 81];      % filter lengths to compare

pb_low   = zeros(size(L_list));   % lower passband edge (rad)
pb_high  = zeros(size(L_list));   % upper passband edge (rad)
pb_width = zeros(size(L_list));   % passband width (rad)

figure;  % new figure for magnitude responses

for i = 1:length(L_list)
    L = L_list(i);
    n = 0:L-1;

    % Designing Hamming bandpass filter for this L
    w_hamming = 0.54 - 0.46*cos(2*pi*n/(L-1));   % Hamming window
    shift     = n - (L-1)/2;                     % center shift
    hL        = w_hamming .* cos(wc*shift);      % BPF impulse response

    % Normalized
    hL = hL ./ max(abs(freqz(hL,1,1024)));

    % Frequency response 
    [HL, wL] = freqz(hL, 1, 1024);
    magHL    = abs(HL);

    % 50% passband 
    peak_mag = max(magHL);
    thresh   = 0.5 * peak_mag;          % 50% of peak
    idx_pass = find(magHL >= thresh);   % indices in passband

    pb_low(i)   = wL(idx_pass(1));      % lower edge (rad)
    pb_high(i)  = wL(idx_pass(end));    % upper edge (rad)
    pb_width(i) = pb_high(i) - pb_low(i);

    % Plot magnitude response for this L
    subplot(3,1,i);
    plot(wL/pi, magHL);
    grid on;
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('|H(e^{j\omega})|');
    title(sprintf('Magnitude Response, L = %d', L));
      hold on;
    yline(thresh, 'r--', '50% level', 'LineWidth', 1.2); 
    plot([pb_low(i) pb_high(i)]/pi, [thresh thresh], 'ro', 'MarkerSize', 6, ...
         'LineWidth', 1.2); 
end

% Store the results in a table 
resultsB = table(L_list.', pb_low.'/pi, pb_high.'/pi, pb_width.'/pi, ...
    'VariableNames', {'L', 'omega_low_over_pi', 'omega_high_over_pi', 'width_over_pi'});

disp('Passband Edges and Width for L = 21, 41, 81:');
disp(resultsB);
