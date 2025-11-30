%Combine 4.1 and 4.2
wc = 0.25*pi;               % desired center frequency
L_list = [21, 41, 81];      % filter lengths to compare

pb_low_rect   = zeros(size(L_list));   % lower passband edge (rad)
pb_high_rect  = zeros(size(L_list));   % upper passband edge (rad)
pb_width_rect = zeros(size(L_list));   % passband width (rad)

pb_low_hamm   = zeros(size(L_list));   % lower passband edge (rad)
pb_high_hamm  = zeros(size(L_list));   % upper passband edge (rad)
pb_width_hamm = zeros(size(L_list));   % passband width (rad)

figure;  % new figure for magnitude responses
titleFont = 16.0;
axesFont = 12.0;
tickFont = 10.0;

for i = 1:length(L_list)
    L = L_list(i);
    n = 0:L-1;

    % Designing bandpass filter for this L
    w_rect = 2/L;   % Rectangular window
    w_hamming = 0.54 - 0.46*cos(2*pi*n/(L-1));   % Hamming window
    shift     = n - (L-1)/2;                     % center shift
    hL_rect        = w_rect .* cos(wc*shift);      % BPF impulse response
    hL_hamm = w_hamming .* cos(wc*shift);      % BPF impulse response


    % Normalized
    hL_rect = hL_rect ./ max(abs(freqz(hL_rect,1,1024)));
    hL_hamm = hL_hamm ./ max(abs(freqz(hL_hamm,1,1024)));

    % Frequency response 
    [HL_rect, wL_rect] = freqz(hL_rect, 1, 1024);
    magHL_rect    = abs(HL_rect);
    
    [HL_hamm, wL_hamm] = freqz(hL_hamm, 1, 1024);
    magHL_hamm    = abs(HL_hamm);

    %% PLOTTING Rectangular
    % 50% passband 
    peak_mag = max(magHL_rect);
    thresh   = 0.5 * peak_mag;          % 50% of peak
    idx_pass_rect = find(magHL_rect >= thresh);   % indices in passband

    pb_low_rect(i)   = wL_rect(idx_pass_rect(1));      % lower edge (rad)
    pb_high_rect(i)  = wL_rect(idx_pass_rect(end));    % upper edge (rad)
    pb_width_rect(i) = pb_high_rect(i) - pb_low_rect(i);

    peak_mag = max(magHL_hamm);
    thresh   = 0.5 * peak_mag;          % 50% of peak
    idx_pass_hamm = find(magHL_hamm >= thresh);   % indices in passband

    pb_low_hamm(i)   = wL_hamm(idx_pass_hamm(1));      % lower edge (rad)
    pb_high_hamm(i)  = wL_hamm(idx_pass_hamm(end));    % upper edge (rad)
    pb_width_hamm(i) = pb_high_hamm(i) - pb_low_hamm(i);

    % Plot magnitude response for this L
    subplot(3,1,i);
    p1 = plot(wL_rect/pi, magHL_rect,'r','LineWidth',1.2);
    hold on;
    p2 = plot(wL_hamm/pi, magHL_hamm,'b','LineWidth',1.2);
    grid on;
    xlabel('Normalized Frequency (\times\pi rad/sample)','FontSize',axesFont);
    ylabel('|H(e^{j\omega})|','FontSize',axesFont);
    title(sprintf('Magnitude Response, L = %d', L),'FontSize',titleFont);
    yline(thresh, 'r--', '50% level', 'LineWidth', 1.2,'Color','k',"FontSize",15.0); 
    %yline(thresh, 'r--', '50% level', 'LineWidth', 1.2); 
    plot([pb_low_rect(i) pb_high_rect(i)]/pi, [thresh thresh], 'ro', 'MarkerSize', 6, ...
         'LineWidth', 1.2);     
    plot([pb_low_hamm(i) pb_high_hamm(i)]/pi, [thresh thresh], 'bo', 'MarkerSize', 6, ...
         'LineWidth', 1.2);     
    % Vertical lines at passband edges
    xline(pb_low_rect(i)/pi,  'r--', 'LineWidth', 1.2);
    xline(pb_high_rect(i)/pi, 'r--', 'LineWidth', 1.2);
    xline(pb_low_hamm(i)/pi,  'b--', 'LineWidth', 1.2);
    xline(pb_high_hamm(i)/pi, 'b--', 'LineWidth', 1.2);
    ax = gca;
    ax.XAxis.FontSize = tickFont;
    ax.YAxis.FontSize = tickFont;
    legend([p1, p2],"Rectangular Window","Hamming Window",'FontSize',15.0,"NumColumns",2);
end

% Store the results in a table 
%resultsB = table(L_list.', pb_low_hamm.'/pi, pb_high_hamm.'/pi, pb_width_hamm.'/pi, ...
%    'VariableNames', {'L', 'omega_low_over_pi', 'omega_high_over_pi', 'width_over_pi'});

resultsB = table(L_list.', pb_low_rect.'/pi, pb_high_rect.'/pi, pb_width_rect.'/pi, ...
            pb_low_hamm.'/pi, pb_high_hamm.'/pi, pb_width_hamm.'/pi, ...
            'VariableNames', {'L','low_rect', 'high_rect', 'width_rect', ... 
            'low_hamm', 'high_hamm', 'width_hamm'});

disp('Passband Edges and Width for L = 21, 41, 81:');
disp(resultsB);
