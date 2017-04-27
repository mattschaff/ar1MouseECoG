%% Test Case 3: Similar Frequency
% 2 waves with some frequency (temporal & spatial)
% check efficacy as lim(d_freq) --> 0
% randomize other parameters (including wave type)

srate = 1000;
x = -1:0.1:1; %21 electrodes 
y = 0;

baseFreq = 50; %Hz
freqVec = [baseFreq+.1, baseFreq+0.25, baseFreq+0.5, baseFreq+1, baseFreq+2, baseFreq+5, baseFreq+10, baseFreq+20];

% define constants for wavelet
freq = 40:0.05:80;
srate = 1000;
windows = 20;
peaks = zeros(1,numel(windows));
c = windows*pi;





for j = 1:length(freqVec)

    wave_array = struct();
    for i=1:3
        wave_array(i).type = 'plane';
        wave_array(i).y_center = ones(1,5000);
        wave_array(i).x_center = ones(1,5000);
        wave_array(i).theta = ones(1,5000);
        wave_array(i).temp_freq = ones(1,5000);
        wave_array(i).spatial_freq = ones(1,5000)*5;
        wave_array(i).amplitude = ones(1,5000);
        wave_array(i).timesteps = [1:5000]; %in s
    end
    wave_array(1).theta = ones(1,5000).*pi;
    wave_array(1).temp_freq = ones(1,5000).*baseFreq;
    wave_array(2).theta = ones(1,5000).*(0.25*pi);
    wave_array(2).temp_freq = ones(1,5000).*freqVec(j);

    %wave_array(3).temp_freq = ones(1,5000).*50;
    %wave_array(3).theta = ones(1,5000).*pi;
    %initialize grid

    x = -1:0.1:1; %21 electrodes 
    y = 0;
    [X, Y] = meshgrid(x, y);
    times = (1:5000).*(1/srate);

    data1 = populate_wave(wave_array(1), X, Y, times);
    data2 = populate_wave(wave_array(2), X, Y, times);
    combined_data = data1 + data2;
    combined_data_rand = combined_data + normrnd(0,.3,size(data1));
    
    
    % Wavelet!
    [wvlt_amp, wvlt_phase] = morletwave(freq,c,squeeze(combined_data_rand),srate,0);

    figure(1) %power spectrum
    clf;
    imagesc(1:5000, freq, squeeze(wvlt_amp(:,:,1)))
    title(['Wavelet Frequency Discrimination ',  num2str(baseFreq), ' Hz vs. ', num2str(freqVec(j)), ' Hz'])
    ylabel('Frequency')
    xlabel('Time')
    colorbar
    caxis([0, 20])
    pause(.001)
    
    saveas(gca, ['/Users/adeetiaggarwal/Documents/ar1MouseECoG/images/waveletsFreqDis', num2str(freqVec(j)), '.jpg'], 'jpg')

%     figure(2) %amplified phase with less noise 
%     clf
%     max_amp = max(wvlt_amp, [], 2);
%     max_amp = repmat(max_amp, [1,5000,1]);
%     better_phase = wvlt_phase.*max_amp;
%     imagesc(1:21, freq, squeeze(better_phase(:,2500,:)))
%     colorbar

%     % ar1
%     %define constants
%     win = numel(x)*10; % twice the number of electrodes
%     n_tp = size(combined_data_rand, 3);
%     data = reshape(combined_data_rand, [] ,n_tp);
%     n_eig = 5;
% 
%     % initialize A - elec x elec x time
%     A = zeros(numel(x), numel(x), (n_tp - win));
%     % initialize eigenvalues and vectors
%     eig_vect = zeros(numel(x), n_eig, (n_tp - win));
%     eig_val = zeros(n_eig, (n_tp - win));
% 
%     
%     parfor i = 1:n_tp - win
%         A_curr = zeros(numel(x), numel(x));
%         fprintf('\n...%d', i)
%         % get time chunk
%         data_chunk = data(:,i:i+win);
% 
%         % demean
%         data_chunk = data_chunk - repmat(mean(data_chunk, 2), 1, size(data_chunk, 2));
%        % y = data_chunk(:,2:end);
%        % x_hat = data_chunk(:, 1:end-1);
% 
%         % get A
%         [w,A_curr,C] = arfit(data_chunk',1,1,'sbc'); % data_n should be a time chunk;
%         %for j = 1:numel(x)
%         %    A_curr(j,:) = y(j,:)'\x_hat';
%         %end
% 
%         A(:,:,i) = A_curr;
%         [vect, val] = eigs(A_curr, n_eig);
%         eig_val(:,i) = diag(val);
%         eig_vect(:,:,i) = vect;
%     end
%     
%     large_vect1 = angle(eig_vect(:, 1, 20));
%     large_vect2 = angle(eig_vect(:, 3, 20));
%     eigFreq1 = (angle(eig_val(1,:)).*srate)./(2*pi);
%     eigFreq3 = (angle(eig_val(3,:)).*srate)./(2*pi);
%     av_ab_eig_val1 = mean(abs(eig_val(1,:)));
%     av_ab_eig_val3 = mean(abs(eig_val(3,:)));
%     
%     %plot
%     figure(1)
%     clf
%     plot(unwrap(eigFreq1))
%     hold on
%     plot(unwrap(eigFreq3), 'r')
%     title(['AR1 Frequency Discrimination', num2str(baseFreq), ' Hz vs. ', num2str(freqVec(j)), ' Hz'])
%     xlabel('Time (samples)')
%     ylabel('Frequency')
%     legend(num2str(av_ab_eig_val1), num2str(av_ab_eig_val3))
%     xlim([1 numel(x)])
%     
%     saveas(gca, ['/Users/adeetiaggarwal/Documents/ar1MouseECoG/images/AR1FreqDisTakeTwo', num2str(freqVec(j)), '.jpg'], 'jpg')


end
