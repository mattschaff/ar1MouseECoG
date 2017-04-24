%% TEST CASE 0: Does it work% 
%make data
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
wave_array(1).theta = ones(1,5000).*0;
wave_array(1).temp_freq = ones(1,5000).*100;
wave_array(2).theta = ones(1,5000).*pi;
wave_array(2).temp_freq = ones(1,5000).*200;
%wave_array(3).temp_freq = ones(1,5000).*50;
%wave_array(3).theta = ones(1,5000).*pi;
%initialize grid
x = -1:0.01:1;
y = 0;
[X, Y] = meshgrid(x, y);
times = (1:5000)*.001;

data1 = populate_wave(wave_array(1), X, Y, times);
data3 = populate_wave(wave_array(3), X, Y, times);
data2 = populate_wave(wave_array(2), X, Y, times);
combined_data = data1 + data2;% + data3;
combined_data_rand = combined_data + normrnd(0,.3,size(data1));

%% FFT
% this function is in the EEGPLOT toolbox
%addpath(genpath([pwd '/powerpectra/eeglab14_0_0b']));
% third input is srate
spectopo(squeeze(combined_data_rand), 0, 1000)

% wavelets
% second input (c) is standard deviation of the temporal Gaussian window times
% pi; This defines the window length. There is a trade off between window
% length and frequency resolution. c = 0 yields the original time-courses
% (a zero width time window); c -->inf yields the Fourier transform. This
% parameter is important, and encompasses a lot of the failure cases
% associated with wavelets
%
% High freq cut off:    srate/2 (Nyquist frequency)
% Low cut off:          1/window size

% define constants
freq = 10:10:400;
srate = 1000;
windows = 100;
peaks = zeros(1,numel(windows));

% visualize peaks at pi*frequency. The change in power is due to how the
% wavelet cuts off parts of the window based on the window size
for i = 1:numel(windows)
    c = windows(i)*pi;
    
    % 5th input is for baseline: no baseline correction here
    [wvlt_amp, wvlt_phase] = morletwave(freq,c,squeeze(combined_data_rand),srate,0);
    
    %peaks(i) = max(wvlt_amp(:,:,1));
    
    clf;
    imagesc(1:5000, freq, squeeze(wvlt_amp(:,:,1)))
    colorbar
    caxis([0, 20])
    pause(.001)
end

%clf;
%plot(windows, peaks)

% phases
clf
subplot(2,1,1)
imagesc(1:21, freq, squeeze(wvlt_phase(:,2500,:)))
colorbar
subplot(2,1,2)
imagesc(1:5000, freq, squeeze(wvlt_amp(:,:,1)))
colorbar

max_amp = max(wvlt_amp, [], 2);
max_amp = repmat(max_amp, [1,5000,1]);
better_phase = wvlt_phase.*max_amp;
imagesc(1:21, freq, squeeze(better_phase(:,2500,:)))
colorbar
% 
% % MOVIE
% % phases
% movie_times = 4000:numel(times);
% figure(1)
% for i = movie_times
%     clf
%     imagesc(1:21, freq, squeeze(wvlt_phase(:,i,:)))
%     title(num2str(i))
%     colorbar
%     caxis([-3 3]);
%     pause(.001)
% end

%% AR1

% ar1
%define constants
win = numel(x)*3; % twice the number of electrodes
n_tp = size(combined_data_rand, 3);
data = reshape(combined_data_rand, [] ,n_tp);
n_eig = numel(x);

% initialize A - elec x elec x time
A = zeros(numel(x), numel(x), (n_tp - win));
% initialize eigenvalues and vectors
eig_vect = zeros(numel(x), n_eig, (n_tp - win));
eig_val = zeros(n_eig, (n_tp - win));

parfor i = 1:n_tp - win
    A_curr = zeros(numel(x), numel(x));
    fprintf('\n...%d', i)
    % get time chunk
    data_chunk = data(:,i:i+win);
    
    % demean
    data_chunk = data_chunk - mean(data_chunk, 2);
   % y = data_chunk(:,2:end);
   % x_hat = data_chunk(:, 1:end-1);
    
    % get A
    [w,A_curr,C] = arfit(data_chunk',1,1,'sbc'); % data_n should be a time chunk;
    %for j = 1:numel(x)
    %    A_curr(j,:) = y(j,:)'\x_hat';
    %end

    A(:,:,i) = A_curr;
    [vect, val] = eig(A_curr);
    eig_val(:,i) = diag(val);
    eig_vect(:,:,i) = vect;
end
 % plot the vector from the largest eignevalues

 %
[~, idx] = sort(abs(eig_val(:,20)));
large_vect1 = angle(eig_vect(:,(idx == 1),20));
large_vect2 = angle(eig_vect(:,(idx == 3),20));

%plot
figure
plot(unwrap(large_vect1))
hold on
plot(unwrap(large_vect2))
xlabel('Time (samples)')
ylabel('Phase')
xlim([1 numel(x)])

%% Test Case 1: Similar Frequency
% 2 waves with some frequency (temporal & spatial)
% check efficacy as lim(d_freq) --> 0
% randomize other parameters (including wave type)
 srate = 1000;
temp_freqs = 2;
base_freq = 200; %in Hz
delta_freq = 10; % how much to change by
k=1;

%make data
wave_array = struct();
for i=1:2
    wave_array(i).type = 'plane';
    wave_array(i).y_center = ones(1,5000);
    wave_array(i).x_center = ones(1,5000);
    wave_array(i).theta = ones(1,5000);
    wave_array(i).temp_freq = ones(1,5000);
    wave_array(i).spatial_freq = ones(1,5000)*5;
    wave_array(i).amplitude = ones(1,5000);
    wave_array(i).timesteps = [1:5000]; %in s
end

% add changing freq over time
wave_array(1).temp_freq = base_freq - delta_freq.*(sin((temp_freqs(k).*(2*pi))./srate*(1:5000)));

%initialize grid
x = -1:0.01:1;
y = 0;
[X, Y] = meshgrid(x, y);
times = (1:5000)*(1/srate);

%% TEST CASE 2: Noisiness
% N waves with varying noise
% check efficacy as lim(sig/noise) --> 0
% randomize other parameters (including wave type)

%% TEST CASE 3: Amplitude (Canonical Power Spectrum)
% N waves with varying amplitudes
% check efficacy as N --> inf
% randomize other parameters (including wave type)

%% TEST CASE 4: Temporal Resolution (changes in parameters)
% 1 wave with varying parameters over time (esp amplitude)
% check efficacy as freq(parameter change) --> inf
% randomize other parameters (including wave type)

srate = 1000;
temp_freqs = 5;
base_freq = 35; %in Hz
delta_freq = 1; % how much to change by
k=1;

%make data
wave_array = struct();
for i=1
    wave_array(i).type = 'plane';
    wave_array(i).y_center = ones(1,5000);
    wave_array(i).x_center = ones(1,5000);
    wave_array(i).theta = ones(1,5000);
    wave_array(i).temp_freq = ones(1,5000);
    wave_array(i).spatial_freq = ones(1,5000)*5;
    wave_array(i).amplitude = ones(1,5000);
    wave_array(i).timesteps = [1:5000]; %in s
end

% add changing freq over time
wave_array(1).temp_freq = base_freq - delta_freq.*(sin((temp_freqs(k).*(2*pi))./srate*(1:5000)));

%initialize grid
x = -1:0.1:1;
y = 0;
[X, Y] = meshgrid(x, y);
times = (1:5000)*(1/srate);

data = populate_wave(wave_array(1), X, Y, times);
combined_data_rand = data;% + normrnd(0,2,size(data));
combined_data = combined_data_rand; %for ar1
yetAnotherData = data;

% % plot wave
% for i = 1:size(data,3)
%     ext = sprintf('%04d',i);
%     imagesc(data(:,:,i));
%     colorbar
%     caxis([-2 2])
%     title(ext);
%     saveas(gca,['/Users/mschaff/Documents/MATLAB/ar1MouseECoG/images/freq_disc_wave/', ext, '.jpg'], 'jpg')
%     pause(0.01);
% end

% wave done, now use wavelet to get phase
freq = 1:.1:100;
srate = 1000;
window = 3;

% visualize peaks at pi*frequency. The change in power is due to how the
% wavelet cuts off parts of the window based on the window size
c = window*pi;

% 5th input is for baseline: no baseline correction here
[wvlt_amp, ~] = morletwave(freq,c,squeeze(combined_data_rand),srate,0);

%peaks(i) = max(wvlt_amp(:,:,1));

% ar1
%define constants
win = numel(x)*3; % twice the number of electrodes
n_tp = size(combined_data_rand, 3);
data = reshape(combined_data_rand, [] ,n_tp);
n_eig = 4;

% initialize A - elec x elec x time
A = zeros(numel(x), numel(x), (n_tp - win));
% initialize eigenvalues and vectors
eig_vect = zeros(numel(x), n_eig, (n_tp - win));
eig_val = zeros(n_eig, (n_tp - win));

parfor i = 1:n_tp - win
    A_curr = zeros(numel(x), numel(x));
    fprintf('\n...%d', i)
    % get time chunk
    data_chunk = data(:,i:i+win);
    
    % demean
    data_chunk = data_chunk - mean(data_chunk, 2);
   % y = data_chunk(:,2:end);
   % x_hat = data_chunk(:, 1:end-1);
    
    % get A
    [w,A_curr,C] = arfit(data_chunk',1,1,'sbc'); % data_n should be a time chunk;
    %for j = 1:numel(x)
    %    A_curr(j,:) = y(j,:)'\x_hat';
    %end

    A(:,:,i) = A_curr;
    [vect, val] = eigs(A_curr, n_eig);
    eig_val(:,i) = diag(val);
    eig_vect(:,:,i) = vect;
end

% raw data plot

%clf
%plot(squeeze(combined_data_rand(1,1,win+1:n_tp)) + 35, 'r')
%hold on
%plot((1:n_tp - win), wave_array(1).temp_freq(win+1:n_tp), 'b')


rawData = squeeze(yetAnotherData(1,1,:));
[~, peakTimes] = findpeaks(rawData);

frequencyTimes = peakTimes;
frequencies = diff(peakTimes);
frequencies(end+1) = frequencies(end);

%scatter(frequencyTimes, srate./frequencies);

trueFrequencies = interp1(frequencyTimes, srate./frequencies, 1:5000);

sigma = 50;
x = -3*sigma:3*sigma;
gaussianFilter = exp(-x.^2/(2*sigma^2));
gaussianFilter = gaussianFilter ./ sum(gaussianFilter);

trueFrequencies = convn(trueFrequencies, gaussianFilter, 'same');

%clf
%plot(trueFrequencies);

%plot
figure(2)
clf
plot((((1:n_tp - win) - win/2)./srate), (angle(eig_val(1,:))./(2*pi))*srate, 'r')
hold on
plot((1:n_tp - win)./srate, trueFrequencies(win+1:n_tp), 'b')
%plot((1:n_tp - win)./srate, (angle(eig_val(3,:))./(2*pi))*srate, 'b')
legend('Eigenvalue 1', 'True Wave')
%ylabel('Frequency (Hz)'); 
xlabel('Time (s)') 
xlim([win, n_tp]./srate);

figure(1)
clf;
imagesc(1:5000, freq, squeeze(wvlt_amp(:,:,1)))
hold on
plot(1:5000, trueFrequencies, 'w')
colorbar
%caxis([0, 20])
pause(.001)

%% TEST CASE 5: Phase Distribution
% N identical waves with varying positions on a 2D grid
% check efficacy as N --> inf
% also check random phase distribution

%% EXPERIMENTAL CASE: Human ECoG data

data = HUP119.data;
srate = HUP119.srate;
freq = 10.^(0:.01:3);
c = 8*pi;

% 5th input is for baseline: no baseline correction here
[wvlt_amp, wvlt_phase] = morletwave(freq,c,data,srate,0);
wvlt_amp_zscore = zscore(wvlt_amp, [], 2);


% phases
figure(1)
clf;
imagesc(1:size(data,1), freq, squeeze(wvlt_phase(:,2500,:)))
colorbar

%% AR1
%define constants
win = size(data,1)*3; % twice the number of electrodes
n_tp = size(data, 2);
%
n_eig = size(data,1);

% initialize A - elec x elec x time
A = zeros(size(data,1), size(data,1), (n_tp - win));
% initialize eigenvalues and vectors
eig_vect = zeros(size(data,1), n_eig, (n_tp - win));
eig_val = zeros(n_eig, (n_tp - win));

parfor i = 1:n_tp - win
    A_curr = zeros(size(data,1), size(data,1));
    fprintf('\n...%d', i)
    % get time chunk
    data_chunk = data(:,i:i+win);
    
    % demean
    data_chunk = data_chunk - mean(data_chunk, 2);
   % y = data_chunk(:,2:end);
   % x_hat = data_chunk(:, 1:end-1);
    
    % get A
    [w,A_curr,C] = arfit(data_chunk',1,1,'sbc'); % data_n should be a time chunk;
    %for j = 1:numel(x)
    %    A_curr(j,:) = y(j,:)'\x_hat';
    %end

    A(:,:,i) = A_curr;
    [vect, val] = eig(A_curr);
    eig_val(:,i) = diag(val);
    eig_vect(:,:,i) = vect;
end

[res, idx] = sort(eig_val(:,20));
large_vect1 = eig_vect(:,(idx == 1),20);
large_vect2 = eig_vect(:,(idx == 2),20);
large_vect3 = eig_vect(:,(idx == 3),20);

%plot
figure
plot(large_vect1)
hold on
plot(large_vect2)
