%% Damping times
%
% This script plots the damping time vs frequency from AR models calculated
% in dynamical_criticality_index.m. In short, these are the equations used.
% lambda is the magnitude of the complex eigen value, nd delta t is the
% sampling rate
%
% lambda =          p*exp(i*phi);
% lambda =          p(cos(phi) + isin(phi))
% freq =            phi/2pi*delta_t
% damping rate =    log(p)/delta_t
%
% For moreinformation on calculation of
% individual damping times and frequencies, see Solovey et al 2012 (monkey
% paper), of the DCI script. This script is most just for plotting.

% @author JStiso March 2017

%% Global Varibles

top_dir = '/Users/tnl/Desktop/C/data/eeg/';
subj = 'HUP121_i';
order = 1; % order od model to use
cond = subj(end); % emergence or inductance

data_dir = [top_dir, subj, '/processed_data/dci/']; % where data is
% loada window size
load([data_dir, 'grid_win']);
load([data_dir, 'srate']); % what the data is sampled at
load([data_dir, 'elecs']);
elecs = elec_info.good;
load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/AR_mod']);
save_dir = [data_dir, 'tw', num2str(win), '/Analysis/'];

%% other variables
ref_num = 59792; % number of time bins to use as a conscious reference distribution: should equal 5 seconds

%% Get data into workable structure
% namely, sorted vector of damping times, and matched sorting for frequency

% initialize
dampr_con = zeros(ref_num, numel(elecs));
freq_con = zeros(size(dampr_con));

% move to matrix for conscious and anesthetized
% conscious
for i = 1:ref_num
    ev = eig_val(:,i);
    a = real(ev); % real part
    b = imag(ev); % imaginary part
    p = abs(ev); % magnitude
    phi = atan2(b,a);
    freq_curr = (phi./(2*pi)).*srate.*10;
    dampr_curr = log2(p).*(srate);
    dampr_con(i,:) = dampr_curr;%AR_mod(i).dampr;
    freq_con(i,:) = freq_curr;%AR_mod(i).freq;
end
% reshape into vectors
dampr_v_con = reshape(dampr_con, 1, []);
freq_v_con = reshape(freq_con, 1, []);

% get rid of negative freq
idx_con = find(freq_v_con > 0);
freq_v_con = freq_v_con(idx_con);
dampr_v_con = dampr_v_con(idx_con);
% get rid of outliers (eigenmode = 0)
idx_con = find(dampr_v_con > -50000);
freq_v_con = freq_v_con(idx_con);
dampr_v_con = dampr_v_con(idx_con);

% plot
clf
scatter(dampr_v_con, freq_v_con, 'b*')
%hold on
%scatter(dampr_v_unc, freq_v_unc, 'r*')
xlabel('Damping Rate (1/s)')
ylabel('Frequency (Hz)')
legend('Conscious', 'Unconscious')
ylim([0, 1024])
xlim([-1000, 0])

%  %histogram
%  clf
%  [N,C] = hist3([freq_v_con', dampr_v_con'], [100 10]);
%  imagesc(dampr_v_con, freq_v_con, N)
%  colorbar
%  colormap('hot')
 