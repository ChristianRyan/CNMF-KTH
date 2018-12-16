clear all
clc
% OLA Parameters for analysis and synthesis
M = 1024; % window length
R = M/4; % hop size: 75% overlap
w = hamming(M, 'periodic'); % Hamming window
w = w./sum(w(1:R:end)); % Normalized Hamming window

%% Generate large concatenated dataset for training
% Dataset provided in this code is for a speaker specific case. I.e. the
% NNMF is trained for one speaker and the test file is from the same
% speaker but unseen to NNMF yet.

flac_files = dir('train_speaker*');
info = audioinfo(flac_files(1).name); % Info about audio. E.g. sampling rate, etc.
train_V = []; %Dataset concatenated in time domain.
% Concatinating dataset.
for i = 1 : numel(flac_files)
    filename = flac_files(i).name;
    audio_file = audioread(filename);
    train_V = [train_V; audio_file];
end

train_V_freq = stft(train_V, M, R, length(train_V), w); % STFT of concatenated dataset
% Possible improvement: Remove silent pauses (0 in time-domain), as there is no
% point in learning spectral bases for no speech.

%% Find suitable I
% Find the decay of singular values and find an appropriate number of Sv:s
% to describe V.

sv = svd(abs(train_V_freq));
power = sv.^2./(cumsum(sv.^2));
plot(power) % Plot the decay of SVs

phonemes = sv.^2/sum(sv.^2);

figure, plot(cumsum(phonemes)) % Plot the power representation of SVs
c_sum = find(cumsum(phonemes) < 0.995); % Threshold is 99.5%, enough SVs then.
I = c_sum(end) + 1;

%% Learn the dictionary
% Beta-divergence parameters
beta = 0.5;
nu = -Inf;
maxiter = 1000;
showcost = true;

% Beta-convolutional NNMF
for seed = 1 : 100 %Convolutional length of 100, T = 100.
    W0{1,seed} = abs(randn(size(train_V_freq, 1), I) + 1i * randn(size(train_V_freq, 1), I)) + eps; %Random initilzation of each codebook
end
H0 = rand(I, size(train_V_freq, 2)) + eps; % Random Initilzation of H.

[W, ~, ~, cost] = betacnmf(abs(train_V_freq), W0, H0, beta, nu, maxiter, showcost); % Codebook training.

%% Normalize codebook W
%Normalzie each column of each vector in W.
for matrix = 1 : length(W)
    for col = 1 : I
        W{1,matrix}(:,col) = W{1,matrix}(:, col) / sum(W{1,matrix}(:, col));
    end
end

disp('Learning is done')

%% Test item
% Reconstruct a test data file.

% Test file and its STFT counterpart.
[test_file, fs] = audioread('Test_speaker.flac');
N = length(test_file);
test_file_freq = stft(test_file,M,R,N,w);

fc = 4000; % Cutoff frequency.
kc = floor(fc/((fs/2)/(M/2+1))); % Cutoff frequency's index in freq. domain.

test_file_freq_low = test_file_freq(1:kc,:); %Narrowband.

%Obtain lower part of trained W.
for i = 1:length(W)
   W0{1,i} = W{1,i}(1:kc,:) 
end

H0 = rand(I, size(test_file_freq_low, 2)) + eps; %Random init of H.

%Update H without altering the already trained W.
[~, H, ~, cost] = betacnmf_update_H(abs(test_file_freq_low), W0, H0, beta, nu, maxiter, showcost); 

%% Reconstruction
V_est = zeros(size(test_file_freq)); %Pre-allocate memory for reconstruction.
% Convolutional Non-negative matrix factorization:
for i = 1:length(W)
   V_est = V_est + W{1,i}*shr(H,i-1);
end

V_est(1:kc,:) = abs(test_file_freq(1:kc,:)); %Narrowband is already transmitted, thus should remove estimation of this part.

% Plot results between true original signal and its CNMF reconstruction.
figure, imagesc([1:1:size(test_file_freq,2)], [1:1:size(test_file_freq,1)]/(size(test_file_freq,1))*(info.SampleRate/2),mag2db(abs(test_file_freq)), [-150, 50]), title('Original Speech Spectrum'), view([0, -90])
ylabel('Frequency [Hz]'), xlabel('Time Bin')
figure, imagesc([1:1:size(V_est,2)], [1:1:size(V_est,1)]/(size(V_est,1))*(info.SampleRate/2),mag2db(abs(V_est)), [-150, 50]), title('Reconstructed Speech Spectrum'), view([0, -90])
ylabel('Frequency [Hz]'), xlabel('Time Bin')

%% Phase reconstruction
phase_original = angle(test_file_freq); % Original phase - pre-transmission.
phase_recon = phase_reconstruction(phase_original, kc, 3); %Reconstruction using cross-correlation technique. Lift = 1, Mirror = 2, Xcorr = 3.

V_est_istft = istft(V_est.*exp(1i*phase_recon), M, R, N); %ISTFT reconstruction into time-domain.

disp('Playing Original Speech Signal...')
player = audioplayer(test_file, fs);
playblocking(player)

disp('Playing Reconstructed Speech Signal...')
player = audioplayer(V_est_istft, fs);
playblocking(player)
