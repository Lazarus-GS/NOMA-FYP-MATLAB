clc;
clear variables;
close all;

% Constants
GRID_SIZE = 10; % 10 meters per block
NUM_BLOCKS = 4;
INNER_BLOCKS = [2, 3, 6, 7]; % Indices of inner 4 blocks
OUTER_BLOCKS = setdiff(1:NUM_BLOCKS^2, INNER_BLOCKS); % Indices of outer 12 blocks
TRANSMIT_POWER_RANGE_DB = -90:5:0;
TRANSMIT_POWER_LINEAR = db2pow(TRANSMIT_POWER_RANGE_DB);
SAMPLE_SIZE = 10^4;
BANDWIDTH = 10^9; % Hz
NOISE_POWER_DB = -174 + 10*log10(BANDWIDTH);
NOISE_POWER_LINEAR = (10^-3) * db2pow(NOISE_POWER_DB);
SNR = TRANSMIT_POWER_RANGE_DB - NOISE_POWER_DB;
PATH_LOSS_EXPONENT = 4;

% Function to calculate the distance from the center
getDistanceFromCenter = @(blockIndex) GRID_SIZE * sqrt((mod(blockIndex-1, NUM_BLOCKS) - 1.5)^2 + (floor((blockIndex-1)/NUM_BLOCKS) - 1.5)^2);

% Main simulation loop
C_noma_sum_avg = zeros(1, length(TRANSMIT_POWER_RANGE_DB));
C_oma_sum_avg = zeros(1, length(TRANSMIT_POWER_RANGE_DB));

% Loop over all possible user position combinations
for nearIdx = INNER_BLOCKS
    for farIdx = OUTER_BLOCKS
        NEAR_USER_DISTANCE = getDistanceFromCenter(nearIdx);
        FAR_USER_DISTANCE = getDistanceFromCenter(farIdx);
        [C_noma, C_oma] = simulateUserPair(NEAR_USER_DISTANCE, FAR_USER_DISTANCE, TRANSMIT_POWER_LINEAR, NOISE_POWER_LINEAR, PATH_LOSS_EXPONENT, SAMPLE_SIZE);
        C_noma_sum_avg = C_noma_sum_avg + C_noma;
        C_oma_sum_avg = C_oma_sum_avg + C_oma;
    end
end

% Average the results
numCombinations = length(INNER_BLOCKS) * length(OUTER_BLOCKS);
C_noma_sum_avg = C_noma_sum_avg / numCombinations;
C_oma_sum_avg = C_oma_sum_avg / numCombinations;

% Plotting
figure;
plot(SNR, C_noma_sum_avg, '-o', 'linewidth', 2); hold on;
plot(SNR, C_oma_sum_avg, '--x', 'linewidth', 2);
xlabel('SNR (dB)');
ylabel('Average Achievable Sum Rate (bps/Hz)');
title('Average Capacity of NOMA and OMA over User Combinations');
legend('NOMA', 'OMA');
grid on;

% Supporting function: Simulate for a specific user pair
function [C_noma_sum, C_oma_sum] = simulateUserPair(d_near, d_far, pt, no, eta, N)
    h_near = (sqrt(d_near^-eta)) * (randn(N, 1) + 1i*randn(N, 1)) / sqrt(2);
    h_far = (sqrt(d_far^-eta)) * (randn(N, 1) + 1i*randn(N, 1)) / sqrt(2);
    g_near = abs(h_near).^2;
    g_far = abs(h_far).^2;

    a_far = 0.75;
    a_near = 1 - a_far;

    % Capacity calculations using element-wise multiplication
    C_noma_far = log2(1 + pt.*a_far.*g_far./(pt.*a_near.*g_far + no));
    C_noma_near = log2(1 + pt.*a_near.*g_near./no);
    C_noma_sum = mean(C_noma_far + C_noma_near);
    
    C_oma_far = 0.5 * log2(1 + pt.*g_far./no);
    C_oma_near = 0.5 * log2(1 + pt.*g_near./no);
    C_oma_sum = mean(C_oma_far + C_oma_near);
end
