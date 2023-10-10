% Initialization
clc;
clear variables;
close all;

% Constants
SEED = 10;
TRANSMIT_POWER_RANGE_DB = -114:5:-54;
TRANSMIT_POWER_LINEAR = db2pow(TRANSMIT_POWER_RANGE_DB);
TARGET_RATES = [1, 2]; % bps/Hz
SAMPLE_SIZE = 10^4;
BANDWIDTH = 10^9; % Hz
NOISE_POWER_DB = -174 + 10*log10(BANDWIDTH);
NOISE_POWER_LINEAR = (10^-3) * db2pow(NOISE_POWER_DB);
SNR = TRANSMIT_POWER_RANGE_DB - NOISE_POWER_DB;
PATH_LOSS_EXPONENT = 4;

% Set random seed for reproducibility
rng(SEED);

% Generate user positions
numPositionPoints = 100; % Number of points on the circle for each user
[nearUserPositions, farUserPositions] = generateUserPositions(numPositionPoints);

% Initialize matrices to store results
C_noma_avg = zeros(numPositionPoints, length(TRANSMIT_POWER_RANGE_DB));
C_oma_avg = zeros(numPositionPoints, length(TRANSMIT_POWER_RANGE_DB));

for i = 1:numPositionPoints
    for j = 1:numPositionPoints
        % Compute distance to the access point (origin)
        NEAR_USER_DISTANCE = norm(nearUserPositions(i, :));
        FAR_USER_DISTANCE = norm(farUserPositions(j, :));

        % Compute Rayleigh fading coefficients
        h_far = computeRayleighCoefficients(FAR_USER_DISTANCE, PATH_LOSS_EXPONENT, SAMPLE_SIZE);
        h_near = computeRayleighCoefficients(NEAR_USER_DISTANCE, PATH_LOSS_EXPONENT, SAMPLE_SIZE);
        
        % Generate noise samples
        w_far = generateNoiseSamples(NOISE_POWER_LINEAR, SAMPLE_SIZE);
        w_near = generateNoiseSamples(NOISE_POWER_LINEAR, SAMPLE_SIZE);
        
        % Generate random binary data and perform BPSK modulation
        data_far = generateBinaryData(SAMPLE_SIZE);
        data_near = generateBinaryData(SAMPLE_SIZE);
        x_far = bpskModulation(data_far);
        x_near = bpskModulation(data_near);
        
        % Compute performance metrics
        [C_noma, C_oma, ~, ~, ~, ~, ~, ~] = computePerformanceMetrics(h_far, h_near, w_far, w_near, x_far(:), x_near(:), TRANSMIT_POWER_LINEAR, TARGET_RATES, SAMPLE_SIZE);
        
        C_noma_avg(i, :) = C_noma_avg(i, :) + C_noma;
        C_oma_avg(i, :) = C_oma_avg(i, :) + C_oma;
    end
end

% Average over all user positions
C_noma_avg = C_noma_avg / numPositionPoints;
C_oma_avg = C_oma_avg / numPositionPoints;

% Plotting
figure;
plot(SNR, mean(C_noma_avg, 1), 'linewidth', 2); hold on;
plot(SNR, mean(C_oma_avg, 1), '--', 'linewidth', 2);
xlabel('SNR (dB)');
ylabel('Achievable sum rate (bps/Hz)');
title('Capacity of NOMA and OMA');
legend('NOMA', 'OMA');
grid on;

% Function definitions
function h = computeRayleighCoefficients(d, eta, N)
    h = (sqrt(d^-eta)) * (randn(N, 1) + 1i*randn(N, 1)) / sqrt(2);
end

function w = generateNoiseSamples(no, N)
    w = sqrt(no) * (randn(1, N) + 1i*randn(1, N)) / sqrt(2);
end

function data = generateBinaryData(N)
    data = randi([0, 1], 1, N);
end

function x = bpskModulation(data)
    x = 2 * data - 1;
end

function [nearUserPos, farUserPos] = generateUserPositions(numPoints)
    % Inner and outer circle radii
    r_inner = 10;
    r_outer = 30;

    % Generate near user positions
    theta_near = linspace(0, 2*pi, numPoints); % Equally spaced angles
    x_near = r_inner * cos(theta_near);
    y_near = r_inner * sin(theta_near);
    nearUserPos = [x_near', y_near'];

    % Generate far user positions
    theta_far = linspace(0, 2*pi, numPoints); % Equally spaced angles
    x_far = r_outer * cos(theta_far);
    y_far = r_outer * sin(theta_far);
    farUserPos = [x_far', y_far'];
end

function [C_noma_sum_m, C_oma_sum_m, SINR_noma, SINR_oma, poutNoma, poutOma, berNoma, berOma] = computePerformanceMetrics(h_far, h_near, w_far, w_near, x_far, x_near, TRANSMIT_POWER_LINEAR, TARGET_RATES, SAMPLE_SIZE)
    % Extract channel gains
    g_far = abs(h_far).^2;
    g_near = abs(h_near).^2;

    % Initialize matrices for storing results
    C_noma_sum_m = zeros(1, length(TRANSMIT_POWER_LINEAR));
    C_oma_sum_m = zeros(1, length(TRANSMIT_POWER_LINEAR));
    SINR_noma = zeros(2, length(TRANSMIT_POWER_LINEAR)); % [far_user, near_user]
    SINR_oma = zeros(2, length(TRANSMIT_POWER_LINEAR));  % [far_user, near_user]
    poutNoma = zeros(2, length(TRANSMIT_POWER_LINEAR));  % [far_user, near_user]
    poutOma = zeros(2, length(TRANSMIT_POWER_LINEAR));   % [far_user, near_user]
    berNoma = zeros(2, length(TRANSMIT_POWER_LINEAR));   % [far_user, near_user]
    berOma = zeros(2, length(TRANSMIT_POWER_LINEAR));    % [far_user, near_user]

    % Power allocation coefficients for NOMA
    a1 = 0.75;
    a2 = 0.25;

    for u = 1:length(TRANSMIT_POWER_LINEAR)
        % NOMA capacity calculation
        C_noma_1 = log2(1 + TRANSMIT_POWER_LINEAR(u)*a1.*g_far./(TRANSMIT_POWER_LINEAR(u)*a2.*g_far + 1)); % far user
        C_noma_2 = log2(1 + TRANSMIT_POWER_LINEAR(u)*a2.*g_near); % near user
        C_noma_sum_m(u) = mean(C_noma_1 + C_noma_2);  % Sum capacity of NOMA

        % OMA capacity calculation
        C_oma_1 = 0.5 * log2(1 + TRANSMIT_POWER_LINEAR(u)*g_far); % far user
        C_oma_2 = 0.5 * log2(1 + TRANSMIT_POWER_LINEAR(u)*g_near); % near user
        C_oma_sum_m(u) = mean(C_oma_1 + C_oma_2); % Sum capacity of OMA

        % SINR calculations
        SINR_noma(1, u) = mean(TRANSMIT_POWER_LINEAR(u)*a1.*g_far./(TRANSMIT_POWER_LINEAR(u)*a2.*g_far + 1)); % NOMA far user
        SINR_noma(2, u) = mean(TRANSMIT_POWER_LINEAR(u)*a2.*g_near); % NOMA near user
        SINR_oma(1, u) = mean(TRANSMIT_POWER_LINEAR(u)*g_far); % OMA far user
        SINR_oma(2, u) = mean(TRANSMIT_POWER_LINEAR(u)*g_near); % OMA near user

        % Outage Probability calculations
        poutNoma(1, u) = sum(C_noma_1 < TARGET_RATES(1)) / SAMPLE_SIZE; % NOMA far user
        poutNoma(2, u) = sum(C_noma_2 < TARGET_RATES(2)) / SAMPLE_SIZE; % NOMA near user
        poutOma(1, u) = sum(C_oma_1 < TARGET_RATES(1)) / SAMPLE_SIZE;   % OMA far user
        poutOma(2, u) = sum(C_oma_2 < TARGET_RATES(2)) / SAMPLE_SIZE;   % OMA near user

        % BER calculations
        % NOMA
        y_far_noma = h_far .* x_far + w_far(:);
        y_near_noma = h_near .* (sqrt(a1)*x_far + sqrt(a2)*x_near) + w_near(:);
        decoded_far_noma = real(y_far_noma) > 0;
        decoded_near_noma = real(y_near_noma - sqrt(a1*TRANSMIT_POWER_LINEAR(u)*decoded_far_noma)) > 0;
        berNoma(1, u) = sum(decoded_far_noma ~= x_far) / SAMPLE_SIZE; % NOMA far user
        berNoma(2, u) = sum(decoded_near_noma ~= x_near) / SAMPLE_SIZE; % NOMA near user

        % OMA
        y_far_oma = h_far .* x_far + w_far(:);
        y_near_oma = h_near .* x_near + w_near(:);
        decoded_far_oma = real(y_far_oma) > 0;
        decoded_near_oma = real(y_near_oma) > 0;
        berOma(1, u) = sum(decoded_far_oma ~= x_far(:)) / SAMPLE_SIZE; % OMA far user
        berOma(2, u) = sum(decoded_near_oma ~= x_near(:)) / SAMPLE_SIZE; % OMA near user
    end
end