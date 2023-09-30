clc; 
clear variables; 
close all;

seed = 10;
rng(seed);

% LDPC parameters
H = dvbs2ldpc(1/2);
K = size(H, 2) - size(H, 1);
ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H, 'MaximumIterationCount', 50, 'DecisionMethod', 'Soft decision');

% transmit power range
Pt = -114:5:-54; % in dB
pt = db2pow(Pt); % in linear scale

rate1 = 1; rate2 = 2; % Target rate of users in bps/Hz

N = 10^2; % Length of bit stream
dgap = [0.1,1,100]; % in meters

un = 7.071;
uf = [21.2132, 15.811, 15.811, 21.2132, 15.811, 15.811, 15.811, 15.811, 21.2132, 15.811, 15.811, 21.2132];

for m = 1:12 % dgap values loop 
    d2 = un; % near user % Distance of users
    d1 = uf(m); % far user
    eta = 4; % Path loss exponent

    % Rayleigh fading coefficients of both users
    h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); % far
    h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); % near

    % Channel gains
    g1 = (abs(h1)).^2;
    g2 = (abs(h2)).^2;
    ga1(m) = mean((abs(h1)).^2);
    ga2(m) = mean((abs(h2)).^2);

    BW = 10^9; % bandwidth
    No = -174 + 10*log10(BW);
    no = (10^-3)*db2pow(No); 

    SNR = Pt - No;

    % Generate noise samples for both users
    w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
    w2 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

    % Generate random binary data for two users
    data1 = randi([0 1],K,1); % Data bits of user 1
    data2 = randi([0 1],K,1); % Data bits of user 2

    % Encode data using LDPC
    encodedData1 = ldpcEncoder(data1);
    encodedData2 = ldpcEncoder(data2);

    % Do BPSK modulation of encoded data
    x1 = 2*encodedData1 - 1;
    x2 = 2*encodedData2 - 1;

    % Power allocation coefficients
    a1 = 0.75; 
    a2 = 1-a1;

    C_noma = zeros(12,length(pt));
    C_oma = zeros(12,length(pt));
    SINR_noma = zeros(12,length(pt));
    SINR_oma = zeros(12,length(pt));

    p = length(Pt);
    pn1 = zeros(12,length(Pt));
    pn2 = zeros(12,length(Pt));
    po1 = zeros(12,length(Pt));
    po2 = zeros(12,length(Pt));

    for u = 1:p
        % capacity
        % NOMA capacity calculation
        C_noma_1 = log2(1 + pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); % far user
        C_noma_12 = log2(1 + pt(u)*a1.*g2./(pt(u)*a2.*g2+no));
        C_noma_2 = log2(1 + pt(u)*a2.*g2/no); % near user  

        C_noma_sum(m,u) = mean(C_noma_1 + C_noma_2); % Sum capacity of NOMA

        % OMA capacity calculation
        C_oma_1 = (1/2)*log2(1 + pt(u)*g1/no); % User 1
        C_oma_12 = (1/2)*log2(1 + pt(u)*g1/no);
        C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no); % User 2

        C_oma_sum(m,u) = mean(C_oma_1 + C_oma_2); % Sum capacity of OMA

        % Average of achievable rates
        R1_av(m,u) = mean(C_noma_1);
        R2_av(m,u) = mean(C_noma_2);

        ga1(m,u) = mean((abs(h1)).^2);
        ga2(m,u) = mean((abs(h2)).^2);

        % Recieved SINR
        % NOMA recieved SINR
        SINR_noma_1(m,u) = mean(pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); % far user 
        SINR_noma_2(m,u) = mean(pt(u)*a2.*g2/no); % near user

        % OMA Recieved SNR
        SINR_oma_1(m,u) = mean(pt(u)*g1/no);
        SINR_oma_2(m,u) = mean(pt(u)*g2/no);

        % Outage Probability
        % NOMA 
        for k = 1:N
            if C_noma_1(k) < rate1
                pn1(u) = pn1(u)+1;
            end
            if (C_noma_12(k) < rate1)||(C_noma_2(k) < rate2)
                pn2(u) = pn2(u)+1;
            end
        end

        poutNoma1(m,u) = pn1(u)/N;
        poutNoma2(m,u) = pn2(u)/N;

        % OMA 
        for ko = 1:N
            if C_oma_1(ko) < rate1
                po1(u) = po1(u)+1;
            end
            if (C_oma_2(ko) < rate2)
                po2(u) = po2(u)+1;
            end
        end

        poutoma1(m,u) = po1(u)/N;
        poutoma2(m,u) = po2(u)/N;

        % BER
        % NOMA
        % Do superposition coding
        x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);
        % Received signals
        y1 = h1'.*x + w1;
        y2 = h2'.*x + w2;

        % Equalize 
        eq1 = y1./h1';
        eq2 = y2./h2';

        % LDPC Decoding
        decodedData1 = zeros(32400, size(eq1, 2));  % Assuming 32400 is the length of the decoded data
        for col = 1:size(eq1, 2)
            release(ldpcDecoder);
            decodedData1(:, col) = ldpcDecoder(eq1(:, col));
        end
        decodedData2 = zeros(32400, size(eq1, 2));
        for col = 1:size(eq2, 2)
            release(ldpcDecoder);
            decodedData2(:, col) = ldpcDecoder(eq2(:, col));
        end

        % Convert decoded data to bits for BER calculation
        x1_hat = double(decodedData1 > 0).';
        x2_hat = double(decodedData2 > 0).';

        % Adjust the size of the decoded data to match the size of the encoded data
        decodedLength = N; % or any other appropriate value
        if length(x1_hat(:, u)) < decodedLength
            x1_hat_padded = [x1_hat(:, u); zeros(decodedLength - length(x1_hat(:, u)), 1)];
        else
            x1_hat_padded = x1_hat(1:decodedLength, u);
        end

        if length(x2_hat(:, u)) < decodedLength
            x2_hat_padded = [x2_hat(:, u); zeros(decodedLength - length(x2_hat(:, u)), 1)];
        else
            x2_hat_padded = x2_hat(1:decodedLength, u);
        end

        % Adjust the size of encodedData1 to match x1_hat_padded
        if length(encodedData1) < length(x1_hat_padded)
            encodedData1_padded = [encodedData1; zeros(length(x1_hat_padded) - length(encodedData1), 1)];
        else
            encodedData1_padded = encodedData1(1:length(x1_hat_padded));
        end

        % Adjust the size of encodedData2 to match x2_hat_padded
        if length(encodedData2) < length(x2_hat_padded)
            encodedData2_padded = [encodedData2; zeros(length(x2_hat_padded) - length(encodedData2), 1)];
        else
            encodedData2_padded = encodedData2(1:length(x2_hat_padded));
        end

        % Now, use the padded matrices in the biterr function
        ber1(m,u) = biterr(encodedData1_padded, x1_hat_padded)/N;
        ber2(m,u) = biterr(encodedData2_padded, x2_hat_padded)/N;


        % OMA
        % Transmission signals
        xoma1 = sqrt(pt(u))*x1;
        xoma2 = sqrt(pt(u))*x2;

        % Received signals
        y1oma = h1'.*xoma1 + w1;
        y2oma = h2'.*xoma2 + w2;

        % equalization
        y1Hat = y1oma./h1';
        y2Hat = y2oma./h2'; 

        % LDPC Decoding for OMA
        decodedData1_oma = zeros(32400, size(y1Hat, 2));
        decodedData2_oma = zeros(32400, size(y2Hat, 2));


        for col = 1:size(y1Hat, 2)
            release(ldpcDecoder);
            decodedData1_oma(:, col) = ldpcDecoder(y1Hat(:, col));
        end

        for col = 1:size(y2Hat, 2)
            release(ldpcDecoder);
            decodedData2_oma(:, col) = ldpcDecoder(y2Hat(:, col));
        end

        % Convert decoded data to bits for BER calculation
        x1_hat_oma = double(decodedData1_oma > 0).';
        x2_hat_oma = double(decodedData2_oma > 0).';

        nErr1(m,u) = size(find([encodedData1 - x1_hat_oma(:, u)]),2);
        nErr2(m,u) = size(find([encodedData2 - x2_hat_oma(:, u)]),2);

        simBer1(m,u) = nErr1(m,u)/N; % simulated ber
        simBer2(m,u) = nErr2(m,u)/N;
    end
end
