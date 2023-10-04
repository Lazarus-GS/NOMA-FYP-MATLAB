clc; 
clear variables; 
close all;

seed = 10;
rng(seed);
total_energy_consumed_OMA = 0;
total_energy_consumed_NOMA = 0;

%transmit power range
Pt = -114:2:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

rate1 = 1; rate2 = 2;       %Target rate of users in bps/Hz

N = 10^4; %Length of bit stream
dgap = [0.1,1,100]; %in meters

un = 7.071;
uf = [21.2132, 15.811, 15.811, 21.2132, 15.811, 15.811, 15.811, 15.811, 21.2132, 15.811, 15.811, 21.2132];

% for jk = 1:4
    
    for m = 1:12 %dgap values loop 
        d2 = un; %near user %Distance of users
        d1 = uf(m); %far user
        eta = 4; %Path loss exponent

        %Rayleigh fading coefficients of both users
        h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); %far
        h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); %near

        %Channel gains
        g1 = (abs(h1)).^2;
        g2 = (abs(h2)).^2;
        ga1(m) = mean((abs(h1)).^2);
        ga2(m) = mean((abs(h2)).^2);
       

        BW = 10^9; %bandwidth
        No = -174 + 10*log10(BW);
        no = (10^-3)*db2pow(No); 

        SNR = Pt - No;

        %Generate noise samples for both users
        w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        w2 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

        %Generate random binary data for two users
        data1 = randi([0 1],1,N);  %Data bits of user 1
        data2 = randi([0 1],1,N);  %Data bits of user 2
 
        %Do BPSK modulation of data
        x1 = 2*data1 - 1;
        x2 = 2*data2 - 1;

        %Pt = 0:2:40;                %Transmit power in dBm
        %pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
        %BW = 10^6;                  %System bandwidth
        %No = -174 + 10*log10(BW);   %Noise power (dBm)
        %no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

        %Power allocation coefficients
        a1 = 0.75; 
        a2 = 1-a1; 
        %a3 = 1-(a1+a2);

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

            %IMPROVED FAIR PA%
            %a1 = epsilon*(no + pt*g1)./(pt*g1*(1+epsilon));
            %a1(a1>1) = 0;
            %a2 = 1 - a1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
            %capacity
            %NOMA capacity calculation
            C_noma_1 = log2(1 + pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); %far user
            C_noma_12 = log2(1 + pt(u)*a1.*g2./(pt(u)*a2.*g2+no));
            C_noma_2 = log2(1 + pt(u)*a2.*g2/no); %near user  


            %gamma_far(u) = mean(C_noma_2);
            C_noma_sum(m,u) = mean(C_noma_1 + C_noma_2);  %Sum capacity of NOMA

            %OMA capacity calculation
            C_oma_1 = (1/2)*log2(1 + pt(u)*g1/no);    %User 1
            C_oma_12 = (1/2)*log2(1 + pt(u)*g1/no);
            C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no);    %User 2

            C_oma_sum(m,u) = mean(C_oma_1 + C_oma_2); %Sum capacity of OMA

            %Average of achievable rates
            R1_av(m,u) = mean(C_noma_1);
            R2_av(m,u) = mean(C_noma_2);
            
            ga1(m,u) = mean((abs(h1)).^2);
            ga2(m,u) = mean((abs(h2)).^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

            %Recieved SINR
            %NOMA recieved SINR

            SINR_noma_1(m,u) = mean(pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); %far user 
            SINR_noma_2(m,u) = mean(pt(u)*a2.*g2/no); %near user

            %OMA Recieved SNR

            SINR_oma_1(m,u) = mean(pt(u)*g1/no);
            SINR_oma_2(m,u) = mean(pt(u)*g2/no);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

            %Outage Probability
            %NOMA 

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


            %OMA 
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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %BER
        %NOMA

            %Do superposition coding
            x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);
            %Received signals
            y1 = h1'.*x + w1;
            y2 = h2'.*x + w2;

            %Equalize 
            eq1 = y1./h1';
            eq2 = y2./h2';

            %AT USER 1--------------------
            %Direct decoding of x1 from y1
            x1_hat = zeros(1,N);
            x1_hat(eq1>0) = 1;

            %Compare decoded x1_hat with data1 to estimate BER
            ber1(m,u) = biterr(data1,x1_hat)/N;

            %----------------------------------

            %AT USER 2-------------------------
            %Direct decoding of x1 from y2
            x12_hat = ones(1,N);
            x12_hat(eq2<0) = -1;

            y2_dash = eq2 - sqrt(a1*pt(u))*x12_hat;
            x2_hat = zeros(1,N);
            x2_hat(real(y2_dash)>0) = 1;

            ber2(m,u) = biterr(x2_hat, data2)/N;
            %-----------------------------------


        %OMA
            %ip = rand(1,N)>0.5; % generating 0,1 with equal probability
            

            %Transmission signals
            xoma1 = sqrt(pt(u))*x1;
            xoma2 = sqrt(pt(u))*x2;

            %Received signals
            y1oma = h1'.*xoma1 + w1;
            y2oma = h2'.*xoma2 + w2;

            % equalization
            y1Hat = y1oma./h1';
            y2Hat = y2oma./h2'; 

            % receiver - hard decision decoding
            ipHat1 = real(y1Hat)>0;
            ipHat2 = real(y2Hat)>0;

            % counting the errors
            nErr1(m,u) = size(find([data1- ipHat1]),2);
            nErr2(m,u) = size(find([data2- ipHat2]),2);

            simBer1(m,u) = nErr1(m,u)/N; % simulated ber
            simBer2(m,u) = nErr2(m,u)/N;
       
        end
        
        SNR = Pt - No;

    end
    
[SNR_grid, uf_grid] = meshgrid(SNR, uf);

figure;
surf(SNR_grid, uf_grid, C_noma_sum, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % NOMA in red
hold on;
surf(SNR_grid, uf_grid, C_oma_sum, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % OMA in blue
xlabel('SNR (dB)');
ylabel('Distance from Access Point (m)');
zlabel('Sum Rate Capacity (bps/Hz)');
title('Sum Rate Capacity for NOMA and OMA over SNR and User Position');
legend('NOMA', 'OMA');
hold off;
    
% Calculate metrics after the loop for OMA and NOMA

spectral_efficiency_OMA = mean(C_oma_sum(:,u)) / BW;
spectral_efficiency_NOMA = mean(C_noma_sum(:,u)) / BW;

data_throughput_OMA = mean(C_oma_sum(:,u));
data_throughput_NOMA = mean(C_noma_sum(:,u));

data_loss_rate_OMA_user1 = simBer1(m,u);
data_loss_rate_OMA_user2 = simBer2(m,u);
data_loss_rate_NOMA_user1 = ber1(m,u);
data_loss_rate_NOMA_user2 = ber2(m,u);

% Display the metrics for OMA and NOMA
disp('Metrics for OMA:');
disp(['Spectral Efficiency: ', num2str(spectral_efficiency_OMA)]);
disp(['Data Throughput: ', num2str(data_throughput_OMA)]);
disp(['Data Loss Rate (User 1): ', num2str(data_loss_rate_OMA_user1)]);
disp(['Data Loss Rate (User 2): ', num2str(data_loss_rate_OMA_user2)]);

disp('Metrics for NOMA:');
disp(['Spectral Efficiency: ', num2str(spectral_efficiency_NOMA)]);
disp(['Data Throughput: ', num2str(data_throughput_NOMA)]);
disp(['Data Loss Rate (User 1): ', num2str(data_loss_rate_NOMA_user1)]);
disp(['Data Loss Rate (User 2): ', num2str(data_loss_rate_NOMA_user2)]);

% Define the position grid
position_grid = {'a1', 'a2', 'a3', 'a4', 
                 'b1', 'b2', 'b3', 'b4',
                 'c1', 'c2', 'c3', 'c4',
                 'd1', 'd2', 'd3', 'd4'};

% Define the near user positions
near_user_positions = {'b2', 'b3', 'c2', 'c3'};

% Define the far user positions based on the uf array
far_user_positions = {'a1', 'a2', 'a3', 'a4', 'b1', 'b4', 'c1', 'c4', 'd1', 'd2', 'd3', 'd4'};

% Define line styles for NOMA and OMA
lineStyles = {'-', '--'}; % Solid for NOMA, Dashed for OMA

% Define markers based on the number of near user positions
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};

% Define a colormap based on the number of far user positions
colors = jet(length(far_user_positions));

% Create a figure
figure;

% Create empty arrays to store legend handles and entries
legendHandles = [];
legendEntries = {};

% Define marker indices to increase the number of markers in a line
markerIndices = 1:1:length(SNR); % Place a marker every 5 data points

% Loop through each near user position
for i = 1:length(near_user_positions)
    % Loop through each far user position
    for j = 1:length(far_user_positions)
        % Extract the sum rate capacity for NOMA and OMA for the current combination
        sum_rate_noma = C_noma_sum(j, :);
        sum_rate_oma = C_oma_sum(j, :);
        
        % Plot the NOMA sum rate capacity vs SNR with solid lines, colormap, specific marker, and increased marker frequency
        hNOMA = plot(SNR, sum_rate_noma, 'LineStyle', lineStyles{1}, 'Marker', markers{i}, 'LineWidth', 1.5, 'Color', colors(j,:), 'MarkerIndices', markerIndices);
        hold on;
        
        % Plot the OMA sum rate capacity vs SNR with dashed lines, colormap, specific marker, and increased marker frequency
        hOMA = plot(SNR, sum_rate_oma, 'LineStyle', lineStyles{2}, 'Marker', markers{i}, 'LineWidth', 1.5, 'Color', colors(j,:), 'MarkerIndices', markerIndices);
        hold on;
        
        % Store handles and entries for the legend
        legendHandles = [legendHandles, hNOMA, hOMA];
        legendEntries{end+1} = ['NOMA: ' near_user_positions{i} ' to ' far_user_positions{j}];
        legendEntries{end+1} = ['OMA: ' near_user_positions{i} ' to ' far_user_positions{j}];
    end
end

% Add labels, title, and legend
xlabel('SNR (dB)');
ylabel('Sum Rate Capacity (bps/Hz)');
title('Sum Rate Capacity vs SNR for Different User Positions');
legend(legendHandles, legendEntries, 'Location', 'eastoutside');

grid on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%
%Sum rate capacity over SNR for different user
% Define a colormap
colors = jet(length(far_user_positions));

markerSpacing = 4; % Adjust this value to change the frequency of markers

for i = 1:length(near_user_positions)
    % Create a new figure for each near user position
    figure;
    
    for j = 1:length(far_user_positions)
        % Extract the sum rate capacity for NOMA and OMA
        sum_rate_noma = C_noma_sum(j, :);
        sum_rate_oma = C_oma_sum(j, :);
        
        % Plot the NOMA sum rate capacity vs SNR with solid lines, colors, and '>' marker
        plot(SNR, sum_rate_noma, 'DisplayName', ['NOMA: ' near_user_positions{i} ' to ' far_user_positions{j}], ...
            'LineStyle', '-', 'LineWidth', 1.5, 'Color', colors(j,:), 'Marker', '>', 'MarkerIndices', 1:markerSpacing:length(SNR));
        hold on;
        
        % Plot the OMA sum rate capacity vs SNR with dashed lines, colors, and 'O' marker
        plot(SNR, sum_rate_oma, 'DisplayName', ['OMA: ' near_user_positions{i} ' to ' far_user_positions{j}], ...
            'LineStyle', '--', 'LineWidth', 1.5, 'Color', colors(j,:), 'Marker', 'O', 'MarkerIndices', 1:markerSpacing:length(SNR));
    end
    
    % Add labels, title, and legend
    xlabel('SNR (dB)');
    ylabel('Sum Rate Capacity (bps/Hz)');
    title(['Sum Rate for Near User at ' near_user_positions{i}]);
    legend('Location', 'eastoutside'); % Place legend outside
    grid on;
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the BER data from the provided code
BER_OMA_near = simBer1;
BER_OMA_far = simBer2;
BER_NOMA_near = ber1;
BER_NOMA_far = ber2;

% Create a figure
figure;

% Plot BER for NOMA - Near User
semilogy(SNR, mean(BER_NOMA_near, 1), '-', 'DisplayName', 'NOMA - Near User', 'LineWidth', 1.5);
hold on;

% Plot BER for NOMA - Far User
semilogy(SNR, mean(BER_NOMA_far, 1), '-', 'DisplayName', 'NOMA - Far User', 'LineWidth', 1.5);

% Plot BER for OMA - Near User
semilogy(SNR, mean(BER_OMA_near, 1), '--', 'DisplayName', 'OMA - Near User', 'LineWidth', 1.5);

% Plot BER for OMA - Far User
semilogy(SNR, mean(BER_OMA_far, 1), '--', 'DisplayName', 'OMA - Far User', 'LineWidth', 1.5);

% Add labels, title, and legend
xlabel('Signal to Noice Ratio (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for Near and Far Users');
legend('Location', 'southwest');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%

%outage probability plot

% Create a figure for the far user
figure;

% Plot outage probability for NOMA - Far User
semilogy(SNR, mean(poutNoma1, 1), '-', 'DisplayName', 'NOMA', 'LineWidth', 1.5);
hold on;

% Plot outage probability for OMA - Far User
semilogy(SNR, mean(poutoma1, 1), '--', 'DisplayName', 'OMA', 'LineWidth', 1.5);

% Add labels, title, and legend
xlabel('SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs SNR for Far User');
legend('Location', 'northeast');
grid on;
hold off;

% Create a figure for the near user
figure;

% Plot outage probability for NOMA - Near User
semilogy(SNR, mean(poutNoma2, 1), '-', 'DisplayName', 'NOMA', 'LineWidth', 1.5);
hold on;

% Plot outage probability for OMA - Near User
semilogy(SNR, mean(poutoma2, 1), '--', 'DisplayName', 'OMA', 'LineWidth', 1.5);

% Add labels, title, and legend
xlabel('SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs SNR for Near User');
legend('Location', 'northeast');
grid on;
hold off;
