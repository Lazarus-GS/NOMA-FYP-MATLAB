clc; 
clear variables; 
close all;

columnNames = {'NearUserPos', 'FarUserPos', 'FarUserDist', 'TransmitPower', 'SNR', 'SumRateFarUser', 'SumRateNearUser', 'BERFarUser', 'BERNearUser', 'OutageProbFarUser', 'OutageProbNearUser', 'MAScheme'};
datasetTable = table('Size', [0 12], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string'}, 'VariableNames', columnNames);

% Define the position grid
position_grid = {'a1', 'a2', 'a3', 'a4', 
                 'b1', 'b2', 'b3', 'b4',
                 'c1', 'c2', 'c3', 'c4',
                 'd1', 'd2', 'd3', 'd4'};

nearUserPositions = {'b2', 'b3', 'c2', 'c3'};
farUserPositions = {'a1', 'a2', 'a3', 'a4', 'b1', 'b4', 'c1', 'c4', 'd1', 'd2', 'd3', 'd4'};

seed = 10;
rng(seed);

Pt = -114:2:-54;	
pt = db2pow(Pt);	

rate1 = 1; rate2 = 2;       

N = 10^4; 

un = 7.071;
uf = [21.2132, 15.811, 15.811, 21.2132, 15.811, 15.811, 15.811, 15.811, 21.2132, 15.811, 15.811, 21.2132];

%Preallocating arrays
C_noma_sum = zeros(length(farUserPositions), length(Pt));
C_oma_sum = zeros(length(farUserPositions), length(Pt));
R1_av = zeros(length(farUserPositions), length(Pt));
R2_av = zeros(length(farUserPositions), length(Pt));
SINR_noma_1 = zeros(length(farUserPositions), length(Pt));
SINR_noma_2 = zeros(length(farUserPositions), length(Pt));
SINR_oma_1 = zeros(length(farUserPositions), length(Pt));
SINR_oma_2 = zeros(length(farUserPositions), length(Pt));
ga1 = zeros(length(farUserPositions), length(Pt));
ga2 = zeros(length(farUserPositions), length(Pt));
poutNoma1 = zeros(length(farUserPositions), length(Pt));
poutNoma2 = zeros(length(farUserPositions), length(Pt));
poutoma1 = zeros(length(farUserPositions), length(Pt));
poutoma2 = zeros(length(farUserPositions), length(Pt));
ber1 = zeros(length(farUserPositions), length(Pt));
ber2 = zeros(length(farUserPositions), length(Pt));
nErr1 = zeros(length(farUserPositions), length(Pt));
nErr2 = zeros(length(farUserPositions), length(Pt));
simBer1 = zeros(length(farUserPositions), length(Pt));
simBer2 = zeros(length(farUserPositions), length(Pt));

for nearIdx = 1:length(nearUserPositions)
    for farIdx = 1:length(farUserPositions)
        for u = 1:length(Pt)
            d2 = un; 
            d1 = uf(farIdx); 
            eta = 4;

            h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); 
            h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); 

            g1 = (abs(h1)).^2;
            g2 = (abs(h2)).^2;

            BW = 10^9; 
            No = -174 + 10*log10(BW);
            no = (10^-3)*db2pow(No); 

            w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
            w2 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

            data1 = randi([0 1],1,N);  
            data2 = randi([0 1],1,N);  

            x1 = 2*data1 - 1;
            x2 = 2*data2 - 1;

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

            C_noma_1 = log2(1 + pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); 
            C_noma_12 = log2(1 + pt(u)*a1.*g2./(pt(u)*a2.*g2+no));
            C_noma_2 = log2(1 + pt(u)*a2.*g2/no); 

            C_noma_sum(farIdx,u) = mean(C_noma_1 + C_noma_2);  

            C_oma_1 = (1/2)*log2(1 + pt(u)*g1/no);   
            C_oma_12 = (1/2)*log2(1 + pt(u)*g1/no);
            C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no);   

            C_oma_sum(farIdx,u) = mean(C_oma_1 + C_oma_2); 

            R1_av(farIdx,u) = mean(C_noma_1);
            R2_av(farIdx,u) = mean(C_noma_2);
            
            ga1(farIdx,u) = mean((abs(h1)).^2);
            ga2(farIdx,u) = mean((abs(h2)).^2);       

            SINR_noma_1(farIdx,u) = mean(pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); 
            SINR_noma_2(farIdx,u) = mean(pt(u)*a2.*g2/no);

            SINR_oma_1(farIdx,u) = mean(pt(u)*g1/no);
            SINR_oma_2(farIdx,u) = mean(pt(u)*g2/no);
      
            for k = 1:N
                if C_noma_1(k) < rate1
                    pn1(u) = pn1(u)+1;
                end
                if (C_noma_12(k) < rate1)||(C_noma_2(k) < rate2)
                    pn2(u) = pn2(u)+1;
                end
            end

            poutNoma1(farIdx,u) = pn1(u)/N;
            poutNoma2(farIdx,u) = pn2(u)/N;

            for ko = 1:N
                if C_oma_1(ko) < rate1
                    po1(u) = po1(u)+1;
                end
                if (C_oma_2(ko) < rate2)
                    po2(u) = po2(u)+1;
                end
            end

            poutoma1(farIdx,u) = po1(u)/N;
            poutoma2(farIdx,u) = po2(u)/N;

            x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);

            y1 = h1'.*x + w1;
            y2 = h2'.*x + w2;

            eq1 = y1./h1';
            eq2 = y2./h2';

            x1_hat = zeros(1,N);
            x1_hat(eq1>0) = 1;

            ber1(farIdx,u) = biterr(data1,x1_hat)/N;

            x12_hat = ones(1,N);
            x12_hat(eq2<0) = -1;

            y2_dash = eq2 - sqrt(a1*pt(u))*x12_hat;
            x2_hat = zeros(1,N);
            x2_hat(real(y2_dash)>0) = 1;

            ber2(farIdx,u) = biterr(x2_hat, data2)/N;
            
            xoma1 = sqrt(pt(u))*x1;
            xoma2 = sqrt(pt(u))*x2;

            y1oma = h1'.*xoma1 + w1;
            y2oma = h2'.*xoma2 + w2;

            y1Hat = y1oma./h1';
            y2Hat = y2oma./h2'; 

            ipHat1 = real(y1Hat)>0;
            ipHat2 = real(y2Hat)>0;

            nErr1(farIdx,u) = size(find([data1- ipHat1]),2);
            nErr2(farIdx,u) = size(find([data2- ipHat2]),2);

            simBer1(farIdx,u) = nErr1(farIdx,u)/N; 
            simBer2(farIdx,u) = nErr2(farIdx,u)/N;

            SNR = Pt(u) - No;

            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, uf(farIdx), Pt(u), SNR, C_oma_sum(farIdx,u), C_oma_sum(nearIdx,u), simBer1(farIdx,u), simBer2(farIdx,u), poutoma1(farIdx,u), poutoma2(farIdx,u), 'OMA'};
            datasetTable = [datasetTable; newRow];
            
            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, uf(farIdx), Pt(u), SNR, C_noma_sum(farIdx,u), C_noma_sum(nearIdx,u), ber1(farIdx,u), ber2(farIdx,u), poutNoma1(farIdx,u), poutNoma2(farIdx,u), 'NOMA'};
            datasetTable = [datasetTable; newRow];
        end
        
        SNR = Pt - No;
    end
end

% Calculate metrics after the loop for OMA and NOMA

spectral_efficiency_OMA = mean(C_oma_sum(:,u)) / BW;
spectral_efficiency_NOMA = mean(C_noma_sum(:,u)) / BW;

data_throughput_OMA = mean(C_oma_sum(:,u));
data_throughput_NOMA = mean(C_noma_sum(:,u));

data_loss_rate_OMA_user1 = simBer1(farIdx,u);
data_loss_rate_OMA_user2 = simBer2(farIdx,u);
data_loss_rate_NOMA_user1 = ber1(farIdx,u);
data_loss_rate_NOMA_user2 = ber2(farIdx,u);

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


[SNR_grid, uf_grid] = meshgrid(SNR, uf);

disp(size(SNR_grid));
disp(size(uf_grid));
disp(size(C_noma_sum));


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


lineStyles = {'-', '--'}; % Solid for NOMA, Dashed for OMA
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
colors = jet(length(farUserPositions));

figure;

legendHandles = [];
legendEntries = {};

markerIndices = 1:1:length(SNR); % Place a marker every 5 data points

for i = 1:length(nearUserPositions)
    for j = 1:length(farUserPositions)
        sum_rate_noma = C_noma_sum(j, :);
        sum_rate_oma = C_oma_sum(j, :);
        
        hNOMA = plot(SNR, sum_rate_noma, 'LineStyle', lineStyles{1}, 'Marker', markers{i}, 'LineWidth', 1.5, 'Color', colors(j,:), 'MarkerIndices', markerIndices);
        hold on;
        hOMA = plot(SNR, sum_rate_oma, 'LineStyle', lineStyles{2}, 'Marker', markers{i}, 'LineWidth', 1.5, 'Color', colors(j,:), 'MarkerIndices', markerIndices);
        hold on;
        
        legendHandles = [legendHandles, hNOMA, hOMA];
        legendEntries{end+1} = ['NOMA: ' nearUserPositions{i} ' to ' farUserPositions{j}];
        legendEntries{end+1} = ['OMA: ' nearUserPositions{i} ' to ' farUserPositions{j}];
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
% Create a new figure before starting the loop
figure;
colors = jet(length(farUserPositions));
% Loop through each near user position
for i = 1:length(nearUserPositions)
    
    % Specify the current subplot
    subplot(2, 2, i); % Assuming you have 4 near user positions, this will create a 2x2 grid of subplots
    
    for j = 1:length(farUserPositions)

        sum_rate_noma = C_noma_sum(j, :);
        sum_rate_oma = C_oma_sum(j, :);
        
        plot(SNR, sum_rate_noma, 'DisplayName', ['NOMA: ' nearUserPositions{i} ' to ' farUserPositions{j}], ...
            'LineStyle', '-', 'LineWidth', 1.5, 'Color', colors(j,:), 'Marker', '>');
        hold on;
        plot(SNR, sum_rate_oma, 'DisplayName', ['OMA: ' nearUserPositions{i} ' to ' farUserPositions{j}], ...
            'LineStyle', '--', 'LineWidth', 1.5, 'Color', colors(j,:), 'Marker', 'O');
    end
    
    % Add labels, title, and legend to each subplot
    xlabel('SNR (dB)');
    ylabel('Sum Rate Capacity (bps/Hz)');
    title(['Sum Rate for Near User at ' nearUserPositions{i}]);
    legend('Location', 'eastoutside'); % Place legend outside
    grid on;
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a figure
figure;

semilogy(SNR, mean(ber1, 1), '-', 'DisplayName', 'NOMA - Near User', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(ber2, 1), '-', 'DisplayName', 'NOMA - Far User', 'LineWidth', 1.5);
semilogy(SNR, mean(simBer1, 1), '--', 'DisplayName', 'OMA - Near User', 'LineWidth', 1.5);
semilogy(SNR, mean(simBer2, 1), '--', 'DisplayName', 'OMA - Far User', 'LineWidth', 1.5);

xlabel('Signal to Noice Ratio (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for Near and Far Users');
legend('Location', 'southwest');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%

%outage probability plot

% Create a figure for the combined outage probability
figure;

semilogy(SNR, mean(poutNoma1, 1), '-', 'DisplayName', 'NOMA - Far User', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(poutoma1, 1), '--', 'DisplayName', 'OMA - Far User', 'LineWidth', 1.5);
semilogy(SNR, mean(poutNoma2, 1), '-', 'DisplayName', 'NOMA - Near User', 'LineWidth', 1.5);
semilogy(SNR, mean(poutoma2, 1), '--', 'DisplayName', 'OMA - Near User', 'LineWidth', 1.5);

xlabel('SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs SNR for Both Near and Far Users');
legend('Location', 'northeast');
grid on;
hold off;
