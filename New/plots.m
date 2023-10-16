load('simdatagrid12.mat');
writetable(datasetTable, 'simulation_dataset2.csv');
Z_data1 = permute(C_noma_sum, [1 3 2]);
Z_data2 = permute(C_oma_sum, [1 3 2]);

size(SNR_grid)
size(near_grid)
size(far_grid)
size(Z_data1)

% ... [Rest of the code remains unchanged]

% Extracting the unique SNR values
uniqueSNRs = unique(datasetTable.SNR);

% Preallocate matrices for storing sum rates for OMA and NOMA
sumRateOMA = zeros(length(nearUserPositions), length(farUserPositions), length(uniqueSNRs));
sumRateNOMA = zeros(size(sumRateOMA));

% Populate the matrices with sum rate values
for i = 1:length(nearUserPositions)
    for j = 1:length(farUserPositions)
        for k = 1:length(uniqueSNRs)
            currentSNR = uniqueSNRs(k);
            omaRow = datasetTable(strcmp(datasetTable.NearUserPos, nearUserPositions{i}) & ...
                                  strcmp(datasetTable.FarUserPos, farUserPositions{j}) & ...
                                  datasetTable.SNR == currentSNR & ...
                                  strcmp(datasetTable.MAScheme, 'OMA'), :);
            nomaRow = datasetTable(strcmp(datasetTable.NearUserPos, nearUserPositions{i}) & ...
                                    strcmp(datasetTable.FarUserPos, farUserPositions{j}) & ...
                                    datasetTable.SNR == currentSNR & ...
                                    strcmp(datasetTable.MAScheme, 'NOMA'), :);
            sumRateOMA(i, j, k) = mean(omaRow.SumRateNearUser);
            sumRateNOMA(i, j, k) = mean(nomaRow.SumRateNearUser);
        end
    end
end

% Plotting the results
[X, Y] = meshgrid(1:length(farUserPositions), 1:length(nearUserPositions));

for k = 1:length(uniqueSNRs)
    figure;
    surf(X, Y, sumRateOMA(:,:,k));
    title(['OMA Sum Rate for SNR = ' num2str(uniqueSNRs(k)) ' dB']);
    xlabel('Far User Position Index');
    ylabel('Near User Position Index');
    zlabel('Sum Rate');
    colorbar;
    
    figure;
    surf(X, Y, sumRateNOMA(:,:,k));
    title(['NOMA Sum Rate for SNR = ' num2str(uniqueSNRs(k)) ' dB']);
    xlabel('Far User Position Index');
    ylabel('Near User Position Index');
    zlabel('Sum Rate');
    colorbar;
end

%writetable(datasetTable, 'simulation_dataset2.csv');

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

% Update the plots
% For the 3D surface plot
[SNR_grid, near_grid, far_grid] = meshgrid(SNR, nearUserDistances, farUserDistances);

figure;
surf(SNR_grid, near_grid, Z_data1, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % NOMA in red
hold on;
surf(SNR_grid, far_grid, Z_data2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % OMA in blue
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
