clc; 
clear variables; 
close all;

% Define the table structure
columnNames = {'NearUserPos', 'FarUserPos', 'NearUserDist', 'FarUserDist', 'TransmitPower', 'SNR', 'SumRateFarUser', 'SumRateNearUser','OutageProbFarUser', 'OutageProbNearUser', 'MAScheme'};
datasetTable = table('Size', [0 11], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string'}, 'VariableNames', columnNames);

% Step 2: Define the Position Grid
gridSize = 12;
blockLength = 10; % Length of a single block

% Generate position labels and coordinates directly
% Define the coordinates for each position in the grid
coordinates = zeros(gridSize, gridSize, 2);
for i = 1:gridSize
    for j = 1:gridSize
        coordinates(i, j, :) = [10*(i-6.5), 10*(j-6.5)];
    end
end

% Generate position labels
positionLabels = cell(gridSize, gridSize);
for i = 1:gridSize
    for j = 1:gridSize
        positionLabels{i, j} = [char(i + 96), num2str(j)];
    end
end

% Step 3: Define Near and Far User Positions
innerBoundary = 4;
outerBoundary = 9;

nearUserPositions = {};
farUserPositions = {};
for i = 1:gridSize
    for j = 1:gridSize
        if (i >= innerBoundary && i <= outerBoundary) && (j >= innerBoundary && j <= outerBoundary)
            nearUserPositions = [nearUserPositions; positionLabels{i, j}];
        else
            farUserPositions = [farUserPositions; positionLabels{i, j}];
        end
    end
end

% Create a logical matrix to determine if a position is near the user
isNearUser = false(gridSize, gridSize);
for i = 1:length(nearUserPositions)
    [row, col] = find(strcmp(positionLabels, nearUserPositions{i}));
    isNearUser(row, col) = true;
end

% Flatten the coordinates and isNearUser matrices
flatCoordinates = reshape(coordinates, [], 2);
flatIsNearUser = isNearUser(:);

% Extract near and far user coordinates using logical indexing
nearUserCoordinates = flatCoordinates(flatIsNearUser, :);
farUserCoordinates = flatCoordinates(~flatIsNearUser, :);

% Calculate distances for near and far users
nearUserDistances = sqrt(sum(nearUserCoordinates.^2, 2));
farUserDistances = sqrt(sum(farUserCoordinates.^2, 2));


seed = 10;
rng(seed);

Pt = -114:2:-54;	
pt = db2pow(Pt);	

rate1 = 1; rate2 = 2;       

N = 10^4; 

BW = 10^9; 
No = -174 + 10*log10(BW);
no = (10^-3)*db2pow(No); 

a1 = 0.75; 
a2 = 1-a1; 

C_noma_sum = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
C_oma_sum = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
SINR_noma_1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
SINR_noma_2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
SINR_oma_1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
SINR_oma_2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
poutNoma1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
poutNoma2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
poutoma1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
poutoma2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
pn1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
pn2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
po1 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));
po2 = zeros(length(nearUserPositions), length(farUserPositions), length(Pt));

for nearIdx = 1:length(nearUserPositions)
    for farIdx = 1:length(farUserPositions)
        for u = 1:length(Pt)
            d_near = nearUserDistances(nearIdx);
            d_far = farUserDistances(farIdx);
            eta = 4;

            h1 = (sqrt(d_far^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); 
            h2 = (sqrt(d_near^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); 
            g1 = (abs(h1)).^2;
            g2 = (abs(h2)).^2;

            w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
            w2 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

            data1 = randi([0 1],1,N);  
            data2 = randi([0 1],1,N);  

            x1 = 2*data1 - 1;
            x2 = 2*data2 - 1;

            C_noma_1 = log2(1 + pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); 
            C_noma_12 = log2(1 + pt(u)*a1.*g2./(pt(u)*a2.*g2+no));
            C_noma_2 = log2(1 + pt(u)*a2.*g2/no); 

            C_noma_sum(nearIdx, farIdx, u) = mean(C_noma_1 + C_noma_2); 

            C_oma_1 = (1/2)*log2(1 + pt(u)*g1/no);   
            C_oma_12 = (1/2)*log2(1 + pt(u)*g1/no);
            C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no);   

            C_oma_sum(nearIdx, farIdx, u) = mean(C_oma_1 + C_oma_2);  

            R1_av = mean(C_noma_1);
            R2_av = mean(C_noma_2);       

            SINR_noma_1(nearIdx, farIdx, u) = mean(pt(u)*a1.*g1./(pt(u)*a2.*g1+no)); 
            SINR_noma_2(nearIdx, farIdx,u) = mean(pt(u)*a2.*g2/no);

            SINR_oma_1(nearIdx, farIdx, u) = mean(pt(u)*g1/no);
            SINR_oma_2(nearIdx, farIdx, u) = mean(pt(u)*g2/no);

            for k = 1:N
                if C_noma_1(k) < rate1
                    pn1(nearIdx, farIdx, u) = pn1(nearIdx, farIdx, u)+1;
                end
                if (C_noma_12(k) < rate1)||(C_noma_2(k) < rate2)
                    pn2(nearIdx, farIdx, u) = pn2(nearIdx, farIdx, u)+1;
                end
            end

            poutNoma1(nearIdx, farIdx, u) = pn1(nearIdx, farIdx,u)/N;
            poutNoma2(nearIdx, farIdx, u) = pn2(nearIdx, farIdx,u)/N;

            for ko = 1:N
                if C_oma_1(ko) < rate1
                    po1(nearIdx, farIdx, u) = po1(nearIdx, farIdx, u)+1;
                end
                if (C_oma_2(ko) < rate2)
                    po2(nearIdx, farIdx, u) = po2(nearIdx, farIdx, u)+1;
                end
            end

            poutoma1(nearIdx, farIdx,u) = po1(nearIdx, farIdx, u)/N;
            poutoma2(nearIdx, farIdx,u) = po2(nearIdx, farIdx, u)/N;

            SNR = Pt(u) - No;
             % Update the dataset table entries
            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, d_near, d_far, Pt(u), SNR, C_oma_sum(nearIdx, farIdx, u), C_oma_sum(nearIdx, farIdx, u), poutoma1(nearIdx, farIdx, u), poutoma2(nearIdx, farIdx, u), 'OMA'};
            datasetTable = [datasetTable; newRow];
            
            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, d_near, d_far, Pt(u), SNR, C_noma_sum(nearIdx, farIdx, u), C_noma_sum(nearIdx, farIdx, u), poutNoma1(nearIdx, farIdx, u), poutNoma2(nearIdx, farIdx, u), 'NOMA'};
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

writetable(datasetTable, 'outage_data.csv');

% Display the metrics for OMA and NOMA
disp('Metrics for OMA:');
disp(['Spectral Efficiency: ', num2str(spectral_efficiency_OMA)]);
disp(['Data Throughput: ', num2str(data_throughput_OMA)]);

disp('Metrics for NOMA:');
disp(['Spectral Efficiency: ', num2str(spectral_efficiency_NOMA)]);
disp(['Data Throughput: ', num2str(data_throughput_NOMA)]);

