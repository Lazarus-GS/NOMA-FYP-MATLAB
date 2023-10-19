clc; 
clear variables; 
close all;

columnNames = {'NearUserPos', 'FarUserPos', 'FarUserDist', 'TransmitPower', 'SNR', 'SumRateFarUser', 'SumRateNearUser', 'BERFarUser', 'BERNearUser', 'OutageProbFarUser', 'OutageProbNearUser', 'MAScheme'};
datasetTable = table('Size', [0 12], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string'}, 'VariableNames', columnNames);

nearUserPositions = {'b2', 'b3', 'c2', 'c3'};
farUserPositions = {'a1', 'a2', 'a3', 'a4', 'b1', 'b4', 'c1', 'c4', 'd1', 'd2', 'd3', 'd4'};

seed = 10;
rng(seed);

Pt = -114:0.5:-54;	
pt = db2pow(Pt);	

rate1 = 1; rate2 = 2;       

N = 10^4; 

un = 7.071;
uf = [21.2132, 15.811, 15.811, 21.2132, 15.811, 15.811, 15.811, 15.811, 21.2132, 15.811, 15.811, 21.2132];

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

            simBer1(farIdx,u) = biterr(data1,ipHat1)/N;
            simBer2(farIdx,u) = biterr(data2,ipHat2)/N;

            SNR = Pt(u) - No;

            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, uf(farIdx), Pt(u), SNR, C_oma_sum(farIdx,u), C_oma_sum(nearIdx,u), ber1(farIdx,u), ber2(farIdx,u), poutoma1(farIdx,u), poutoma2(farIdx,u), 'OMA'};
            datasetTable = [datasetTable; newRow];
            
            newRow = {nearUserPositions{nearIdx}, farUserPositions{farIdx}, uf(farIdx), Pt(u), SNR, C_noma_sum(farIdx,u), C_noma_sum(nearIdx,u), simBer1(farIdx,u), simBer2(farIdx,u),  poutNoma1(farIdx,u), poutNoma2(farIdx,u), 'NOMA'};
            datasetTable = [datasetTable; newRow];
        end
    end
end

writetable(datasetTable, 'simulation_dataset1.csv');
