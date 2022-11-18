clc; 
clear variables; 
close all;

%transmit power range
Pt = -114:5:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

N = 10^4;
dgap = [0.1,1,100]; %in meters

for i = 1:3 %dgap values loop 
d3 = 2;      %near user %Distance of users
d2 = d3+dgap(i); %far user
eta = 4;                    %Path loss exponent

%Rayleigh fading coefficients of both users
h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h3 = (sqrt(d3^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Channel gains
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;

BW = 10^9; %bandwidth
No = -174 + 10*log10(BW);
no = (10^-3)*db2pow(No);

%Pt = 0:2:40;                %Transmit power in dBm
%pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
%BW = 10^6;                  %System bandwidth
%No = -174 + 10*log10(BW);   %Noise power (dBm)
%no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

%Power allocation coefficients
a2 = 0.75; 
a3 = 1-a2; 
%a3 = 1-(a1+a2);

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));

for u = 1:2
    
    %NOMA capacity calculation
    C_noma_2 = log2(1 + pt(u)*a3.*g3./(pt(u)*a2.*g2+no)); %far user                  
    C_noma_3 = log2(1 + pt(u)*a2.*g2/no); %near user  
    
    
    %gamma_far(u) = mean(C_noma_2);
    C_noma_sum(u) = mean(C_noma_2 + C_noma_3);  %Sum capacity of NOMA
    
    %OMA capacity calculation
    C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no);    %User 2
    C_oma_3 = (1/2)*log2(1 + pt(u)*g3/no);    %User 3
    
    C_oma_sum(u) = mean(C_oma_2 + C_oma_3); %Sum capacity of OMA
    
    deltaSignal = abs(C_noma_sum - C_oma_sum);
    percentageDifference = deltaSignal ./ C_oma_sum; % Percent by element.
    meanPctDiff = mean(percentageDifference); % Average percentage over all elements.

    disp(percentageDifference)
    %disp(meanPctDiff*100)

    %percentError = 0.03696863;
    fprintf('Mean Percentage Difference: %0.2f%%', meanPctDiff*100);


    
end

SNR = Pt - No;

end


