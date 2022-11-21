clc; 
clear variables; 
close all;

%transmit power range
Pt = -114:5:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

rate1 = 1; rate2 = 2;       %Target rate of users in bps/Hz

N = 10^4;
dgap = [0.1,1,100]; %in meters

for i = 1:3 %dgap values loop 
d2 = 2;      %near user %Distance of users
d1 = d2+dgap(i); %far user
eta = 4;                    %Path loss exponent

%Rayleigh fading coefficients of both users
h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); %far
h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2); %near

%Channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

BW = 10^9; %bandwidth
No = -174 + 10*log10(BW);
no = (10^-3)*db2pow(No);

%Pt = 0:2:40;                %Transmit power in dBm
%pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
%BW = 10^6;                  %System bandwidth
%No = -174 + 10*log10(BW);   %Noise power (dBm)
%no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

%Power allocation coefficients
a1 = 0.75; 
a2 = 1-a1; 
%a3 = 1-(a1+a2);

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));
SINR_noma = zeros(1,length(pt));
SINR_oma = zeros(1,length(pt));

p = length(Pt);
p1 = zeros(1,length(Pt));
p2 = zeros(1,length(Pt));

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
    C_noma_sum(u) = mean(C_noma_1 + C_noma_2);  %Sum capacity of NOMA

    %OMA capacity calculation
    C_oma_1 = (1/2)*log2(1 + pt(u)*g1/no);    %User 1
    C_oma_2 = (1/2)*log2(1 + pt(u)*g2/no);    %User 2
    
    C_oma_sum(u) = mean(C_oma_1 + C_oma_2); %Sum capacity of OMA
    
    %Average of achievable rates
    R1_av(u) = mean(C_noma_1);
    R2_av(u) = mean(C_noma_2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    %Recieved SINR
    %NOMA recieved SINR
    
    SINR_noma_1 = pt(u)*a1.*g1./(pt(u)*a2.*g1+no); %far user 
    SINR_noma_2 = pt(u)*a2.*g2/no; %near user
    
    %OMA Recieved SNR
    
    SINR_oma_1 = pt(u)*g1/no;
    SINR_oma_2 = pt(u)*g2/no;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    
    %Outage Probability
    %NOMA 
    
    for k = 1:N
        if C_noma_1(k) < rate1
            p1(u) = p1(u)+1;
        end
        if (C_noma_12(k) < rate1)||(C_noma_2(k) < rate2)
            p2(u) = p2(u)+1;
        end
    end
    
end

SNR = Pt - No;
figure (1);

plot(SNR,C_noma_sum,'color',rand(1,3),'linewidth',2); hold on; grid on;
plot(SNR,C_oma_sum,'--','color',rand(1,3),'linewidth',2);
xlabel('SNR (dB)');
ylabel('Achievable sum rate (bps/Hz)');


legend('NOMA - 0.1m', 'NOMA - 1m','NOMA - 100m','OMA - 0.1m', 'OMA - 1m','OMA - 100m');
title('Capacity of NOMA');
ylim([0 max(C_noma_sum)+1]);
ylim([0 max(C_oma_sum)+1]);
hold on ;

end