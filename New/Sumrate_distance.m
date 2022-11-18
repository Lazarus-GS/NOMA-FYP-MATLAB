clc; 
clear variables; 
close all;

%transmit power range
Pt = -114:5:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

N = 10^4;
%dgap = 1; %in meters

d1 = 302; d2 = 202; d3 = 102; d4 = 2;	%Distances
a1 = 0.49; a2 = 0.33; a3 = 0.17; a4 = 0.01;	%Power allocation coefficients

eta = 4;	%Path loss exponent

%Rayleigh fading coefficients of both users
h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h3 = (sqrt(d3^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h4 = (sqrt(d4^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;

BW = 10^9; %bandwidth
No = -174 + 10*log10(BW);
no = (10^-3)*db2pow(No);

dbstop if error

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));
C_noma_sum = zeros(1,length(pt));

for u = 1:length(pt)
    
    %NOMA capacity calculation
    C_noma_1 = log2(1 + pt(u)*a1.*g1./(pt(u)*a2.*g1 + pt(u)*a3.*g1 + pt(u)*a4.*g1 + no)); %User 1
    C_noma_2 = log2(1 + pt(u)*a2.*g2./(pt(u)*a3.*g2 + pt(u)*a4.*g2 + no));                   %User 2
    C_noma_3 = log2(1 + pt(u)*a3.*g3/(pt(u)*a4.*g3 + no));
    C_noma_4 = log2(1 + pt(u)*a4.*g4/no);
    
    %gamma_far(u) = mean(C_noma_2);
    %C_noma_sum(u) = mean(C_noma_1 + C_noma_2 + C_noma_3 + C_noma_4);  %Sum capacity of NOMA
    
    %OMA capacity calculation
    C_oma_1 = (1/4)*log2(1 + pt(u)*g1/no);    %User 1
    C_oma_2 = (1/4)*log2(1 + pt(u)*g2/no);    %User 2
    C_oma_3 = (1/4)*log2(1 + pt(u)*g3/no);    %User 3
    C_oma_4 = (1/4)*log2(1 + pt(u)*g4/no);
    
    %C_oma_sum(u) = mean(C_oma_1 + C_oma_2 + C_oma_3 + C_oma_4); %Sum capacity of OMA
end    


SNR = Pt - No;
figure (1);

%plot(SNR,C_noma_1,'color',rand(1,3),'linewidth',2); hold on; grid on;
%plot(SNR,C_noma_2,'color',rand(1,3),'linewidth',2); hold on; grid on;
%plot(SNR,C_noma_3,'color',rand(1,3),'linewidth',2); hold on; grid on;
%plot(SNR,C_noma_4,'color',rand(1,3),'linewidth',2); hold on; grid on;
plot(SNR,C_oma_1,'--','color',rand(1,3),'linewidth',2); hold on; grid on;
plot(SNR,C_oma_2,'--','color',rand(1,3),'linewidth',2); hold on; grid on;
plot(SNR,C_oma_3,'--','color',rand(1,3),'linewidth',2); hold on; grid on;
plot(SNR,C_oma_4,'--','color',rand(1,3),'linewidth',2); hold on; grid on;
xlabel('SNR (dB)');
ylabel('Achievable sum rate (bps/Hz)');

legend(handlevec,'NOMA - 0.1m', 'NOMA - 1m','NOMA - 100m','OMA - 0.1m', 'OMA - 1m','OMA - 100m');
title('Capacity of NOMA');
ylim([0 max(C_noma)+1]);
ylim([0 max(C_oma)+1]);
hold on ;







