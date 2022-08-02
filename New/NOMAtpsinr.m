clc; 
clear variables; 
close all;

%SNR range
Pt = -114:5:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

N = 10^4;
dgap = [0.1,1,100];%in meters

for i = 1:3 %dgap values loop 
d3 = 2;      %near user %Distance of users
d2 = d3+dgap(i);%far user
eta = 4;                    %Path loss exponent

%Rayleigh fading coefficients of both users
h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h3 = (sqrt(d3^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Channel gains
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;

BW = 10^9;
No = -174 + 10*log10(BW);
no = (10^-3)*db2pow(No);

%Power allocation coefficients
a1 = 0.75; 
a2 = 0.1825; 
a3 = 1-(a1+a2);

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));

for u = 1:length(pt)
    
    %NOMA capacity calculation
    C_noma_2 = log2(1 + pt(u)*a3.*g3./(pt(u)*a2.*g2+no)); %far user                  
    C_noma_3 = log2(1 + pt(u)*a2.*g2/no); %near user  
    
    
    %gamma_far(u) = mean(C_noma_2);
    C_noma_sum(u) = mean(C_noma_2 + C_noma_3);  %Sum capacity of NOMA
    
    %OMA capacity calculation
    C_oma_2 = (1/3)*log2(1 + pt(u)*g2/no);    %User 2
    C_oma_3 = (1/3)*log2(1 + pt(u)*g3/no);    %User 3
    
    C_oma_sum(u) = mean(C_oma_2 + C_oma_3); %Sum capacity of OMA
    
end

SNR = Pt - No;
figure (1);
C = {'o' ; 'k' ; 'b'} ;
h = plot(SNR,C_noma_sum,'color',rand(1,3),'linewidth',2); hold on; grid on;
u = plot(SNR,C_oma_sum,'--','color',rand(1,3),'linewidth',2);
xlabel('SNR (dB)');
ylabel('Achievable sum rate (bps/Hz)');
legend('NOMA - 0.1m', 'NOMA - 1m','NOMA - 100m','OMA - 0.1m', 'OMA - 1m','OMA - 100m','Orientation','vertical');
%legend2 = ('NOMA - 0.1m', 'NOMA - 1m','NOMA - 100m','OMA - 0.1m', 'OMA - 1m','OMA - 100m');
title('Capacity of NOMA');
ylim([0 max(C_noma_sum)+1]);
ylim([0 max(C_oma_sum)+1]);
hold on ;
end



