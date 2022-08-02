%ref: https://ecewireless.blogspot.com/2020/09/noma-vs-oma...
%-capacity-comparison.html

clc; clear variables; close all;

%SNR range
Pt = -114:5:-54;	%in dB
pt = db2pow(Pt);	%in linear scale

N = 10^4;
dgap = [1,100];

for i = 1:2
d3 = 2;      %Distance of users
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
a1 = 0.75; a2 = 0.1825; a3 = 1-(a1+a2);

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));

for u = 1:length(pt)
    
    %NOMA capacity calculation
    C_noma_2 = log2(1 + pt(u)*a3.*g3./(pt(u)*a2.*g2+no));                   %User 2
    C_noma_3 = log2(1 + pt(u)*a2.*g2/no);                   %User 3
    gamma_far(u) = mean(C_noma_2);
    C_noma_sum(u) = mean(C_noma_2 + C_noma_3);  %Sum capacity of NOMA
    
end

SNR = Pt - No
figure (1);
plot(SNR,C_noma_sum,'linewidth',2); hold on; grid on;
xlabel('SNR (dB)');
ylabel('Achievable sum rate (bps/Hz)');
legend('NOMA');
title('Capacity of NOMA');
ylim([0 max(C_noma_sum)+1]);
hold on ;
end