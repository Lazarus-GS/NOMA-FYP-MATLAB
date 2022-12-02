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

C_noma = zeros(1,length(pt));
C_oma = zeros(1,length(pt));
SINR_noma = zeros(1,length(pt));
SINR_oma = zeros(1,length(pt));

p = length(Pt);
pn1 = zeros(1,length(Pt));
pn2 = zeros(1,length(Pt));
po1 = zeros(1,length(Pt));
po2 = zeros(1,length(Pt));

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
    C_oma_12 = (1/2)*log2(1 + pt(u)*g1/no);
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
            pn1(u) = pn1(u)+1;
        end
        if (C_noma_12(k) < rate1)||(C_noma_2(k) < rate2)
            pn2(u) = pn2(u)+1;
        end
    end
    
    poutNoma1 = pn1/N;
    poutNoma2 = pn2/N;
    
    
    %OMA
    
    C_SISO = zeros(1,length(SNR));
    

    rateth = 1;
    ITER = 10;%number of trials

    p_siso_iter  = zeros(ITER,length(SNR));
    

    for i = 1:ITER
    count = 0;
    count2 = 0;
    count_siso = 0;
    
    
    for K = 1:length(SNR)
       for b = 1:N  %over each bit    
            h_SISO = (randn +1i*randn)/sqrt(2);

            C_SISO(K) = C_SISO(K) + log2(1+ SNR(K)*norm(h_SISO)^2);
            

            %outage analysis
            
            siso(b) = log2(1+ SNR(K)*norm(h_SISO)^2);
            
            
            if(siso(b)<rateth)
                count_siso = count_siso+1;
            end    
       end  

       p_siso(K)  = count_siso/N;  
       
       count = 0;
       count2 = 0;
       count_siso = 0;
       

    end

    %average capacity
    C_SISO = (C_SISO/N);

    p_siso_iter(ITER,:) = p_siso_iter(ITER,:) + p_siso;
    end

    %average capacity over ITER iterations
    C_SISO = (C_SISO/ITER);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%BER
%NOMA

    %Do superposition coding
    x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);
    %Received signals
    y1 = h1.*x + w1;
    y2 = h2.*x + w2;
    
    %Equalize 
    eq1 = y1./h1;
    eq2 = y2./h2;
    
    %AT USER 1--------------------
    %Direct decoding of x1 from y1
    x1_hat = zeros(1,N);
    x1_hat(eq1>0) = 1;
    
    %Compare decoded x1_hat with data1 to estimate BER
    ber1(u) = biterr(data1,x1_hat)/N;
    
    %----------------------------------
    
    %AT USER 2-------------------------
    %Direct decoding of x1 from y2
    x12_hat = ones(1,N);
    x12_hat(eq2<0) = -1;
    
    y2_dash = eq2 - sqrt(a1*pt(u))*x12_hat;
    x2_hat = zeros(1,N);
    x2_hat(real(y2_dash)>0) = 1;
    
    ber2(u) = biterr(x2_hat, data2)/N;
    %-----------------------------------

    
%OMA


end





end

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