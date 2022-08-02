% ----------------------------------------------------------------------- %
% -------------- Effect of Power Allocation Coefficients ---------------- %
% ------------------------- Author : T.Benn Roshnan --------------------- %
% ----------------------------------------------------------------------- %

clc;
clear all;
close all;

P = 1 ; % Total Power
a = 0.05:0.05:0.45 ; % Power allocation Factor   %% ------------??????????? check here
E1 = a.*P ; % Near user Power (user-1)
E2 = (1-a).*P ; % Far user Power (user-2)

N = 1e6 ; % No. of (bits)''SYMBOLS'' transmitted 
SNR_dB = 10:10:40;% SNR in decibels
SNR_lin = 10.^(SNR_dB./10) ; % Conversion of SNR to linear form  

ha = 0.5012; % for the channel allocation process => reflecting (- 3 dB)

BER_FU = zeros(length(SNR_lin),length(a));
BER_NU = zeros(length(SNR_lin),length(a));
BER_FU_th = zeros(length(SNR_lin),length(a));
BER_NU_th = zeros(length(SNR_lin),length(a));

m_NU = randi([0 1],1,N); % generating the NU bits 
m_FU = randi([0 1],1,N); % generating the FU bits

h = (randn(1,N) + 1j*randn(1,N)); 
h_1 = sqrt(1/2)*h; % Channel Coefficients for near user ?????????
h_2 = sqrt(ha)*sqrt(1/2)*h; % Channel Coefficients for far user

noise_1 =sqrt(1/2)*(randn(1,N) + 1j*randn(1,N));
noise_2 =sqrt(1/2)*(randn(1,N) + 1j*randn(1,N));

sd = 1./sqrt(SNR_lin);

for m=1:length(SNR_lin)
    for alpha=1:length(a)
    
        m_FU_mod = sqrt(E2(alpha))*BPSK_mod(m_FU,N); % modulating FU bits as BPSK symbols
        m_NU_mod = sqrt(E1(alpha))*BPSK_mod(m_NU,N); % modulating NU bits as QPSK symbols

        m_SC = m_NU_mod + m_FU_mod; % Superposition Coding
    
        rec_NU = m_SC .* h_1 + sd(m) *noise_1; % received message NU at receiver
        rec_FU = m_SC .* h_2 + sd(m) *noise_2; % received message FU at receiver
    
        det_user2 = rec_FU./h_2; % detected symbols at user 2
        det_user1 = rec_NU./h_1; % detected symbols at user 0
    
        det_m_FU = BPSK_demod (det_user2,N); % detected symbols at user 2 decoded as BPSK symbols

    
    %-----------------------------------------------------------------------------------------------------------%
    %------------------------------------ SIC Process at Near User ---------------------------------------------%
    %-----------------------------------------------------------------------------------------------------------%
        det_m_FU_at_U1 = BPSK_demod(det_user1,N); % detected symbols at user 1 decoded as BPSK symbols and fU symbols get detected  
        det_symb_U1_NU = det_user1 - (sqrt(E2(alpha))*BPSK_mod(det_m_FU_at_U1,N)); % near user symbols get detected as QPSK symbols
        det_m_NU = BPSK_demod(det_symb_U1_NU,N); % NU bits get retrieved
    %-----------------------------------------------------------------------------------------------------------%
    
    
    %-----------------------------------------------------------------------------------------------------------%
    %----------------------------------- Bit/Symbol Error rate calculation -------------------------------------%
    %-----------------------------------------------------------------------------------------------------------%
        
        BER_FU(m,alpha) = sum(det_m_FU ~= m_FU)/(N);
        
        BER_NU(m,alpha) = sum(det_m_NU ~= m_NU)/(N); 
    
    %-----------------------------------------------------------------------------------------------------------%
    end
end


%--------------------------------------------------------------------------------------------------------------%
%----------------------------------------Theoreticalexpressions------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------%

for k=1:length(SNR_dB)
    
    No = (sd(k)^2); 
    
    avgamma_a = ((sqrt(2.*E2) + sqrt(2.*E1)).^2)*ha./(No);
    avgamma_b = ((sqrt(2.*E2) - sqrt(2.*E1)).^2)*ha./(No); 
    
    avgamma_c = 2.*E1.*1./(No); 
    avgamma_d = ((sqrt(2.*E2) + sqrt(2.*E1)).^2)*1./(No); 
    avgamma_e = ((2*sqrt(2.*E2) + sqrt(2.*E1)).^2)*1./(No); 
    avgamma_f = ((2*sqrt(2.*E2) - sqrt(2.*E1)).^2)*1./(No); 
    avgamma_g = ((sqrt(2.*E2) - sqrt(2.*E1)).^2).*1./(No);

    BER_FU_th(k,:) = 0.25 .* ( 2 - sqrt(avgamma_a./(2+avgamma_a)) - sqrt(avgamma_b./(2+avgamma_b))); % BER of the Far user 
    BER_NU_th(k,:) = (0.5 .* ( 1 - sqrt(avgamma_c./(2+avgamma_c))))  +  0.25 .* (  sqrt(avgamma_d./(2+avgamma_d)) + sqrt(avgamma_f./(2+avgamma_f)) - sqrt(avgamma_e./(2+avgamma_e)) - sqrt(avgamma_g./(2+avgamma_g)) ); 

end

%--------------------------------------------------------------------------------------------------------------%



%----------------------------------------------------------------------------------------------------------------%
%------------------------------------------------- Plotting the Graph -------------------------------------------%
%----------------------------------------------------------------------------------------------------------------%

% for k = 1:length(SNR_lin)
%     semilogy (a,BER_FU_th(k,:),a,SER_NU_th(k,:));grid on;hold on
% end    

semilogy(a,BER_FU,'p-m');grid on;hold on 
semilogy(a,BER_NU,'d-r');hold on;
semilogy(a,BER_FU_th,'k'); hold on;
semilogy(a,BER_NU_th,'k'); hold on;
xlabel('Power Allocation Coefficient (a)');
ylabel('BER');
legend('simulated DL-NOMA for FU','','','','simulated DL-NOMA for NU','','','','analytical','','',''); 
title ('Effect of the Power Allocation Coefficient for NOMA Down Link and Up Link ');
%----------------------------------------------------------------------------------------------------------------%