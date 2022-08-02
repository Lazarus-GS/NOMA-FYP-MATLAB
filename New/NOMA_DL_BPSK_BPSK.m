% ----------------------------------------------------------------------- %
% ------------------- BER Comaprison for NOMA Down Link ----------------- %
% ---------------------------- (NU:BPSK,FU:BPSK) ------------------------ %
% ------------------------- Author : T.Benn Roshnan --------------------- %
% ----------------------------------------------------------------------- %

tic;

clc;
clear all;
close all;

N = 1e5 ; % No. of (bits)''SYMBOLS'' transmitted 
SNR_dB = 0:2:40; % SNR in decibels
SNR_lin = 10.^(SNR_dB./10) ; % Conversion of SNR to linear form  

P = 1 ; % Total Power
a = 0.4 ; % Power allocation Factor   %% ------------??????????? check here

E1 = a*P ; % Near user Power (user-1)
E2 = (1-a)*P ; % Far user Power (user-2)

ha_dB = -6;
ha = 10^(ha_dB/10); % for the channel allocation process => reflecting (- 3 dB)

BER_FU = zeros(1,length(SNR_lin));
BER_NU = zeros(1,length(SNR_lin));

m_NU = randi([0 1],1,N); % generating the NU bits 
m_FU = randi([0 1],1,N); % generating the FU bits

m_FU_mod = sqrt(E2)*BPSK_mod(m_FU,N); % modulating FU bits as BPSK symbols
m_NU_mod = sqrt(E1)*BPSK_mod(m_NU,N); % modulating NU bits as QPSK symbols

m_SC = m_NU_mod + m_FU_mod; % Superposition Coding

h = (randn(1,N) + 1j*randn(1,N)); 
h_1 = sqrt(1/2)*h; % Channel Coefficients for near user 
h_2 = sqrt(ha)*sqrt(1/2)*h; % Channel Coefficients for far user

noise_1 = sqrt(1/2)*(randn(1,N) + 1j*randn(1,N));
noise_2 = sqrt(1/2)*(randn(1,N) + 1j*randn(1,N));

sd = sqrt(P./SNR_lin);

for k=1:length(SNR_lin)

    rec_NU = m_SC .* h_1 + sd(k) *noise_1; % received message NU at receiver
    rec_FU = m_SC .* h_2 + sd(k) *noise_2; % received message FU at receiver
    
    det_user2 = rec_FU./h_2; % detected symbols at user 2
    det_user1 = rec_NU./h_1; % detected symbols at user 0
    
    det_m_FU = BPSK_demod (det_user2,N); % detected symbols at user 2 decoded as BPSK symbols

    
    %-----------------------------------------------------------------------------------------------------------%
    %------------------------------------ SIC Process at Near User ---------------------------------------------%
    %-----------------------------------------------------------------------------------------------------------%
    det_m_FU_at_U1 = BPSK_demod(det_user1,N); % detected symbols at user 1 decoded as BPSK symbols and fU symbols get detected  
    det_symb_U1_NU = det_user1 - (sqrt(E2)*BPSK_mod(det_m_FU_at_U1,N)); % near user symbols get detected as QPSK symbols
    det_m_NU = BPSK_demod(det_symb_U1_NU,N); % NU bits get retrieved
    %-----------------------------------------------------------------------------------------------------------%
    
    
    %-----------------------------------------------------------------------------------------------------------%
    %----------------------------------- Bit/Symbol Error rate calculation -------------------------------------%
    %-----------------------------------------------------------------------------------------------------------%
    
    BER_FU(k) = sum(det_m_FU ~= m_FU)/(N);
    BER_NU(k) = sum(det_m_NU ~= m_NU)/(N);

    %-----------------------------------------------------------------------------------------------------------%
    
    
end

    No = (sd.^2); 
    avgamma_a = ((sqrt(2*E2) + sqrt(2*E1))^2)*ha./(No);
    avgamma_b = ((sqrt(2*E2) - sqrt(2*E1))^2)*ha./(No); 
    
    avgamma_c = 2*E1*1./(No); 
    avgamma_d = ((sqrt(2*E2) + sqrt(2*E1))^2)*1./(No); 
    avgamma_e = ((2*sqrt(2*E2) + sqrt(2*E1))^2)*1./(No); 
    avgamma_f = ((2*sqrt(2*E2) - sqrt(2*E1))^2)*1./(No); 
    avgamma_g = ((sqrt(2*E2) - sqrt(2*E1))^2)*1./(No); 
    


%--------------------------------------------------------------------------------------------------------------%
%-------------------------------------------Checking Process---------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------%

scatterplot(det_symb_U1_NU); grid on; title('det symb U1 NU');
% scatterplot(m_NU_mod) ; grid on; title('m NU mod');
% scatterplot(m_FU_mod) ; grid on; title('m FU mod');
% scatterplot(m_SC) ; grid on; title('m SC');
% scatterplot(rec_NU./h_1); grid on; title('rec NU');
% scatterplot(rec_FU./h_2); grid on; title('rec FU');
% det_symb_U1_NU
% scatterplot(det_user1); grid on;title('det user1'); 
% scatterplot(det_user2); grid on;title('det user2');
% scatterplot(det_m_FU); grid on; title('det m FU');
% scatterplot(det_m_FU_at_U1); grid on; title('det m FU at U1');
scatterplot(det_symb_U1_NU); grid on; title('det symb U1 NU');

%--------------------------------------------------------------------------------------------------------------%


%--------------------------------------------------------------------------------------------------------------%
%----------------------------------------Theoreticalexpressions------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------%
 BER_FU_th = 0.25 .* ( 2 - sqrt(avgamma_a./(2+avgamma_a)) - sqrt(avgamma_b./(2+avgamma_b))); % BER of the Far user 
 BER_NU_th = (0.5 .* ( 1 - sqrt(avgamma_c./(2+avgamma_c))))  +  0.25 .* (  sqrt(avgamma_d./(2+avgamma_d)) + sqrt(avgamma_f./(2+avgamma_f)) - sqrt(avgamma_e./(2+avgamma_e)) - sqrt(avgamma_g./(2+avgamma_g)) ); 

%--------------------------------- Theoretical Expression for OMA ---------------------------------------------%
 BER_OMA_th = 0.5.*(1-sqrt(SNR_lin./(2+SNR_lin))) ;
%--------------------------------------------------------------------------------------------------------------%


%----------------------------------------------------------------------------------------------------------------%
%------------------------------------------------- Plotting the Graph -------------------------------------------%
%----------------------------------------------------------------------------------------------------------------%

%  semilogy(SNR_dB,BER_FU,'*-',SNR_dB,BER_NU,'*-',SNR_dB,BER_FU_th,SNR_dB,BER_NU_th,SNR_dB,BER_OMA_th,'*-');grid on;
% semilogy(SNR_dB,BER_NU,'o-',SNR_dB,BER_NU_th,'*-');grid on;
% semilogy(SNR_dB,BER_FU,'o-',SNR_dB,BER_FU_th,'*-');grid on;
% semilogy(SNR_dB,BER_FU_th,'o-',SNR_dB,SER_NU_th,'<-');grid on;
% semilogy(SNR_dB,BER_FU,'p-',SNR_dB,BER_NU,'o-',SNR_dB,BER_OMA_th,'k-');grid on;


xlabel('Transmit SNR(dB)');
ylabel('BER');
legend('simulated DL-NOMA for FU','simulated DL-NOMA for NU','analytical DL-NOMA for FU','analytical DL-NOMA for NU','OMA'); 
title ('BER Comaprison for NOMA Down Link ')

%----------------------------------------------------------------------------------------------------------------%

toc
