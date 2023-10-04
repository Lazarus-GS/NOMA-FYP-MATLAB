clc;
close all;
clear all;

mT = 2;
mR = 3;
N = 10000; %number of bits
SNRdB = [-15:30];%SNR range in dB
SNR = 10.^(SNRdB/10);%SNR in linear scale
C_SISO = zeros(1,length(SNR));
C_SIMO = zeros(1,length(SNR));
C_MISO = zeros(1,length(SNR));
C_MIMO = zeros(1,length(SNR));
p_simo = zeros(1,length(SNR));
p_miso = zeros(1,length(SNR));

rateth = 1;
ITER = 10;%number of trials

p_simo_iter  = zeros(ITER,length(SNR));
p_miso_iter  = zeros(ITER,length(SNR));
p_siso_iter  = zeros(ITER,length(SNR));
p_mimo_iter  = zeros(ITER,length(SNR));

for i = 1:ITER
count = 0;
count2 = 0;
count_siso = 0;
count_mimo = 0;
p_simo = zeros(1,length(SNR));
p_miso = zeros(1,length(SNR));
for K = 1:length(SNR)
   for b = 1:N  %over each bit    
        h_SISO = (randn +1i*randn)/sqrt(2);
        h_SIMO = (randn(mR,1)+1i*randn(mR,1))/sqrt(2);
        h_MISO = (randn(1,mT)+1i*randn(1,mT))/sqrt(2);
        h_MIMO = (randn(mR,mT)+1i*randn(mR,mT))/sqrt(2);

        C_SISO(K) = C_SISO(K) + log2(1+ SNR(K)*norm(h_SISO)^2);
        C_SIMO(K) = C_SIMO(K) + log2(1+ SNR(K)*norm(h_SIMO)^2);
        C_MISO(K) = C_MISO(K) + log2(1+ SNR(K)*norm(h_MISO)^2/mT);
        C_MIMO(K) = C_MIMO(K) + log2(abs(det(eye(mR)+SNR(K)*h_MIMO*h_MIMO'/mT)));

        %outage analysis
        s(b) = log2(1+ SNR(K)*norm(h_SIMO)^2);
        g(b) = log2(1+ SNR(K)*norm(h_MISO)^2/mT);
        siso(b) = log2(1+ SNR(K)*norm(h_SISO)^2);
        mimo(b) = log2(abs(det(eye(mR)+SNR(K)*h_MIMO*h_MIMO'/mT)));
       
        if(s(b)<rateth)
            count = count+1;
        end    
        if(g(b)<rateth)
            count2 = count2+1;
        end
        if(siso(b)<rateth)
            count_siso = count_siso+1;
        end    
        if(mimo(b)<rateth)
            count_mimo = count_mimo+1;
        end
   end  
   
   p_simo(K)  = count/N;  
   p_miso(K)  = count2/N;
   p_siso(K)  = count_siso/N;  
   p_mimo(K)  = count_mimo/N;
   count = 0;
   count2 = 0;
   count_siso = 0;
   count_mimo = 0;
   
end

%average capacity
C_SISO = (C_SISO/N);
C_SIMO = (C_SIMO/N);
C_MISO = (C_MISO/N);
C_MIMO = (C_MIMO/N);

p_simo_iter(ITER,:) = p_simo_iter(ITER,:) + p_simo;
p_miso_iter(ITER,:) = p_miso_iter(ITER,:) + p_miso;
p_siso_iter(ITER,:) = p_siso_iter(ITER,:) + p_siso;
p_mimo_iter(ITER,:) = p_mimo_iter(ITER,:) + p_mimo;
end

%average capacity over ITER iterations
C_SISO = (C_SISO/ITER);
C_SIMO = (C_SIMO/ITER);
C_MISO = (C_MISO/ITER);
C_MIMO = (C_MIMO/ITER);

figure(1)
grid on;
plot(SNRdB,C_SISO,'r',SNRdB,C_SIMO,'b',SNRdB,C_MISO,'m',SNRdB,C_MIMO,'k')
legend('SISO','SIMO','MISO','MIMO');
xlabel('SNR in dB');
ylabel('Capacity (b/s/Hz)');
title('Capacity Vs. SNR');

figure(2)
grid on;
plot(SNRdB,mean(p_simo_iter),'r',SNRdB,mean(p_miso_iter),...
    'b--',SNRdB,mean(p_siso_iter),'g--',SNRdB,mean(p_mimo_iter),'k-');
legend('SIMO','MISO','SISO','MIMO');
xlabel('SNR in dB');
ylabel('Outage (b/s/Hz)');
title('Outage Vs. SNR');