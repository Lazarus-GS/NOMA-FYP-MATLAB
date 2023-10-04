clc
clear all;
close all;
n=10^6;%no_of_samples
i=randi([0,1],1,n);%generates random integers 0's and 1's
i1=2*i-1;%bpsk modulation  i.e, mapping 1's as 1 and 0's as -1

%rayleigh fading channel
a=randn(1,n);%generates samples of size 1xn which are gaussian distributed
b=randn(1,n);%generates samples of size 1xn which are gaussian distributed
rc=1/sqrt(2)*(sqrt(a.^2+b.^2));%rayleigh channel

for l=0:1:20
snr=10^((l/10));    % SNR values in absolute scale
sdev=sqrt(0.5/snr); % standard deviation of noise calculated from SNR
N=random('norm',0,sdev,[1,n]);% generation of noise sequence
yrc=rc.*i1+N;       %signal received through rayleigh and awgn channel


YR=(yrc>=0); %baseband detection from Rayleigh,AWGN channel
ErrorR=sum((xor(YR,i)));% no of errors in detected signal


ber_R(l+1)=ErrorR/n;%simulated BER for AWGN,rayleigh channel

berthR(l+1)=0.5*(1-sqrt(snr/(snr+1)));%theoretical bit error rate of
% rayleigh,awgn channel
p=((1-2*ber_R(l+1))^2)/(4*(ber_R(l+1)-(ber_R(l+1)^2)));

outage(l+1)=1-exp(-3.16/p);%simulated ber computed from rayleigh channel
outageT(l+1)=1-exp(-3.16/snr);%theoretical ber computed from rayleigh
%channel

end
% scatterplot(yawgn)
% scatterplot(yrc)
figure
%comparison plot generation of theoretical and calculated SNR
q=0:1:20 ;

semilogy(q,berthR(q+1),'+r');
hold on
semilogy(q,ber_R(q+1),'g-');
axis([0 20 10^-5 1]);
%title('BER PERFORMANCE OF BPSK MODULATION SCHEME IN AWGN,RAYLEIGH CHANNEL')
xlabel('SNR')
ylabel('BIT ERROR RATE')
%figure
% plot(q,outage(q+1),'-');
% hold on
% plot(q,outageT(q+1),'+g');
%title('OUTAGE PERFORMANCE OF BPSK MODULATION IN AWGN,RAYLEIGH CHANNEL')
xlabel('SNR')
ylabel('OUTAGE PROBABILITY')

    
    
