clc
clear all;
close all;
n=10^6;%no_of_samples
i=randi([0,1],1,n);%generates random integers 0's and 1's
i1=2*i-1;%bpsk modulation  i.e, mapping 1's as 1 and 0's as -1
%scatterplot(i1);%scatterplot of bpsk modulation scheme
%rayleigh fading channel
a=randn(1,n);%generates samples of size 1xn which are gaussian distributed
b=randn(1,n);%generates samples of size 1xn which are gaussian distributed
rc=1/sqrt(2)*(sqrt(a.^2+b.^2));%rayleigh channel
for l=0:1:20
snr=10^((l/10));    % SNR values in absolute scale
sdev=sqrt(0.5/snr); % standard deviation of noise calculated from SNR
N=random('norm',0,sdev,[1,n]);% generation of noise sequence
yrc=rc.*i1+N;       %signal received through rayleigh and awgn channel
yawgn=i1+N;         %signal received through awgn channel
Yb=(yawgn>=0);% baseband signal detection from awgn channel
YR=(yrc>=0); %baseband detection from Rayleigh,AWGN channel
ErrorR=sum((xor(YR,i)));% no of errors in detected signal
ErrorA=sum((xor(Yb,i)));

ber_A(l+1)=ErrorA/n;% simulated BER for awgn channel
ber_R(l+1)=ErrorR/n;%simulated BER for AWGN,rayleigh channel

berthR(l+1)=0.5*(1-sqrt(snr/(snr+1)));%theoretical bit error rate of
% rayleigh,awgn channel
p=((1-2*ber_R(l+1))^2)/(4*(ber_R(l+1)-(ber_R(l+1)^2)));

outage(l+1)=1-exp(-3.16/p);%simulated ber computed from rayleigh channel
outageT(l+1)=1-exp(-3.16/snr);%theoretical ber computed from rayleigh
%channel
berthA(l+1)=0.5*erfc(sqrt(2*snr));%theoretical ber computed for awgn channel
end
% scatterplot(yawgn)
% scatterplot(yrc)
figure
%comparison plot generation of theoretical and calculated SNR
q=0:1:20 ;

semilogy(q,ber_A(q+1),'-^b');
hold on
semilogy(q,berthA(q+1),'->y');
hold on
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
%RICIAN FADING
k1=10; %Rician factor
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
Nr2=randn(1,length(i1))*sigma+mean;
Ni2=randn(1,length(i1))*sigma;
%To generate the Rician Random Variable
No3=sqrt(Nr2.^2+Ni2.^2); %Rician fading coefficient
for k=0:1:20
    snrl=10^(k/10);%convert the SNR in dB value
    Np=1/snrl;%To generate the noise power
    sd=sqrt(Np/2);% standard deviation of guassian noise
    No=random('Normal',0,sd,1,length(i1)); %Generates Gaussian noise
    t1=i1.*No3+No; % s means transmitted signal...please take the value as u have taken
    z1=t1./No3;
 op1=(z1>0); % threshold detection
    Berr(k+1)=sum(xor(op1,i))/n; % observed BER
    BerTr(k+1)=.5*erfc(sqrt(k1*snrl/(k1+snrl)));% theoretical BER
    end;
    %figure;
    k=0:1:20;
    semilogy(k,Berr(k+1),'-*');
    hold on;
    semilogy(k,BerTr(k+1),'-<');
    % axis([0 10 10^-5 1]);
% %     plot(k,Berr(k),'-');
%     hold on;
%     plot(k,BerTr(k),'+g');
title('BER PERFORMANCE OF BPSK MODULATION SCHEME IN AWGN RAYLEIGH RICIAN FADING')
xlabel('SNR')
ylabel('BIT ERROR RATE')
hleg=legend('BER AWGN','BERth AWGN','BER rayleigh','BERth rayleigh','BER Rician','BERth Rician')
    
    
    
