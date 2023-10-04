%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share — to copy, distribute and transmit the work
% to Remix — to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating 4-PAM transmission and reception and compare the 
% simulated and theoretical symbol error probability

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 16 February 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% symbol error rate for 4-PAM modulation

clear
N = 10^5; % number of symbols
alpha4pam = [-3 -1 1 3]; % 4-PAM alphabets
Es_N0_dB = [-3:20]; % multiple Eb/N0 values
ipHat = zeros(1,N);
for ii = 1:length(Es_N0_dB)
ip = randsrc(1,N,alpha4pam);
s = (1/sqrt(5))*ip; % normalization of energy to 1
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

% demodulation
r = real(y); % taking only the real part

ipHat(find(r< -2/sqrt(5))) = -3;
ipHat(find(r>= 2/sqrt(5))) = 3;
ipHat(find(r>=-2/sqrt(5) & r<0)) = -1;
ipHat(find(r>=0 & r<2/sqrt(5))) = 1;

nErr(ii) = size(find([ip- ipHat]),2); % couting the number of errors
end

simBer = nErr/N;
theoryBer = 0.75*erfc(sqrt(0.2*(10.^(Es_N0_dB/10))));
close all
figure
semilogy(Es_N0_dB,theoryBer,'b.-');
hold on
semilogy(Es_N0_dB,simBer,'mx-');
axis([-3 20 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 4-PAM modulation')

