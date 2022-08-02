clc;
clear all;

%=============================================================================
%=============================================================================
%           Non Orthogonal Multiple Access (NOMA) Simulation
%       Version 1: Encode only: Transmitter side, no modulation
%=============================================================================
% This code was based on codes made by Simith (SMT) and Thanh Nguyen (TNY)
% published in:
% Website:  MathWorks (R)   -   www.mathworks.com
%   Forum:  MATLAB Answers [TM]   -   /matlabcentral/answers/
%   Title:  I am struggling with a code.the code for signal transmition in
%           Non Orthogonal Multiple Access. Please help me.
%   Asked:  simith on 15 Mar 2017   -   SMT
%  Answer:  Thanh Nguyen on 21 Apr 2017   -   TNY
%=============================================================================
%=============================================================================

%=============================================================================
%   Variable Names:
%     Variable names were changed according to a notation 'similar' to the
%     'Hungarian Notation' and structured programming for better 
%     understanding of NOMA techniques since many MATLAB Users are not use
%     to the meaning of equation variables related to NOMA mathematics.
%
%                +--------- new NOMA variable name
%                |
%            +---+---+
%             VarName_XXX
%                    +-+-+
%                      |
%                      +--- original SMT/TNY variable name,
%                           '_00' if the variable isn't in SMT/TNY codes
%=============================================================================

% n0.0f bits for transmit signal: 4 bits only - REMOVE '/25' for 100 bits
% as original SMT/TNY codes;
        TxBits_n = 100/25;

% distance from Base Station (BS) to Primary User, to User1, to User2 - SMT
% Assuming maximum distance as 10, for attenuation calculation purposes
        DstBStoPUser_dp = 3;
        DstBStoUser1_ds1 = 2;
        DstBStoUser2_ds2 = 4;
        MaxDsttoUser_00 = 10;

      SumSquareDst_00=DstBStoPUser_dp^2+DstBStoUser1_ds1^2+DstBStoUser2_ds2^2;

% signal from BS - Power allocation already applied on signal(???) - SMT:
%       xp=rand(1,n)>0.5;
%       xs1=rand(1,n)>0.2;
%       xs2=rand(1,n)>0.3;
%
% power allocation for Primary User, User1 and User 2 - TNY:
%       pp=0.5;
%       p1=0.3; inverse from SMT: p1<->p2
%       p2=0.2; inverse from SMT: p1<->p2
%
% NOMA corrected: farther away from Base Station (BS), more allocated power;
% Assuming the Base Station (BS) has total power of 1 and will allocate
% power for each User as proportional to squared distance: Pwr ~ Dst^2
%
        TotPwrBS_00 = 1.0;

          % previously suggested: 0.3
          PwrPUser_pp = TotPwrBS_00*(DstBStoPUser_dp^2)/SumSquareDst_00;
          % previously suggested: 0.2
          PwrUser1_p1 = TotPwrBS_00*(DstBStoUser1_ds1^2)/SumSquareDst_00;
          % previously suggested: 0.5
          PwrUser2_p2 = TotPwrBS_00*(DstBStoUser2_ds2^2)/SumSquareDst_00;

%=============================================================================
    %%%Create random binary messages/signals from Base Station (BS)
%=============================================================================
% signal of 'n' bits from BS to Primary User, User1 and User2 based on
% with power allocation already applied??? - SMT:
%       xp=rand(1,n)>0.5;
%       xs1=rand(1,n)>0.2;
%       xs2=rand(1,n)>0.3;
% signal stream of Primary User, User1 and User2 - TNY:
% 'rand' generates any number between [0 and 1], NOT binary, noise included???
%       xp=rand(1,n);
%       xs1=rand(1,n);
%       xs2=rand(1,n);
%
% Correct (actual) binary messages of 'TxBits_n' bits length:
% equal probability of 0 and 1 in every bit;
%   'rand'  generates numbers in [0 to 1], uniformally distributed;
%   Mean is 0.5
%
        SgnPUser_xp = rand(1,TxBits_n) > 0.5;
        SgnUser1_xs1 = rand(1,TxBits_n) > 0.5;
        SgnUser2_xs2 = rand(1,TxBits_n) > 0.5;

%=============================================================================
    %%%Superposition Encoding
%=============================================================================
% Direct sum of signals: incorrect  - SMT: X=xp+xs1+xs2;
% NOMA: Power-domain Multiplexing, sum of products signal*sqrt(power) - TNY:
%
        Enc_X = sqrt(PwrPUser_pp)*SgnPUser_xp;
        Enc_X = sqrt(PwrUser1_p1)*SgnUser1_xs1 + Enc_X;
        Enc_X = sqrt(PwrUser2_p2)*SgnUser2_xs2 + Enc_X;

%=============================================================================
    %%%Received signals for all Users
%=============================================================================
% Adding Gaussian Noise: use 'randn' instead of 'rand':
%   'randn' generates numbers in [-Inf,+Inf], normally distributed (Gaussian);
%   Mean is zero, but with strong concetration in [-1 to +1];
%   'rand'  generates numbers in [0 to 1], uniformally distributed;
%   Mean is 0.5
%
% NOMA: Additive White Gaussian Noise (AWGN) with ZERO MEAN and double-side
%  power spectral density, N0/2.
% Noise variation on time N(t) (addition): different for every bit;
%
% Since 'randn' concetrates in [-1 to +1] or even larger and a bit is only
% [0 or 1], and Power<=1 the Signal-to-Noise Ratio (SNR) could be too low.
% So a constant to reduce Noise level is necessary.
%
        NoiseReduc_0 = 10; 
        NoisePUser_N = randn(1,TxBits_n)/NoiseReduc_0;
        NoiseUser1_N = randn(1,TxBits_n)/NoiseReduc_0;
        NoiseUser2_N = randn(1,TxBits_n)/NoiseReduc_0;

% Channel Attenuation Gain (multiplier): different for every User/channel,
% no variation on time. Attenuation is inversely proportional to the power
% and directly proportional to squared distance.
% As the allocated power is proportional to squared distance: Pwr ~ Dst^2,
% it makes all Attenuations become a "boring" constant (0.71 in this case).
% So a random System Loss inversely proportional to power was included to
% increase the unpredictability of simulation.
%
  SyLosPUser_00 = 1 + rand/(10*PwrPUser_pp);
  AtnGnPUser_00 = TotPwrBS_00/PwrPUser_pp * 1/SyLosPUser_00;
  AtnGnPUser_00 = AtnGnPUser_00*(DstBStoPUser_dp^2)/(MaxDsttoUser_00^2);
  AtnGnPUser_00 = 1 - AtnGnPUser_00;

    SyLosUser1_00 = 1 + rand/(10*PwrUser1_p1);
    AtnGnUser1_00 = TotPwrBS_00/PwrUser1_p1 * 1/SyLosUser1_00;
    AtnGnUser1_00 = AtnGnUser1_00*(DstBStoUser1_ds1^2)/(MaxDsttoUser_00^2);
    AtnGnUser1_00 = 1 - AtnGnUser1_00;

    SyLosUser2_00 = 1 + rand/(10*PwrUser2_p2);
    AtnGnUser2_00 = TotPwrBS_00/PwrUser2_p2 * 1/SyLosUser2_00;
    AtnGnUser2_00 = AtnGnUser2_00*(DstBStoUser2_ds2^2)/(MaxDsttoUser_00^2);
    AtnGnUser2_00 = 1 - AtnGnUser2_00;

%=============================================================================
    %%%Signal received by each User
%=============================================================================
% SMT: bit-lenght iteraction adding a Noise, increasing per bit - incorrect;
% Adding signal/(squared distance) - incorrect;
%       for s=0.1:0.01:1    % noise variation
%           N=s*randn(1,n);
%           ys1=xs1/ds1.^2+N+xp/dp.^2+xs2/ds2.^2
%       end
%
% NOMA: the Superposition Encoding (Enc_X) containing the messages to all
% Users already calculated above is affected by Attenuation (multiplier) and
% Noise (additive) just one time only as in Linear equation below:
%
%         Yk(t) = X(t).Gk + Wk(t)     where
%
% Yk(t) = superimposed signal received by User[k];
% X(t)  = superimposed signal with all Users messages as transmitted by BS;
% Gk    = channel attenuation gain for the link between BS and User[k];
% Wk(t) = additive White Gaussian Noise (AWGN) at the User[k] with
%         mean ZERO and density N0;
% 
        RxSgnPUser_ysp = Enc_X*AtnGnPUser_00 + NoisePUser_N;
        RxSgnUser1_ys1 = Enc_X*AtnGnUser1_00 + NoiseUser1_N;
        RxSgnUser2_ys2 = Enc_X*AtnGnUser2_00 + NoiseUser2_N;