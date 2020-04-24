%% Main Program Information
% Main Scope: Channel estimation
% Modulation: OFDM systems
% Channel: Rayleigh multipath fading channel
% Algorithms: Least Square (LS) and Minimum Mean Square Errors (MMSE)
% Environment: MATLAB R2019b

%% Parameter Setting and Initialization
clc;clear;close all;

Nfft = 2048; % number of FFT points
Ng = 512;
Nofdm = 2560; 
Nsym = 100;
Nps = 4; %Pilot Spacing
Np = Nfft/Nps; %Number of pilots per OFDM symbol
Nbps = 4;
itr = 10; % Number of Iterations

% Reception at the receiver end for various SNR
snr_list = 5:2:25;  % SNR range
nn = length(snr_list);  %Number of SNR observations

%Initialising MSE arrays for different SNR observations:
z_linin = zeros(1,nn);
z_splin = zeros(1,nn);
z_lindft = zeros(1,nn);
z_spldft = zeros(1,nn);
zmmse = zeros(1,nn);
zmmse_dft = zeros(1,nn);

%% Channel code
dopjakes = doppler('Jakes');

nsd1 = [0.05, 0.1];
dopgauss1 = doppler('BiGaussian', ...
    'NormalizedStandardDeviations', nsd1,...
    'NormalizedCenterFrequencies', [-0.8, 0.4],...
    'PowerGains', [1,0.1].*sqrt(2*pi*nsd1.^2));

nsd2 = [0.1, 0.15];
dopgauss2 = doppler('BiGaussian',...
    'NormalizedStandardDeviations', nsd2,...
    'NormalizedCenterFrequencies', [0.7, -0.4],...
    'PowerGains', [1,0.1^1.5].*sqrt(2*pi*nsd2.^2));

fd=130; %Maximum Doppler Shift, hertz
ts=(7/64)*10^-6; %Sampling Time, second

%Rayleigh Channel - Multifading Channel
chan = comm.RayleighChannel('SampleRate', 1/ts,...
    'MaximumDopplerShift', fd,...
    'PathDelays', [0.0 0.2 0.5 1.6 2.3 5.0] * 1e-6,...
    'AveragePathGains', [-3 0 -2 -6 -8 -10],...
    'DopplerSpectrum', {dopjakes,dopjakes,dopjakes,...
    dopgauss1,dopgauss2,dopgauss2},...
    'NormalizePathGains', 0);

%% Main Program
for t=1:itr
    %% Error Cumulation Vectors
    z = zeros(1,nn); 
    z1 = zeros(1,nn); 
    z2 = zeros(1,nn);
    z3 = zeros(1,nn);
    zm = zeros(1,nn);
    zms = zeros(1,nn);

    %% Baseband Signal
    % 4 - QAM - Modulation Scheme 
    M = 4;
    % hmod = qammod('M',M, 'symOrder','gray');
    Es = 1; 
    A = sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
    for nsym=1:Nsym
        Xp = 2 * (randn(1,Np)>0)-1; % Pilot sequence generation
        msgint = randi(1, Nfft-Np, M); % bit generation
        dat_ser = A * qammod(msgint, M, 'gray');
    end

    % serial to parllel conversion
    dat_par = dat_ser.';

    % Pilot Insertion - Comb Type Arrangement
    counter = 0; 
    loc = [];
    X = zeros(1,Nfft);
    for i = 1:Nfft
        if mod(i,Nps) == 1
            X(i) = Xp(floor(i/Nps)+1); 
            loc = [loc i]; 
            counter = counter+1;
        else
            X(i) = dat_par(i-counter);
        end
    end

    % inverse discret Fourier transform (IFFT)
    X_ifft = ifft(X,Nfft);

    % Adding Cyclic Prefix - Guard interval
    guard = X_ifft(:,end-511:end); % this is the Cyclic Prefix part to be appended.
    ofdm_par = [guard X_ifft];

    % parallel to serial - Generation of the OFDM Signal 
    ofdm = ofdm_par.';

    %% Signal Transmission
    %Passing OFDM Signal through the created Channel 
    chan_op = chan(ofdm);
    XFG = fft(chan_op);

    %Channel Covariance Matrix for MMSE estimation  
    x_test = randn(6,1);
    X_test = fft(x_test,Nfft);
    y_test = chan(x_test);
    Y_test = fft(y_test,Nfft);
    H1 = Y_test' / X_test';
    h=ifft(H1,6);
    H=fft(h,Nfft);
    ch_length=length(h);

    %% Received Signal Processing
    %Reception at the receiver end for various SNR
    for snr_cnt = 1:nn
        %% Parameters
        SNR = snr_list(snr_cnt);
        
        %% Noise Setting (AWGN Modelling)
        n1 = ones(2560,1);
        %Just to ensure that the function awgn adds 'complex gaussian noise'..
        n1 = n1 * 0.000000000000000001i;
        noise = awgn(n1,SNR);
        variance = var(noise);
        N = fft(noise);
        
        %% Received Signal
        Y_rec=XFG+N;
        y_ser=ifft(Y_rec);
        y_par=y_ser.';              % serial to parallel conversion
        y = y_par(Ng+1:Nofdm);      % guard interval
        Y = fft(y);                 % FFT

        %% channel estimation

        %LS Estimator with Linear Interpolator 
        H_est = LS_CE(Y,Xp,loc,Nfft,Nps,'linear');
        err = (H-H_est)*(H-H_est)';
        z(snr_cnt) = err/(Nfft*Nsym);

        %LS Estimator with Linear Interpolator - DFT Based
        h_est = ifft(H_est);
        h_DFT = h_est(1:ch_length);
        H_DFT = fft(h_DFT,Nfft); 
        err = (H-H_DFT)*(H-H_DFT)';
        z3(snr_cnt) = err/(Nfft*Nsym);

        %LS Estimator with Spline Cubic Interpolator
        H_est = LS_CE(Y,Xp,loc,Nfft,Nps,'spline');
        err=(H-H_est)*(H-H_est)';
        z1(snr_cnt) = err/(Nfft*Nsym);

        %LS Estimator with Spline Cubic Interpolator - DFT Based
        h_est = ifft(H_est);
        h_DFT = h_est(1:ch_length);
        H_DFT = fft(h_DFT,Nfft); 
        err=(H-H_DFT)*(H-H_DFT)';
        z2(snr_cnt) = err/(Nfft*Nsym);  

        % MMSE Estimator
        H_est = MMSE_CE(Y,Xp,loc,Nfft,Nps,h,SNR);
        err=(H-H_est)*(H-H_est)';
        zm(snr_cnt) = err/(Nfft*Nsym);

        % MMSE Estimator - DFT Based
        h_est = ifft(H_est);
        h_DFT = h_est(1:ch_length);
        H_DFT = fft(h_DFT,Nfft); 
        err=(H-H_DFT)*(H-H_DFT)';
        zms(snr_cnt) = err/(Nfft*Nsym);
    end

    z_linin = z_linin + z;
    z_splin = z_splin + z1;
    z_lindft = z_lindft + z2;
    z_spldft = z_spldft + z2;
    zmmse = zmmse + zm;
    zmmse_dft = zmmse_dft + zms;
end

%% Draw Figures
figure(1)
semilogy(snr_list,(1/itr)*z_linin,'r+:', snr_list,(1/itr)*z_splin,'bo:', snr_list,(1/itr)*z_lindft,'--xg', snr_list,(1/itr)*z_spldft,'--sc');
legend('LS - Linear Interpolation','LS - Spline Cubic Interpolation','LS - Linear Interpolation(DFT)','LS - Spline Cubic Interpolation(DFT)');
xlabel('SNR');
ylabel('MSE');
grid on
hold on

figure(2)
semilogy(snr_list,(1/itr)*zmmse,'r+:',snr_list,(1/itr)*zmmse_dft,'bo:');
xlabel('SNR');
ylabel('MSE');
legend('MMSE','MMSE - DFT Based');
hold on
grid on
 

figure(3)
semilogy(snr_list,(1/itr)*z_linin,'r+:', snr_list,(1/itr)*z_splin,'bo:', snr_list, (1/itr)*zmmse,'--xg');
xlabel('SNR');
ylabel('MSE');
legend('LS - Linear Interpolation','LS - Spline Cubic Interpolation','MMSE');
hold on
grid on 