function [H_LS] = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt)
% LS channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% N = FFT size
% Nps = Pilot spacing
% int_opt = Interpolation Methods (the same as the interp function)
% output:
% H_LS = LS Channel estimate
Np= Nfft/Nps; k = 1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k);% LS channel estimation
% Interpolation
H_LS = interpolate(LS_est,pilot_loc,Nfft,int_opt);
%  H_LS=zeros(2048,2048);
%  for i=1:2048
%      H_LS(i,i)=H_i(i);
%  end
% % H_LS=LS_est;