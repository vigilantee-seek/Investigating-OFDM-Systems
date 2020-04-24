function [result] = ssl(Y,Xp,pilot_loc,Nfft,Nps,int_opt,sslm)
% Semi-supervised learning based channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% N = FFT size
% Nps = Pilot spacing
% int_opt = "LS" or "MMSE"
% sslm = "self-learning" (the only supported method in the current version)
% output:
% result = SSL Channel estimate
if int_opt == "LS"
    Np=Nfft/Nps; k=1:Np;
    LS_est(k) = Y(pilot_loc(k))./Xp(k);% LS channel estimation
    if sslm == "self-learning"
        
    end
elseif int_opt == "MMSE"
    
end
