function [output, S, fposDer, fnegDer]= kalman_spike(Th, input, S, fposDer, fnegDer)

%% Kalman recursive algorithm with an extension for spikes detection & correction

% input                input variable, i.e. raw signal at time point t
% output               filtered value
% Th                   despiking threshold
% S                    structure containing covariance matrices Q, R and state vector x
% spikeLength          length of the spike, considered unique
% fposDer              counting variables used to detect the subsequent positive or negative spikes, i.e. if spikeLength > 1 
% fnegDer 

% for Kalman filter start up see the reference

% Reference: 
% Koush, Y., M. Zvyagintsev, M. Dyck, K. A. Mathiak, K. Mathiak, 2012.
% Signal quality and Bayesian signal processing in neurofeedback based on real-time fMRI, NeuroImage 59, 478-489.

% Written by Yury Koush, 2011-2012, UKA Aachen, Germany, EPFL Lausanne, Switzerland
% Contact: yury.koush@epfl.ch, yurykoush@gmail.com

%% Kalman algorithm
A = 1; 
H = 1;
I = 1;
S.x = A*S.x;                
S.P = A*S.P*A' + S.Q;        
K = S.P*H'*pinv(H*S.P*H' + S.R);
tmp_x = S.x;
tmp_p = S.P;
delta = K*(input - H*S.x);
S.x = S.x + delta;
S.P = (I - K*H)*S.P;
spikeLength = 1; % considered a spike with an unique length

%% Spikes correction before the Kalman 'output' is assigned:
% - covariance matrice (P) and state variables (x) are taken from previous iteration (tmp_x & tmp_p) and not
%   updated if the spike is detected
% - sign of the spike is defined by sign of delta

% spike removal
if abs(delta) < Th 
    output = H*S.x;
    fnegDer = 0;
    fposDer = 0;
else
    if delta > 0
        if  fposDer < spikeLength
            output = H*tmp_x;
            S.x = tmp_x; 
            S.P = tmp_p;
            fposDer = fposDer + 1;
        else
            output = H*S.x; 
            fposDer = 0;
        end
    else
        if  fnegDer < spikeLength
            output = H*tmp_x;
            S.x = tmp_x; 
            S.P = tmp_p;
            fnegDer = fnegDer + 1;
        else
            output = H*S.x; 
            fnegDer = 0;
        end
    end
end
