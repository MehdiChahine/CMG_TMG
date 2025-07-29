clear all;close all;clc;
addpath('./lib/')
addpath('./Utility/')
%% Introduction
% This script provides an example of how to use the proposed sampler in the
% case of BOX Type Constaints STC (only)
% The choosen example is the unit ell_1 ball, i.e., ||X||_1 <= 1

%% Number of sampled
Tbin = 5000; % Number of Burn-in samples
Te = 5000;   % Number of desired samples

%% Restricted Domaine
% Unit ell_1 ball, i.e., ||X||_1 <= 1
N = 2;   % Dimension of the target distribution

% Linear ineqaulity matrix and bounds such that c <= Rx <= d
R = [1 1;1 -1;];
c = -ones(N,1);
d = ones(N,1);


%% Precision matrix and mean vector
Q = eye(N);     % Precision matrix
mu = ones(N,1); % Mean vector

%% Call CMG_TMG_STC
X = CMG_TMG_BTC(Q,mu,R,c,d,Te,Tbin);

%% Display Scatter Plot
figure;hold on;grid on;box on;
plot(X(1,:),X(2,:),'b.')
xlabel('x_1')
ylabel('x_2')
xlim([-1.05 1.05])
ylim([-1.05 1.05])