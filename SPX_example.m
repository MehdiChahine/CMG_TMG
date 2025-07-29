clear all;close all;clc;
addpath('./lib/')
addpath('./Utility/')
%% Introduction
% This script provides an example of how to use the proposed sampler in the
% case of Single Type Constaints STC (only)
% The choosen example is the unit simplex i.e., X_n >=0 and sum_n X_n <= 1

%% Number of sampled
Tbin = 5000; % Number of Burn-in samples
Te = 5000;   % Number of desired samples

%% Restricted Domaine
% Unit simplex, i.e., X_n >=0 and sum_n X_n <= 1
N = 2;   % Dimension of the target distribution

% Linear ineqaulity matrix and bounds such that c <= Rx
R = [eye(N);-ones(1,N)];
c = [zeros(N,1);-1];


%% Precision matrix and mean vector
Q = eye(N);     % Precision matrix
mu = ones(N,1); % Mean vector

%% Call CMG_TMG_STC
X = CMG_TMG_STC(Q,mu,R,c,Te,Tbin);

%% Display Scatter Plot
figure;hold on;grid on;box on;
plot(X(1,:),X(2,:),'b.')
xlabel('x_1')
ylabel('x_2')
xlim([-0.05 1.05])
ylim([-0.05 1.05])