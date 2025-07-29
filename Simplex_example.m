clear all;close all;clc;
addpath('./lib/')
% Required samples & Burn-in samples
Tbin = 5000;
Te = 5000;

% Sampling Truncated Gaussian restrincted to the unit simplex in R^N
N = 2;
% Upper and lower bounds : R c d
R = [eye(N);-ones(1,N)];
c = [zeros(N,1);-1];


%
Sigma = eye(N);
Q = inv(Sigma);
mu = ones(N,1);


X = CMG_TMG_STC(Q,mu,R,c,Te,Tbin);


figure;hold on;grid on;box on;
plot(X(1,:),X(2,:),'b.')
xlabel('x_1')
xlabel('x_2')