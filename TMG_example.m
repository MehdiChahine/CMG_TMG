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
R0 = [1 1;1 -1;];
c0 = -ones(N,1);

R1 = [0 1];
c1 = -1;
d1 = 1;


%% Precision matrix and mean vector
Q = eye(N);     % Precision matrix
mu = zeros(N,1); % Mean vector

%% Call CMG_TMG_STC
X = CMG_TMG(Q,mu,R0,c0,R1,c1,d1,Te,Tbin);

%% Display Scatter Plot
figure;hold on;grid on;box on;
plot(X(1,:),X(2,:),'b.')
xlabel('x_1')
ylabel('x_2')
xlim([-1.05 1.05])
ylim([-1.05 1.05])

%%
%Constraints display
x = -0.1:.01:10;
y = -1.:.01:1.;

l1 = -y + c0(1);
l2 = y - c0(2);

l3 = -1*ones(size(x));
l4 = 1*ones(size(x));
plot(y,l1,'k')
plot(y,l2,'k')
plot(x,l3,'k')
plot(x,l4,'k')
ylim([-1.2,1.2])
xlim([-1.2,1.5])

