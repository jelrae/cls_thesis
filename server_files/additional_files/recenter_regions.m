%% Set up paths
addpath('/home/jordan/neuro_thesis/neuro_thesis/cls_thesis/gc_hierarchies/')
addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
%addpath('D:/MATLAB/mvgc_v1.0')
startup
addpath('/home/jordan/common/matlab/fieldtrip-20210411')
ft_defaults

format short;
clear all;
close all;clc;
disp('starting\n')

monkey = 'pele';

if strcmp(monkey, 'kurt')
    load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_gc_one_v_all.mat')
else
    load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_gc_one_v_all.mat')
end

%% Split the data into 3 buckets for theta beta and gamma

%theta = 4-7; p tol = 4
%beta = 16-31 ; p tol = 5
%gamma = 32 -100 ; p tol = 10

df = fnq/(fres+1);

theta_bottom = floor(4/df);
theta_top = ceil(7/df);

beta_bottom = floor(21/df);
beta_top = ceil(26/df);

gamma_bottom = floor(42/df);
gamma_top = ceil(90/df);

theta_x_range = 1:1:8;
beta_x_range = 1:1:10;
gamma_x_range = 1:1:20;

x_range = (1:1:fres+1)./(fres+1);
x_range = x_range.*fnq;

test_forward = gc_forward{1};
test_backward = gc_backward{1};

test_theta = test_forward(:,theta_bottom:theta_top);
test_theta_range = x_range(theta_bottom:theta_top);
test_theta_matrix = repmat(test_theta_range, size(test_theta,1),1);

[M_t,I_t] = max(test_theta,[],2);

p = polyfit(test_theta_matrix', test_theta', 2);

fitted = polyval(p, test_theta_range);

% dtheta = []
% 
% for i = 2:length(test_theta_range)
%     dtheta(end+1,:) = 
% end

plot(test_theta_range, fitted);

% test_beta = test_forward(:,beta_bottom:beta_top);
% test_beta_range = x_range(beta_bottom:beta_top);
% test_beta_matrix = repmat(test_beta_range, size(test_beta,1),1);
% 
% [M_b,I_b] = max(test_beta,[],2);
% 
% p = polyfit(test_beta, test_beta_matrix, 2);