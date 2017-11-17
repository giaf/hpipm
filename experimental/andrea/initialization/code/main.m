close all
clear all
clc

import casadi.*

B_STRATEGY = 'MONOTONE'; % barrier strategy, possible values: {'MONOTONE', 'MEHROTRA'}
% B_STRATEGY = 'MEHROTRA';
MAX_ITER = 100;
TAU0 = 1e-1;
MAX_LS_IT = 100;
TH = 1e-16;

tol = 1e-1;
term_tol = 1e-6;
kappa = 0.3;

A = [-0.5, 1];
b = -0.1;

nc = 4;
nx = 2;
ne = 1;

C = [ 1 -1; ...
      1  1; ...
     -1  1; ...
     -1 -1];
 
shift = 1;

d = [1 0.5 1 0.5].';

Q = [10 0; 0 1];
q = [0.0; -0.];

x   = MX.sym('x', nx, 1);
l   = MX.sym('l', ne, 1);
mu  = MX.sym('mu', nc, 1);
s   = MX.sym('s', nc, 1);
tau = MX.sym('tau', 1, 1);

f = 1/2*[x - shift].'*Q*[x - shift] + q.'*[x - shift];
g = A*x - b;
h = C*x - d;

rf_exp = [  jacobian([f + l.'*g + mu.'*h],[x]).'; ...
            g; ...
            h + s; ...
            diag(s)*mu-ones(nc,1)*tau];
        
rf = Function('rf', {x, l, mu, s, tau}, {rf_exp});        
j_rf = Function('j_rf', {x, l, mu, s}, {jacobian(rf_exp, [x; l; mu; s])});  

% Newton method

clear x l mu s tau
max_iter = 100;

alpha = [];
tau = TAU0;

iter.x = zeros(nx,MAX_ITER);
iter.l = zeros(ne,MAX_ITER);
iter.mu = zeros(nc,MAX_ITER);
iter.s = zeros(nc,MAX_ITER);

iter.l(:,1) = 0;
iter.x(:,1) = [0; 0];

if sum(-C*iter.x(:,1) + d > 0) == nc
    iter.s(:,1) = -C*iter.x(:,1) + d;
    iter.mu(:,1) = tau./iter.s(:,1);
else
    iter.s(:,1) = sqrt(tau)*ones(nc,1);
    iter.mu(:,1) = sqrt(tau)*ones(nc,1);
end

for i = 1:MAX_ITER
    
    % compute step-szie
    x = iter.x(:,i);
    l = iter.l(:,i);
    mu = iter.mu(:,i);
    s = iter.s(:,i);
    
    % compute search direction
    rhs = full(rf(x, l, mu, s, tau));
    
    err_s = norm(rhs(1:nx));
    err_e = norm(rhs(nx+1:nx+ne));
    err_i = norm(rhs(nx+ne+1:nx+ne+nc));
    err_c = norm(rhs(nx+ne+nc+1:nx+ne+nc+nc));
    
    format shortE
    fprintf("it = %5.e   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n", i, err_s,  err_e,  err_i,  err_c, tau, alpha);
    
    err = norm(rhs);
    
    if strcmp(B_STRATEGY, 'MONOTONE')
        if err_s < tol && err_e < tol && err_i < tol && err_c < tol
            tau = kappa*tau;
        end
    end
    
    if err_s < term_tol && err_e < term_tol && err_i < term_tol && err_c < term_tol && tau < term_tol
        fprintf('\n')
        display(['-> solution found in ', num2str(i), ' iterations.']);
        fprintf('\n')
        break
    end
    
    J = full(j_rf(x, l, mu, s));
    step = -J\rhs;
    
    % extract search direction componenents
    dx = step(1:nx);
    dl = step(nx+1:nx+ne);
    dmu = step(nx+ne+1:nx+ne+nc);
    ds = step(nx+ne+nc+1:end);
    
    % do line-search
    alpha = 1;
    for ls_iter = 1:MAX_LS_IT
        t_s  = s  + alpha*ds;
        t_mu = mu + alpha*dmu;
        if isempty(t_s(t_s < TH)) && isempty(t_mu(t_mu < TH))
            break
        end
        if ls_iter == MAX_LS_IT
            display('-> line-search failed! Exiting.')
            i = MAX_ITER;
        end
        alpha = 0.5*alpha;
    end
    
    iter.x(:,i+1)  = x  + alpha*dx;
    iter.l(:,i+1)  = l  + alpha*dl;
    iter.mu(:,i+1) = mu + alpha*dmu;
    iter.s(:,i+1)  = s  + alpha*ds;
    
    if i == MAX_ITER
        display('-> maximum number of iterations reached! Exiting.')
        i = MAX_ITER;
    end
end

iter.x = iter.x(:,1:i);

figure() 
hold all
plotregion(-C, -d, [], [], 'b')
x1 = linspace(-10,10, 100);
x2 = 1/A(2)*(b - A(1)*x1);
plot(x1, x2, 'r')
syms x1 x2
f = 1/2*[ [x1; x2] - shift ].'*Q*[ [x1; x2] - shift ] + q.'*[ [x1; x2] - shift ];
ezcontour(f, [-1, 1], [-1, 1])
title('')
grid on

plot(iter.x(1,:), iter.x(2,:), '-')
plot(iter.x(1,end), iter.x(2,end), 'red*')
