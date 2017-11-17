close all
clear all
clc

import casadi.*

A = [-0.1, 1];
b = -0.1;

nc = 4;
nx = 2;
ne = 1;

C = [ 1 -1; ...
      1  1; ...
     -1  1; ...
     -1 -1];
 
shift = 0;

d = [1 0.05 1 0.05].' + [ shift; shift; -shift; -shift ];

Q = diag([0.001, 100]);
q = [0.5; -0.2];

x   = MX.sym('x', nx, 1);
l   = MX.sym('l', ne, 1);
mu  = MX.sym('mu', nc, 1);
s   = MX.sym('s', nc, 1);
tau = MX.sym('tau', 1, 1);

f = 1/2*x.'*Q*x + q.'*x;
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

tau = 1e1;
max_ls_it = 100;
th = 1e-8;

tol = 1e0;
term_tol = 1e-6;
kappa = 0.3;
alpha = [];

iter.x = zeros(nx,max_iter);
iter.l = zeros(ne,max_iter);
iter.mu = sqrt(tau)*ones(nc,max_iter);
iter.s = sqrt(tau)*ones(nc,max_iter);

iter.l(:,1) = 10;
iter.x(:,1) = [1;-100];
iter.mu = 0.1*ones(nc,max_iter);

for i = 1:max_iter
    
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
    if err < tol
        tau = kappa*tau;
    end
    
    if err < term_tol && tau < term_tol
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
    for ls_iter = 1:max_ls_it
        t_s  = s  + alpha*ds;
        t_mu = mu + alpha*dmu;
        if isempty(t_s(t_s < th)) && isempty(t_mu(t_mu < th))
            break
        end
        if ls_iter == max_ls_it
            error('-> line-search failed!')
        end
        alpha = 0.5*alpha;
    end
    
    iter.x(:,i+1)  = x  + alpha*dx;
    iter.l(:,i+1)  = l  + alpha*dl;
    iter.mu(:,i+1) = mu + alpha*dmu;
    iter.s(:,i+1)  = s  + alpha*ds;
    
    if i == max_iter
        error('-> maximum number of iterations reached!')
    end
end

iter.x = iter.x(:,1:i);

figure() 
plot(iter.x(1,:), iter.x(2,:), '-')
hold all
plot(iter.x(1,end), iter.x(2,end), 'red*')
plotregion(-C, -d, [], [], 'b')
x1 = linspace(-10,10, 100);
x2 = 1/A(2)*(b - A(1)*x1);
plot(x1, x2, 'r')
syms x1 x2
f = 1/2*[x1; x2].'*Q*[x1; x2] + q.'*[x1; x2];
ezcontour(f, [-1 1], [-1 1])
title('')
grid on