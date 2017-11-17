close all
clear all
clc
addpath('C:\Users\skc\Documents\MATLAB\casadi-matlabR2014b-v3.0.0')
addpath('C:\Users\skc\Desktop\hpipm\hpipm\experimental\andrea\initialization\code\plotregion')
import casadi.*
load benchmark

nQP = size(H,1);
C_full = A;
n_ieq = nc;
n_eq = ne;
num_pass = 0;
grad = g;
clear g A nc ne
for num_prob = 1: nQP

    nc = 2*(n_ieq{num_prob,1}+nv{num_prob,1});
    nx = nv{num_prob,1};
    ne = n_eq{num_prob,1};

    A = [];
    b = [];
    lbC = [];
    ubC = [];
    C_ieq = [];
    for i = 1: size(C_full{num_prob},1)
        if lbA{num_prob}(i) == ubA{num_prob}(i)
           A = [A;C_full{num_prob}(i,:)];
           b = [b;lbA{num_prob}(i)];
        else
           C_ieq = [C_ieq;C_full{num_prob}(i,:)];
           lbC = [lbC;lbA{num_prob}(i)];
           ubC = [ubC;ubA{num_prob}(i)];
        end
    end
    
    C  = [C_ieq;-C_ieq;eye(nx);-eye(nx)];
        
    if size(A) ~= [ne,nx]
        display('Dimension of Ax are wrong');
    end
    
    if size(C,1) ~= [nc,nx]
        display('Dimension of Cx are wrong');
    end

    shift = 0;

    d = [ubC;-lbC;ubx{num_prob,1};-lbx{num_prob,1}] + [ shift*ones(n_ieq{num_prob,1},1); -shift*ones(n_ieq{num_prob,1},1);shift*ones(nx,1); -shift*ones(nx,1) ];

    Q = H{num_prob,1};
    q = grad{num_prob,1};

    x   = MX.sym('x', nx, 1);
    l   = MX.sym('l', ne, 1);
    mu  = MX.sym('mu', nc, 1);
    s   = MX.sym('s', nc, 1);
    tau = MX.sym('tau', 1, 1);

    f = 1/2*x.'*Q*x + q.'*x;
    if isempty(A)
        g = [];
    else
        g = A*x - b;
    end
    h = C*x - d;
    if ne > 0
        Lag = f + l.'*g + mu.'*h;
    else
        Lag = f + mu.'*h;
    end
    rf_exp = [  jacobian(Lag,x).'; ...
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
    iter.x(:,1) = 0.5 * (lbx{num_prob,1} + ubx{num_prob,1});
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

%        format shortE
%        fprintf('it = %5.e   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', i, err_s,  err_e,  err_i,  err_c, tau, alpha);

        err = norm(rhs);
        if err < tol
            tau = kappa*tau;
        end

        if err < term_tol && tau < term_tol
        format shortE
        fprintf(' num_QP = %5.e   it = %5.e   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
        num_pass = num_pass + 1;
%             fprintf('\n')
%             display(['-> solution found in ', num2str(i), ' iterations.']);
%             fprintf('\n')
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
            format shortE
            fprintf(' num_QP = %5.e   it = %5.e   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
          %  error('-> maximum number of iterations reached!')
        end
end


end
fprintf('Number of QP = %d,  Number of solved QP  = %d, ratio = %5.e', num_prob, num_pass, num_pass/num_prob);

% iter.x = iter.x(:,1:i);
% 
% figure()
% plot(iter.x(1,:), iter.x(2,:), '-')
% hold all
% plot(iter.x(1,end), iter.x(2,end), 'red*')
% plotregion(-C, -d, [], [], 'b')
% x1 = linspace(-10,10, 100);
% x2 = 1/A(2)*(b - A(1)*x1);
% plot(x1, x2, 'r')
% syms x1 x2
% f = 1/2*[x1; x2].'*Q*[x1; x2] + q.'*[x1; x2];
% ezcontour(f, [-1 1], [-1 1])
% title('')
% grid on