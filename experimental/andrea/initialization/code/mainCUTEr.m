close all
clear all
clc
addpath('C:\Users\skc\Documents\MATLAB\casadi-matlabR2014b-v3.0.0')
addpath('C:\Users\skc\Desktop\hpipm\hpipm\experimental\andrea\initialization\code\plotregion')
import casadi.*

% B_STRATEGY = 'MONOTONE'; % barrier strategy, possible values: {'MONOTONE', 'MEHROTRA'}
B_STRATEGY = 'MEHROTRA';
MAX_ITER = 100;
TAU0 = 1e6;
MAX_LS_IT = 100;
TH = 1e-16;
X0 = [1; 1];
PRINT_LEVEL = 1; % barrier strategy, possible values: {1, 2}
DTB = 0.1; % distance to boundaries

tol = 1e-1;
term_tol = 1e-6;
kappa = 0.3;

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
    
    alpha = [];
    tau = TAU0;
    
    iter.x  = zeros(nx,MAX_ITER);
    iter.l  = zeros(ne,MAX_ITER);
    iter.mu = zeros(nc,MAX_ITER);
    iter.s  = zeros(nc,MAX_ITER);
    
    iter.l(:,1) = 10;
    iter.x(:,1) = 0.5 * (lbx{num_prob,1} + ubx{num_prob,1});
    iter.mu = 0.1*ones(nc,MAX_ITER);
    
    iter.l(:,1) = 0;
    
    if sum(-C*iter.x(:,1) + d > DTB) == nc
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
        
        if PRINT_LEVEL == 2
            fprintf('it = %3.f   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', i, err_s,  err_e,  err_i,  err_c, tau, alpha);
        end
        
        err = norm(rhs);
        
        if strcmp(B_STRATEGY, 'MONOTONE')
            if err_s < tol && err_e < tol && err_i < tol && err_c < tol
                tau = kappa*tau;
            end
        end
        
        if err_s < term_tol && err_e < term_tol && err_i < term_tol && err_c < term_tol && tau < term_tol
            fprintf(' num_QP = %3.f   it = %3.f   solved = 1   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
            num_pass = num_pass + 1;
            break
        end
        
        J = full(j_rf(x, l, mu, s));
        
        if strcmp(B_STRATEGY, 'MEHROTRA')
            % compute duality measure
            dm = s.'*mu/nc;
            % solve affine system
            aff_rhs = full(rf(x, l, mu, s, 0));
            aff_step = -J\aff_rhs;
            % extract mu and s step
            aff_mu_step = aff_step(nx+ne+1:nx+ne+nc);
            aff_s_step  = aff_step(nx+ne+nc+1:end);
            % do line-search for aff
            aff_alpha = 1;
            for ls_iter = 1:MAX_LS_IT
                t_s  = s  + aff_alpha*aff_s_step;
                t_mu = mu + aff_alpha*aff_mu_step;
                if isempty(t_s(t_s < TH)) && isempty(t_mu(t_mu < TH))
                    break
                end
                if ls_iter == MAX_LS_IT
                    fprintf(' num_QP = %3.f   it = %3.f   solved = 0   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
                    break;
                    i = MAX_ITER;
                end
                aff_alpha = 0.5*aff_alpha;
            end
            % compute affine duality measure
            aff_dm = t_s.'*t_mu/nc;
            
            % compute centering parameter
            sigma = (aff_dm/dm)^3;
            
            % update tau
            tau = sigma*dm;
            % update complementarity residuals (build aggregated rhs)
            rhs = full(rf(x, l, mu, s, tau));
            rhs((nx+ne+nc+1:end)) = rhs(nx+ne+nc+1:end) + diag(aff_s_step)*aff_mu_step;
        end
        
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
                fprintf(' num_QP = %5.f   it = %5.f   solved = 0   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
                break;
                i = MAX_ITER;
            end
            alpha = 0.5*alpha;
        end
        
        iter.x(:,i+1)  = x  + alpha*dx;
        iter.l(:,i+1)  = l  + alpha*dl;
        iter.mu(:,i+1) = mu + alpha*dmu;
        iter.s(:,i+1)  = s  + alpha*ds;
        
        if i == MAX_ITER
            format shortE
            fprintf(' num_QP = %5.f   it = %5.f   solved = 0   err_s = %5.e    err_e = %5.e    err_i = %5.e    err_c =    %5.e    tau = %5.e    alpha = %5.e\n', num_prob, i, err_s,  err_e,  err_i,  err_c, tau, alpha);
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