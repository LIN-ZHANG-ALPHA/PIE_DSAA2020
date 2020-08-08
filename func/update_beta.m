

function [alpha,beta_new] = update_beta(A,z,L,inv_tmp,opts)

%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| A*beta - z ||_2^2 + \lambda || z ||_1 + \lambda beta' *L* \lambda
%   s.t \alpha = \beta

% for ADMM
if ~isfield(opts, 'eta'),                 opts.epsilon =1e-4; end  %
if ~isfield(opts, 'tau'),                 opts.tau     =1.1; end  %
if ~isfield(opts, 'rho'),                 opts.rho     =1e-6; end  %
if ~isfield(opts, 'max_rho'),             opts.max_rho =1e10; end  %
if ~isfield(opts, 'max_max_iter_beta'),   opts.max_max_iter_beta = 100; end  %

rho     = opts.rho;
epsilon = opts.epsilon ;
tau     = opts.tau;
max_rho = opts.max_rho ;

% max_iter = opts.max_max_iter_beta;

QUIET    = 1;
MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
%% initialization
[m, n]  = size(A);

x       =  zeros(n,1);
y       =  zeros(n,1); % Langarian multiplier
u       =  zeros(n,1); % euqals to y_k/rho
alpha   =  zeros(n,1);  % s.t. alpha = beta

obj_k   = 10^10;

% pre-compute
Atz     = A'*z;
[LA UA] = factor(A, rho);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
        'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

%% main loop

for k = 1: MAX_ITER
    
    %++++++++++++++++++++ beta-update
     q = Atz + rho*(alpha - u);    % temporary value
%     if( m >= n )    % if skinny
%         beta_new = UA \ (LA \ q);
%     else            % if fat
%         beta_new = q/rho - (A'*(UA \ ( LA \ (A*q) )))/rho^2;
%     end
%     
    beta_new = inv_tmp*q;    
    
    %++++++++++++++++++++ alpha-update with relaxation
    % alpha_k_1  = soft_thresholding(beta_k_1 + y_k/rho, opts.lambda_1/rho);
    % alpha_k_1  = shrinkage(beta_k_1 + y_k/rho, opts.lambda_1/rho);
    alpha_old = alpha;
    alpha     = shrinkage(beta_new + u, opts.lambda_1/rho);
    
    % u-update: lagarian multiplier
    u      = u + (beta_new - alpha);
    
    % ++++++++++++++++++++ update penalty parameter
    % rho        =   min(tau*rho, max_rho);
    
    % diagnostics, reporting, termination checks
    history.objval(k)  = get_obj(z, A, L, beta_new, opts);
    
    history.r_norm(k)  = norm(beta_new-alpha);
    history.s_norm(k)  = norm(-rho*(alpha - alpha_old));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-alpha));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
    
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end
    if (history.r_norm(k) < history.eps_pri(k) && ...
            history.s_norm(k) < history.eps_dual(k))
        break;
    end
    
    
    % update variables
    beta_old  =  beta_new;
    %alpha_old =  alpha;
    % obj_k   =  obj_k_1;
end

beta_output = beta_old;
end

function obj = get_obj(z, A,L, beta, opts)
obj =   0.5*sum((A*beta - z).^2) + opts.lambda_1*norm(beta,1) + opts.lambda_2* beta'*L*beta;
end

function z = shrinkage(x, kappa)
z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L U] = factor(A, rho)
[m, n] = size(A);
if ( m >= n )    % if skinny
    L = chol( A'*A + rho*speye(n), 'lower' );
else            % if fat
    L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end