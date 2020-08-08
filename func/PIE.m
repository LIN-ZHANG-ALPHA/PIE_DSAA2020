

function [completed_x,beta_output,detect_periods_vector] = PIE(x,opts)


% for period dictionary
if ~isfield(opts, 'Dictionary_type'),   opts.Dictionary_type   = 'Ramanujan'; end  %
if ~isfield(opts, 'Pmax_LP')&&~isfield(opts, 'Pmax'),         opts.Pmax_LP   = 50; end  %

if ~isfield(opts, 'A'),            opts.A = Create_Dictionary(opts.Pmax_LP,size(x,1),opts.Dictionary_type); end  %
if ~isfield(opts, 'Penalty_type'), opts.Penalty_type   = 'square'; end 


if ~isfield(opts, 'halting_Thr'), opts.halting_Thr = 1e-3; end  % for algo  
if ~isfield(opts, 'max_iter'),    opts.max_iter = 100; end  % for algo 

if ~isfield(opts, 'visual'),    opts.visual = 0; end  % for algo 

if ~isfield(opts, 'rho'),      opts.rho =1e-5; end 

if isfield(opts, 'Pmax')&&~isfield(opts, 'Pmax_LP'), opts.Pmax_LP = opts.Pmax; end

Pmax        = opts.Pmax_LP;
halting_Thr = opts.halting_Thr; % for outer loop

% Penalty Vector Calculation
penalty_vector = [];
for i=1:Pmax
    k=1:i;k_red=k(gcd(k,i)==1);k_red=length(k_red);
    penalty_vector=cat(1,penalty_vector,i*ones(k_red,1));
end

penalty_vector = penalty_vector.^2;
H_inv          =  diag(1./penalty_vector); % inverse of diagnal matrix


Adj = get_adj(penalty_vector);
D   = diag(sum(Adj, 2));
L   = D-Adj; % Laplacian matrix


QUIET    = 0;
t_start = tic;

%% pre-compute
A       =  opts.A; % periods dictionary
A_hat   =  A* H_inv;  % H is diagnal penlty matrix
xlen    =  length(x);
one_vec =  ones(xlen,1);

mask    =  double(logical(x)); % mask of x: 1/0, binary vector
tmp1    = one_vec + 2*opts.lambda_0*mask;

lwx = opts.lambda_0 * mask.*x;

tmp2     =  A_hat'*A_hat + 2*opts.lambda_2*L + opts.rho*eye(size(L));
inv_tmp2 =  tmp2^(-1);

%% initalization
z_t   = x; 
obj_t = 10^10;

%% main loop from here

for iter =  1: opts.max_iter
    
    %% update beta : sparse group lasso
    
    [alpha,beta_t_1] = update_beta(A_hat,z_t,L,inv_tmp2,opts);
    
    
    %% update z
    
    z_t_1 = (2*lwx + A_hat*beta_t_1)./tmp1;
    
    
    %% check covergence
    
    obj_t_1 = get_object(A_hat,x,mask,z_t_1,beta_t_1,L,opts);
    obj_resi = abs(obj_t_1 - obj_t) /obj_t;
    
    if obj_resi < halting_Thr
        break
    end
    
    if mod(iter,100) ==0
        disp(['Iter:',num2str(iter),'::','obj_res = ',num2str(abs(obj_resi))]);
    end
    
    
    
    %% update variable
    obj_t   = obj_t_1;
    z_t     = z_t_1;
    beta_t  = beta_t_1;
end

completed_x = z_t;
beta_output = beta_t;

s =  H_inv * alpha;
% s =  beta_output;

if ~QUIET
    toc(t_start);
end



%% get periods

% AA = A*H_inv;

energy_s = 0.*[1:Pmax];

current_index_end = 0;
for i=1:Pmax
    k_orig = 1:i;k=k_orig(gcd(k_orig,i)==1);
    current_index_start = current_index_end + 1;
    current_index_end = current_index_end + size(k,2);
    
    for j=current_index_start:current_index_end
        
        energy_s(i) = energy_s(i)+((abs(s(j)))^2);
        
        periods_zone{i} = [current_index_start,current_index_end];
        
        %subdictionary_energy(i) = norm(AA(:,j)*T(j));
        
    end
    
end

energy_s(1) = 0;
% figure, stem(energy_s,'linewidth',3,'color',[0 0 0]);
% title('LP incomplete norm minimization');
% xlabel('Period');
% ylabel('Strength');

% energy_periods_output =  energy_s;

detect_periods_vector =  energy_s./norm(energy_s);
% figure, stem(detect_periods_vector,'linewidth',3,'color',[0 0 0]);
% title('LP incomplete norm minimization');
% xlabel('Period');
% ylabel('Strength');

end

function obj = get_object(A_hat,x,w,z,beta,L,opts)

obj = 0.5*(norm(z-A_hat*beta))^2+ opts.lambda_0*norm((x-z).*w)+ opts.lambda_1*norm(beta,1)+ opts.lambda_1* beta'*L*beta;

end


