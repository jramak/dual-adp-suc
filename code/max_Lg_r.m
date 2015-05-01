function [best_cost,best_idx,opt_cost_vec,time,tot_time,L_g_cost,...
    z_sol_avg] = max_Lg_r(filename,maxIters,rho0,stepVec)

% load in from data file
irgdx(filename);

% need to assign data for the sake of the parfor loop
numPts = numPts; 
q = q;
c = c;
c_bar = c_bar;
h_bar = h_bar;
Lu = Lu;
Ld = Ld;
p_s = p_s;
Ru = Ru;
Rd = Rd;
D_s_ub = D_s_ub;
D_s_idx_ub = D_s_idx_ub;
q_min = q_min;
q_max = q_max;

numGens = length(Lu) - 2; % number of generators not including buy/sell
T = 168;
quadPts = length(p_s); % number of quadrature points in distribution
% power level step size for each generator
Pstep = (Pow(2,1:numGens) - Pow(1,1:numGens))';

% initialize lambda via heuristic approach
slopes = inf(numGens+1,1);
for i = 1:(numGens+1)
    slopes(i) = (c(i,numPts(i)) + c_bar(i) + h_bar(i))/q_max(i);
end
[Y,I] = sort(slopes);
slopes_sorted = slopes(I);
q_max_sorted = q_max(I);
cumsum_q_max_sorted = cumsum(q_max_sorted);

lambda_init = inf(size(x_s));
for i = 1:T
    for j = 1:quadPts
        tmp = find(cumsum_q_max_sorted-x_s(i,j)>0);
        if numel(tmp) == 0
            error('something wrong');
        end
        lambda_init(i,j) = slopes_sorted(tmp(1));
    end
end

% initialization for time-dependent lambda (not history dependent)
lambda = lambda_init*p_s;

nPts = 100; % number of function evaluation points for each generator
z_sol_avg = nan(T,numGens+2);
% store value function and on/off decisions of relaxed problem
% for each generator
V1 = cell(numGens,1);
V2 = cell(numGens,1);
V = nan(2,Lu(end)+Ld(end),T+1); % for buy / sell generators
uGen1 = cell(numGens,1);
uGen2 = cell(numGens,1);
uGen = nan(2,Lu(end)+Ld(end),T); % for buy / sell generators
L_g_cost = nan(numGens+2,1); % decomposed costs from Lagrangian relaxation
grad = nan(size(lambda)); % stochastic subgradient

% keep track of best cost information when optimizating Lagrangian
best_idx = 1;
best_cost = -inf;
best_lambda = nan(size(lambda)); 
opt_cost_vec = nan(maxIters,1);
best_L_g_cost = nan(numGens+2,1);
best_grad = nan(size(lambda));
best_V1 = cell(numGens,1);
best_uGen1 = cell(numGens,1);
best_V2 = cell(numGens,1);
best_uGen2 = cell(numGens,1);
best_V = nan(2,Lu(end)+Ld(end),T+1);
best_uGen = nan(2,Lu(end)+Ld(end),T+1);
for i = 1:numGens
    V1{i} = nan(Lu(i),nPts,T+1);
    V2{i} = nan(Ld(i),T+1);
    uGen1{i} = nan(Lu(i),nPts,T);
    uGen2{i} = nan(Ld(i),T);
    % store the value function and on/off decisions for best lambda
    best_V1{i} = nan(Lu(i),nPts,T+1);
    best_uGen1{i} = nan(Lu(i),nPts,T);
    best_V2{i} = nan(Ld(i),T+1);
    best_uGen2{i} = nan(Ld(i),T);
end

% poolobj = parpool('local',8);

time = nan(numGens+2,1);

ttot = tic;
for k = 1:maxIters
    % all generators excluding buy and sell
%     parfor i = 1:numGens
    for i = 1:numGens
        tind = tic;
        [L_g_cost(i),z_sol_avg(:,i),V1{i},uGen1{i},V2{i},uGen2{i}] = ...
            Lg_r(q_min(i),q_max(i),c_bar(i),h_bar(i),Lu(i),Ld(i),lambda,...
            T,Rd(i),Ru(i),nPts,Pow(:,i),Feval(:,i),Pstep(i));
        time(i) = toc(tind);
    end
    tind = tic;
    % for buy generator
    [L_g_cost(numGens+1),z_sol_avg(:,numGens+1),V(1,:,:),uGen(1,:,:)]...
        = Lg(q(numGens+1,1:numPts(numGens+1))',...
        c(numGens+1,1:numPts(numGens+1))',c_bar(numGens+1),...
        h_bar(numGens+1),Lu(numGens+1),Ld(numGens+1),lambda,T);
    time(numGens+1) = toc(tind);

    tind = tic; 
    % for sell generator
    [L_g_cost(numGens+2),z_sol_avg(:,numGens+2),V(2,:,:),uGen(2,:,:)]...
        = Lg(q(numGens+2,1:numPts(numGens+2))',...
        c(numGens+2,1:numPts(numGens+2))',c_bar(numGens+2),...
        h_bar(numGens+2),Lu(numGens+2),Ld(numGens+2),-lambda,T);
    z_sol_avg(:,numGens+2) = -1*z_sol_avg(:,numGens+2);
    time(numGens+2) = toc(tind); 
    
    opt_cost_vec(k) = sum(L_g_cost) + sum(lambda.*(x_s*p_s));
    
    if (opt_cost_vec(k) > best_cost)
        best_cost = opt_cost_vec(k);
        best_idx = k;
        best_lambda = lambda; 
        best_V1 = V1;
        best_uGen1 = uGen1;
        best_V2 = V2;
        best_uGen2 = uGen2;
        best_V = V; 
        best_uGen = uGen;
        best_L_g_cost = L_g_cost;
        best_grad = grad;
    end
    
    grad = x_s*p_s - sum(z_sol_avg,2);
    lambda = lambda + stepVec(k).*rho0.*grad;  
%     % some other possible choices for subgradient step
%     lambda = lambda + alpha^k.*rho0.*grad;
%     lambda = lambda + (1/k).*rho0.*grad;
%     lambda = lambda + (1/sqrt(k)).*rho0.*grad;
%     lambda = lambda + rho0.*grad;
end
tot_time=toc(ttot);

% delete(poolobj);

save(['max_Lg_r_',filename]);
