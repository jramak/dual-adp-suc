function [cost,z_sol,V1,uGen1,V2,uGen2] = Lg_r(qmin,qmax,c_bar,h_bar,...
    Lu,Ld,lambda,T,Rd_i,Ru_i,nPts,P,F,Pstep,gen_state,z_prev)
% all inputs are in column form (if they are vectors)

% number of states = Lu + Ld
% first Lu correspond to min-up time and generator being on
% Lu+1 to Lu+Ld correspond to min-down time (gen being off)
% T = 168; % time periods

% pre-compute matrices for DP runs
% pow_s = cell(Lu,1);
idx_s = cell(Lu,1);
P_mat_s = cell(Lu,1);
F_mat_s = cell(Lu,1);
for i = 1:(Lu-1)
    tmp_idx = (P <= min(qmax,qmin+i*Ru_i)+eps);
    idx_s{i} = find(tmp_idx);
    p = P(tmp_idx);
    tmp_P_mat = repmat(P,1,length(p));
    tmp_F_mat = repmat(F,1,length(p));
    for j = 1:length(p)
        % include one discretization point above p(j)+Ru
        % and also one discretization point below p(j)-Rd
        % this ensures algorithm outputs true lower bound
        tmp_idx2 = ((p(j)-Rd_i-Pstep) > P) | (P > (p(j)+Ru_i+Pstep));
        tmp_P_mat(tmp_idx2,j) = inf;
        tmp_F_mat(tmp_idx2,j) = inf;
    end
    P_mat_s{i} = tmp_P_mat;
    F_mat_s{i} = tmp_F_mat;
end
% for case when state is Lu, any previous power level is possible
idx_s{Lu} = (1:nPts)';
p = P;
tmp_P_mat = repmat(P,1,length(p));
tmp_F_mat = repmat(F,1,length(p));
for j = 1:length(p)
    tmp_idx2 = ((p(j)-Rd_i-Pstep) > P) | (P > (p(j)+Ru_i+Pstep));
    tmp_P_mat(tmp_idx2,j) = inf;
    tmp_F_mat(tmp_idx2,j) = inf;
end
P_mat_s{Lu} = tmp_P_mat;
F_mat_s{Lu} = tmp_F_mat;

% for case when generator can be turned off only if within ramp down
% contraint
can_off_flag = (P <= (qmin+Rd_i)+eps);

% store value function and policy in two appropriate parts:
% PART1
% for states l = 1, 2, ..., Lu, need previous power levels for state
% PART2
% for states l = Lu+1, Lu+2, ..., Lu+Ld, previous power level 0
% PART1
V1 = inf(Lu,nPts,T+1); % value function storage
V1(:,:,T+1) = 0; % terminal condition
uGen1 = nan(Lu,nPts,T); % store generator on/off decisions
zPow1 = nan(Lu,nPts,T); % amount of power produced
zIdx1 = nan(Lu,nPts,T); % index of amount power produced
% PART2
V2 = inf(Ld,T+1);
V2(:,T+1) = 0;
uGen2 = nan(Ld,T); % store generator on/off decisions
zPow2 = nan(Ld,T);
zIdx2 = nan(Ld,T);

% precompute for l=Lu+Ld case
tmp_idx_n = idx_s{1};
tmp_P = P(tmp_idx_n);
tmp_F = F(tmp_idx_n);
% backwards recursion
for t = T:-1:1
    % PART1
    % for states l = 1, 2, ..., Lu, need previous power levels for state

    % for l = 1, 2, 3, ..., Lu-1, generator must stay on
    for i = 1:(Lu-1)
        idx = idx_s{i};
        uGen1(i,idx,t) = 1; % generator is on
        tmp_V = repmat(V1(i+1,:,t+1)',1,length(idx));
        [V1(i,idx,t),zIdx1(i,idx,t)] = min(F_mat_s{i} - lambda(t).*P_mat_s{i} + c_bar + tmp_V);
        zPow1(i,idx,t) = P(zIdx1(i,idx,t));
    end

    % for l = Lu, generator can turn off or stay on
    idx = idx_s{Lu};
    uGen1(Lu,idx,t) = 1; % assume generator is turned on
    tmp_V = repmat(V1(Lu,:,t+1)',1,length(idx));
    [V1(Lu,idx,t),zIdx1(Lu,idx,t)] = min(F_mat_s{Lu} - lambda(t).*P_mat_s{Lu} + c_bar + tmp_V);
    zPow1(Lu,idx,t) = P(zIdx1(Lu,idx,t));
    % turn off case
    tmp_flag = (V2(1,t+1) < V1(Lu,idx,t)) & can_off_flag';
    uGen1(Lu,tmp_flag,t) = 0; % generator is turned off
    V1(Lu,tmp_flag,t) = V2(1,t+1);
    zPow1(Lu,tmp_flag,t) = 0;
    zIdx1(Lu,tmp_flag,t) = nan;

    % PART2
    % for states l = Lu+1, Lu+2, ..., Lu+Ld, previous power level 0

    % l = Lu+1, Lu+2, ..., Lu+Ld-1, generator must stay off
    uGen2(1:(Ld-1),t) = 0; % generator must stay off
    V2(1:(Ld-1),t) = V2(2:Ld,t+1);
    zPow2(1:(Ld-1),t) = 0;

    % l = Lu+Ld, generator could turn on or stay off
    % same as when first turned on: l = 1
    % turn on case:
    uGen2(Ld,t) = 1;
    [V2(Ld,t),zIdx2(Ld,t)] = min(tmp_F - lambda(t).*tmp_P + c_bar + h_bar + V1(1,tmp_idx_n,t+1)');
    zPow2(Ld,t) = tmp_P(zIdx2(Ld,t));
    % stay off if
    if V2(Ld,t+1) < V2(Ld,t)
        uGen2(Ld,t) = 0;
        V2(Ld,t) = V2(Ld,t+1);
        zPow2(Ld,t) = 0;
        zIdx2(Ld,t) = nan;
    end
end

% move forward in time to generate solution
cost = V2(Ld,1); % Lagrangian evaluation
z_sol = nan(size(lambda));
%z_idx_prev = nan;
%gen_state = Lu+Ld; % initial start state
[~, z_idx_prev] = min(abs(zPow1(gen_state,:,1) - z_prev))
myeps = 0.01;
for t = 1:T
    % PART1
    % for states l = 1, 2, ..., Lu, need previous power levels for state
    % zPow1 = nan(Lu,nPts,T,numDemands);

    % PART2
    % for states l = Lu+1, Lu+2, ..., Lu+Ld, previous power level 0
    % zPow2 = nan(Ld,T);

    if gen_state < Lu
        % generator must stay on
        z_sol(t) = zPow1(gen_state,z_idx_prev,t);
        assert(z_sol(t) > 0);
        z_idx_prev = zIdx1(gen_state,z_idx_prev,t);
        gen_state = gen_state + 1;
    elseif (abs(gen_state - Lu) <= myeps)
        % generator could stay on or turn off
        z_sol(t) = zPow1(gen_state,z_idx_prev,t);
        z_idx_prev = zIdx1(gen_state,z_idx_prev,t);
        if (abs(z_sol(t)) <= myeps)
            gen_state = gen_state + 1;
        end
    elseif ((Lu < gen_state) && (gen_state < (Lu+Ld)))
        % generator must stay off
        z_sol(t) = zPow2(gen_state-Lu,t);
        z_idx_prev = zIdx2(gen_state-Lu,t);
        gen_state = gen_state + 1;
        assert(abs(z_sol(t)) <= myeps);
    elseif (abs(gen_state - (Lu+Ld)) <= myeps)
        % generator can stay off or turn on
        z_sol(t) = zPow2(gen_state-Lu,t);
        z_idx_prev = zIdx2(gen_state-Lu,t);
        if z_sol(t) > myeps
            gen_state = 1;
        end
    else
        disp('something is wrong with the forward calculation');
    end
end
