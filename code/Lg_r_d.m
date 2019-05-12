function [cost,z_sol,z_sol_std,V1,uGen1,V2,uGen2] = Lg_r_d(q,c,qmin,...
    qmax,c_bar,h_bar,Lu,Ld,ps,lambda,T,Rd_i,Ru_i,nPts,paths,P,F,...
    numDemands,Pstep,gen_state_init,z_init)
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
zPow1 = nan(Lu,nPts,T,numDemands); % amount of power produced
zIdx1 = nan(Lu,nPts,T,numDemands); % index of amount power produced
% PART2
V2 = inf(Ld,T+1);
V2(:,T+1) = 0;
uGen2 = nan(Ld,T); % store generator on/off decisions
zPow2 = nan(Ld,T,numDemands);
zIdx2 = nan(Ld,T,numDemands);

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
        V1(i,idx,t) = 0;
        for j = 1:numDemands
            [tmp_out,zIdx1(i,idx,t,j)] = min(F_mat_s{i} - lambda(t,j).*P_mat_s{i} + c_bar + tmp_V);
            V1(i,idx,t) = V1(i,idx,t) + ps(j).*tmp_out;
            zPow1(i,idx,t,j) = P(zIdx1(i,idx,t,j));
        end
    end

    % for l = Lu, generator can turn off or stay on
    idx = idx_s{Lu};
    uGen1(Lu,idx,t) = 1; % assume generator is turned on
    tmp_V = repmat(V1(Lu,:,t+1)',1,length(idx));
    V1(Lu,idx,t) = 0;
    for j = 1:numDemands
        [tmp_out,zIdx1(Lu,idx,t,j)] = min(F_mat_s{Lu} - lambda(t,j).*P_mat_s{Lu} + c_bar + tmp_V);
        V1(Lu,idx,t) = V1(Lu,idx,t) + ps(j).*tmp_out;
        zPow1(Lu,idx,t,j) = P(zIdx1(Lu,idx,t,j));
    end
    % turn off case
    tmp_flag = (V2(1,t+1) < V1(Lu,idx,t)) & can_off_flag';
    uGen1(Lu,tmp_flag,t) = 0; % generator is turned off
    V1(Lu,tmp_flag,t) = V2(1,t+1);
    zPow1(Lu,tmp_flag,t,:) = 0;
    zIdx1(Lu,tmp_flag,t,:) = nan;

    % PART2
    % for states l = Lu+1, Lu+2, ..., Lu+Ld, previous power level 0

    % l = Lu+1, Lu+2, ..., Lu+Ld-1, generator must stay off
    uGen2(1:(Ld-1),t) = 0; % generator must stay off
    V2(1:(Ld-1),t) = V2(2:Ld,t+1);
    zPow2(1:(Ld-1),t,:) = 0;

    % l = Lu+Ld, generator could turn on or stay off
    % same as when first turned on: l = 1
    % turn on case:
    uGen2(Ld,t) = 1;
    V2(Ld,t) = 0;
    for j = 1:numDemands
        [tmp_out,zIdx2(Ld,t,j)] = min(tmp_F - lambda(t,j).*tmp_P + c_bar + h_bar + V1(1,tmp_idx_n,t+1)');
	V2(Ld,t) = V2(Ld,t) + ps(j).*tmp_out;
        zPow2(Ld,t,j) = tmp_P(zIdx2(Ld,t,j));
    end
    % stay off if
    if V2(Ld,t+1) < V2(Ld,t)
        uGen2(Ld,t) = 0;
        V2(Ld,t) = V2(Ld,t+1);
        zPow2(Ld,t,:) = 0;
        zIdx2(Ld,t,:) = nan;
    end
end

% move forward in time to generate solution for gradient computation
D_idx = gendist(ps',T,paths); % demand indices T x numsamplepaths
cost = V2(Ld,1); % Lagrangian evaluation
z_sol_mat = nan(length(ps),T,paths); % solution for every sample path
                             % numDemands x T x numsamplepaths
for k = 1:paths
    gen_state = gen_state_init; % initial start state
    if gen_state > Lu
        z_idx_prev = nan;
    else
        [~, z_idx_prev] = min(abs(P - z_init));
    end
    for t = 1:T
        % PART1
        % for states l = 1, 2, ..., Lu, need previous power levels for state
        % zPow1 = nan(Lu,nPts,T,numDemands);

        % PART2
        % for states l = Lu+1, Lu+2, ..., Lu+Ld, previous power level 0
        % zPow2 = nan(Ld,T);

        if gen_state < Lu
            % generator must stay on
            z_sol_mat(:,t,k) = zPow1(gen_state,z_idx_prev,t,:);
            z_idx_prev = zIdx1(gen_state,z_idx_prev,t,D_idx(t,k));
            gen_state = gen_state + 1;
%             assert(all(z_sol_mat(:,t,k) > 0));
%             assert(~isnan(z_idx_prev));
        elseif gen_state == Lu
            % generator could stay on or turn off
            z_sol_mat(:,t,k) = zPow1(gen_state,z_idx_prev,t,:);
            z_idx_prev = zIdx1(gen_state,z_idx_prev,t,D_idx(t,k));
            if z_sol_mat(D_idx(t,k),t,k) == 0
                gen_state = gen_state + 1;
            end
        elseif ((Lu < gen_state) && (gen_state < Lu+Ld))
            % generator must stay off
            z_sol_mat(:,t,k) = zPow2(gen_state-Lu,t,:);
            z_idx_prev = zIdx2(gen_state-Lu,t,D_idx(t,k));
            gen_state = gen_state + 1;
%             assert(all(z_sol_mat(:,t,k) == 0));
%             assert(isnan(z_idx_prev));
        elseif (gen_state == (Lu+Ld))
            % generator can stay off or turn on
            z_sol_mat(:,t,k) = zPow2(gen_state-Lu,t,:);
            z_idx_prev = zIdx2(gen_state-Lu,t,D_idx(t,k));
            if z_sol_mat(D_idx(t,k),t,k) > 0
                gen_state = 1;
%                 assert(~isnan(z_idx_prev));
            else
%                 assert(zPow2(gen_state-Lu,t,D_idx(t,k)) == 0);
%                 assert(isnan(z_idx_prev));
            end
        else
            disp('something is wrong with the forward calculation');
        end
    end
end
z_sol = mean(z_sol_mat,3)';
z_sol_std = std(z_sol_mat,0,3)';
