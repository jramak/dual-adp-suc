function [ub_cost,ub_time,optca_mat,modelstat_mat,...
    tot_run_time] = compute_ub_fv_mip(filename)

load(['max_Lg_r_d_',filename]);

irgdx(filename,'D_s_ub','D_s_idx_ub');

% compute upper bound from policy
[~,numSamples] = size(D_s_ub);

% solution for every sample path = numDemands x T x numsamplepaths
% store solution for every sample path
u_sol_mat_ub = nan(T,numGens,numSamples); 
z_sol_mat_ub = nan(T,numGens+2,numSamples); 

% store ub cost and optimization time
ub_cost = nan(T,numSamples);
ub_time = nan(T,numSamples);
optca_mat = nan(T,numSamples);
modelstat_mat = nan(T,numSamples);

% get rid of infinities in value function by setting to a large number
for i = 1:numGens
    tmpV1 = V1{i};
    tmpV1(isinf(tmpV1)) = 1e50;
    V1{i} = tmpV1; 
    
    tmpV2 = V2{i};
    tmpV2(isinf(tmpV2)) = 1e50;
    V2{i} = tmpV2; 
end

myeps = 0.01;
myeps2 = 0.1;

tic; 
for k = 1:numSamples
    % initialize all states to be Lu+Ld
    states = (Lu(1:numGens)+Ld(1:numGens)).*ones(numGens,1); 
    u_prev = zeros(numGens,1);
    z_prev = zeros(numGens,1);
    Voff = zeros(numGens,1);
    Von = zeros(numGens,nPts);
    
    for t = 1:T
        % reset lower and upper bounds on generator on/off decisions
        u_fx_lo = zeros(numGens,1);
        u_fx_hi = ones(numGens,1);
        
        % compute future value function for use in optimization
        % and also set gen on/off constraint
        for i = 1:numGens
            if (states(i) < (Lu(i) - myeps))
                % generator must stay on in next period
                u_fx_lo(i) = 1;
                u_fx_hi(i) = 1;
                
                tmpV = V1{i};
                Von(i,:) = tmpV(states(i)+1,:,t+1);
                Voff(i) = 0; % optimization should not choose this
            elseif (abs(states(i) - Lu(i)) <= myeps)
                tmpV = V1{i};
                Von(i,:) = tmpV(states(i),:,t+1);
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)+1-Lu(i),t+1);
            elseif (((Lu(i) + myeps) < states(i)) && ...
                    (states(i) < (Lu(i)+Ld(i)-myeps)))
                u_fx_lo(i) = 0;
                u_fx_hi(i) = 0;
                
                Von(i,:) = 0; % optimization should not choose this
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)+1-Lu(i),t+1);
            elseif (abs(states(i) - (Lu(i)+Ld(i))) <= myeps)
                tmpV = V1{i};
                Von(i,:) = tmpV(1,:,t+1);
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)-Lu(i),t+1);
            else
                error('something wrong with states(i)');
            end
        end
        
        % store data and run optimization
        iwgdx('ubdata.gdx','u_prev',[u_prev;0;0],'z_prev',[z_prev;0;0],...
            'u_fx_lo',[u_fx_lo;1;1],'u_fx_hi',[u_fx_hi;1;1],'Von',Von,...
            'Voff',Voff,'dj',x_s(t,:)');
        [status,result]=system(['gams ub_fv_mip --probdata=',...
            filename,' --subprob=',num2str(numGens+2),' --nPts=',num2str(nPts),' mip=cplex lo=3']);
        irgdx('ubout.gdx');
        
        assert(abs(optcr) <= myeps2);
        assert((abs(modelstatus - 8) <= myeps) | (abs(modelstatus - 1) <= myeps));
        
        % store results from optimization
        ub_cost(t,k) = tot_cost(D_s_idx_ub(t,k));
        ub_time(t,k) = time; 
        optca_mat(t,k) = optca;
        modelstat_mat(t,k) = modelstatus;
        u_prev = u_l(1:numGens);
        u_sol_mat_ub(t,:,k) = u_l(1:numGens);
        
        flag_rel = logical(u_l(1:numGens));
        
        % set previous power level of generators
        z_prev(~flag_rel) = 0;
        z_prev(flag_rel) = z_l([flag_rel;false;false],D_s_idx_ub(t,k));
        z_sol_mat_ub(t,:,k) = z_l(:,D_s_idx_ub(t,k)); 
        
        % do a state update for each generator    
        for i = 1:numGens
            if (states(i) < (Lu(i) - myeps))
                % check to make sure generator stayed on
                assert(abs(u_l(i)-1) <= myeps);
                states(i) = states(i) + 1;
            elseif (abs(states(i) - Lu(i)) <= myeps)
                if (abs(u_l(i)) <= myeps)
                    states(i) = states(i) + 1;
                else
                    assert(abs(u_l(i)-1) <= myeps);
                end
            elseif (((Lu(i) + myeps) < states(i)) && ...
                    (states(i) < (Lu(i)+Ld(i)-myeps)))
                % generator must have stayed off
                assert(abs(u_l(i)) <= myeps);
                states(i) = states(i) + 1; 
            elseif (abs(states(i) - (Lu(i)+Ld(i))) <= myeps)
                if (abs(u_l(i) - 1) <= myeps)
                    states(i) = 1;
                else
                    assert(abs(u_l(i)) <= myeps);
                end
            else
                error('something wrong with states(i)');
            end
        end       
    end
end
tot_run_time = toc;
