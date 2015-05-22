function [ub_cost,ub_time,optca_mat,modelstat_mat,...
    tot_run_time] = compute_ub_fv(filename)

load(['max_Lg_r_d_',filename]);

irgdx(filename,'D_s_ub');

% compute upper bound from policy
[~,numSamples] = size(D_s_ub);

% solution for every sample path = numDemands x T x numsamplepaths
% store solution for every sample path
u_sol_mat_ub = nan(T,numGens,numSamples); 
z_sol_mat_ub = nan(T,numGens+2,numSamples); 

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
    z_prev_idx = nan(numGens,1);
    Voff = zeros(numGens,1);
    Von = zeros(numGens,nPts);
    
    for t = 1:T
        v_l = zeros(numGens,1);
        u_l = nan(numGens,1);
             
        % compute future value function for use in optimization  
        for i = 1:numGens
            if (states(i) < (Lu(i) - myeps))
                % generator must stay on in next time period
                tmpV = V1{i};
                Von(i,:) = tmpV(states(i)+1,:,t+1);
                Voff(i) = 0;

                u_l(i) = 1;
                states(i) = states(i) + 1; 
            elseif (abs(states(i) - Lu(i)) <= myeps)
                tmpV = V1{i};
                Von(i,:) = tmpV(states(i),:,t+1);
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)+1-Lu(i),t+1);
                
                tmp = uGen1{i};
                u_l(i) = tmp(states(i),z_prev_idx(i),t);

                if (abs(u_l(i)) <= myeps) % if turned off
                    % check to make sure generator can be turned off!
                    % generator can be turned off only if below ramping
                    if ((z_prev(i) - Rd(i) - q_min(i)) > 0)
                        u_l(i) = 1; 
                    else
                        states(i) = states(i) + 1;
                    end
                else
                    assert(abs(u_l(i)-1) <= myeps);
                end
            elseif (((Lu(i) + myeps) < states(i)) && ...
                    (states(i) < (Lu(i)+Ld(i)-myeps)))
                % generator must stay off
                Von(i,:) = 0; % optimization should not choose this
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)+1-Lu(i),t+1);
                
                u_l(i) = 0;

                states(i) = states(i) + 1;
            elseif (abs(states(i) - (Lu(i)+Ld(i))) <= myeps)
                tmpV = V1{i};
                Von(i,:) = tmpV(1,:,t+1);
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)-Lu(i),t+1);
                
                tmp = uGen2{i};
                u_l(i) = tmp(states(i)-Lu(i),t);
                if (abs(u_l(i) - 1) <= myeps)
                    v_l(i) = 1; % turn on flag for generator
                    states(i) = 1;
                else
                    assert(abs(u_l(i)) <= myeps);
                end
            else
                error('something wrong with states(i)');
            end
        end

        % store data and run optimization
        iwgdx('ubdata.gdx','u_prev',[u_prev;0;0],'z_prev',[z_prev;0;0],...
            'u_l',[u_l;1;1],'v_l',[v_l;1;1],'d',D_s_ub(t,k),...
            'Von',Von,'Voff',Voff);
        [status,result]=system(['gams ub_fv --probdata=',filename,...
            ' --subprob=',num2str(numGens+2),' --nPts=',num2str(nPts),' mip=cplex lo=3']);
        irgdx('ubout.gdx');
        
        assert(abs(optcr) <= myeps2);
        assert((abs(modelstatus - 8) <= myeps) | (abs(modelstatus - 1) <= myeps));
        
        ub_cost(t,k) = tot_cost;
        ub_time(t,k) = time; 
        optca_mat(t,k) = optca;
        modelstat_mat(t,k) = modelstatus;
        u_prev = u_l;
        u_sol_mat_ub(t,:,k) = u_l;
        z_sol_mat_ub(t,:,k) = z_l;
        z_prev = z_l(1:numGens);
        
        flag_rel = logical(u_l);
        flag_idx = find(flag_rel);
        assert(all(abs(z_l(~flag_rel)) <= myeps));
        
        % make sure z_prev_idx is in the right range
        z_prev_idx(~flag_rel) = nan;
        z_prev_idx(flag_rel) = round((z_prev(flag_rel) - ...
            q_min([flag_rel;false;false]))./Pstep(flag_rel) + 1);
        for i = flag_idx'
            if (Pow(z_prev_idx(i),i) > z_hi(i))
                z_prev_idx(i) = z_prev_idx(i) - 1;
                assert(Pow(z_prev_idx(i),i) <= z_hi(i));
            elseif (Pow(z_prev_idx(i),i) < z_lo(i))
                z_prev_idx(i) = z_prev_idx(i) + 1;
                assert(Pow(z_prev_idx(i),i) >= z_lo(i));
            end
        end
    end
end
tot_run_time = toc;
