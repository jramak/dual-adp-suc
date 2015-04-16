function [ub_cost,ub_time] = compute_ub_fv(filename)

load(['max_Lg_r_d_',filename]);
%load(filename);

% z = 2.575; % for 99%
% z = 1.96; % for 95%

% z.*z_sol_std./(sqrt(paths))
irgdx(filename,'D_s_ub');

% compute upper bound from policy
[~,numSamples] = size(D_s_ub);

% store ub cost for every sample path
cost_ub = nan(numSamples,1);
% solution for every sample path = numDemands x T x numsamplepaths
% store solution for every sample path
u_sol_mat_ub = nan(T,numGens,numSamples); 
% u_sol_mat_ub(:,(numGens+1):(numGens+2),:) = 1;
% z_sol_mat_ub = nan(T,numGens+2,numSamples); 

ub_cost = nan(T,numSamples);
ub_time = nan(T,numSamples);

% for k = 1:numSamples
tic; 
for k = 1:numSamples
    z_prev = zeros(numGens,1);
    z_prev_idx = nan(numGens,1);
    states = (Lu(1:numGens)+Ld(1:numGens)).*ones(numGens,1); 
%     for t = 1:T
    for t = 1:T
        % for the normal generators
        v_i = zeros(numGens,1);
        
%         fv_i = zeros(numGens,1); % keep track of which need future cost
%         sumV = 0;
        % generate future value function matrix
%         fV = nan(numGens,nPts); 
        Voff = zeros(numGens,1);
%         Von = nan(numGens,nPts);
        Von = zeros(numGens,nPts);
        
        for i = 1:numGens
            % do a state update
            if states(i) < Lu(i)
                % generator must stay on
                tmp = uGen1{i};
                u_sol_mat_ub(t,i,k) = tmp(states(i),z_prev_idx(i),t);
%                try
                assert(u_sol_mat_ub(t,i,k) == 1); % ensure generator on
%                catch
%                    disp('something wrong');
%                end
                % state update
                states(i) = states(i) + 1; 
                tmpV = V1{i};
                Von(i,:) = tmpV(states(i),:,t+1);
            elseif states(i) == Lu(i)
                tmp = uGen1{i};
                u_sol_mat_ub(t,i,k) = tmp(states(i),z_prev_idx(i),t);
                % state update
                if u_sol_mat_ub(t,i,k) == 0
                    states(i) = states(i) + 1;
                    tmpV = V2{i};
                    Voff(i) = tmpV(states(i)-Lu(i),t+1);
                else
%		try
                    assert(u_sol_mat_ub(t,i,k) == 1); 
                    tmpV = V1{i};
                    Von(i,:) = tmpV(states(i),:,t+1);
%		catch
%			disp(i);
%			disp(t);
%			disp(k);
%			disp(u_sol_mat_ub(t,i,k));
%			disp(states);
%			disp(z_prev_idx);
%		end
                end
            elseif ((Lu(i) < states(i)) && (states(i) < Lu(i)+Ld(i)))
                % generator must stay off
                tmp = uGen2{i};
                u_sol_mat_ub(t,i,k) = tmp(states(i)-Lu(i),t);
                assert(u_sol_mat_ub(t,i,k) == 0); % ensure generator off
                % state update
                states(i) = states(i) + 1;
                tmpV = V2{i};
                Voff(i) = tmpV(states(i)-Lu(i),t+1);
            elseif states(i) == Lu(i)+Ld(i)
                tmp = uGen2{i};
                tmpV = V2{i};
                u_sol_mat_ub(t,i,k) = tmp(states(i)-Lu(i),t);
                if u_sol_mat_ub(t,i,k) == 1
                    v_i(i) = 1; % turn on flag for generator
                    states(i) = 1;
                    tmpV = V1{i};
                    Von(i,:) = tmpV(states(i),:,t+1);
                else
                    assert(u_sol_mat_ub(t,i,k) == 0);
                    tmpV = V2{i};
                    Voff(i) = tmpV(states(i)-Lu(i),t+1);
                end
            else
                error('something wrong with states(i)');
            end
            
%             if states(i) <= Lu(i)
%                 tmpV = V1{i};
%                 fv_i(i) = 1;
%                 fV(i,:) = tmpV(states(i),:,t);
%             else
%                 tmpV = V2{i};
%                 fv_i(i) = 0; 
% %                 try 
% %                 sumV = sumV + tmpV(states(i)-Lu(i),t);
% %                 catch
% %                     disp('something is wrong');
% %                 end
%             end
        end
        % sum value function for buy and sell generators
%         sumV = sumV + V(1,1,t) + V(2,1,t);

        u_i = u_sol_mat_ub(t,:,k)';
        z_lo = [max(q_min(1:numGens).*u_i,z_prev-Rd(1:numGens)-...
            (1-u_i).*q_min(1:numGens));q_min((numGens+1):(numGens+2))];
        z_hi = [min(q_max(1:numGens).*u_i,z_prev+Ru(1:numGens)+...
            v_i.*q_min(1:numGens));q_max((numGens+1):(numGens+2))];
        
%         z_idx_lo = zeros(numGens,1);
%         z_idx_hi = zeros(numGens,1); 
        z_idx_lo = ones(numGens,1);
        z_idx_hi = ones(numGens,1); 
        for i = find(u_i)'
            tmp = find(z_lo(i) <= Pow(:,i));
            assert(~isempty(tmp));
%             z_idx_lo(i) = max(tmp(1)-1,1);
            z_idx_lo(i) = tmp(1);
            tmp = find(Pow(:,i) <= z_hi(i));
            assert(~isempty(tmp));
%             z_idx_hi(i) = min(tmp(end)+1,nPts);
            z_idx_hi(i) = tmp(end);
        end
        
%         iwgdx('ubdata.gdx','z_lo',z_lo,'z_hi',z_hi,...
%             'u_i',[u_i;1;1],'v_i',[v_i;1;1],'d',D_s_ub(t,k),...
%             'fv_i',fv_i,'sumV',sumV,'Plevels',Plevels','fV',fV,...
%             'states',[states;1;1],'z_idx_lo',z_idx_lo,'z_idx_hi',z_idx_hi);
        iwgdx('ubdata.gdx','z_lo',z_lo,'z_hi',z_hi,...
            'u_i',[u_i;1;1],'v_i',[v_i;1;1],'d',D_s_ub(t,k),...
            'Plevels',Pow','Von',Von,'Voff',Voff,...
            'states',[states;1;1],'z_idx_lo',[z_idx_lo;1;1],...
            'z_idx_hi',[z_idx_hi;1;1]);
        [status,result]=system(['gams subgrad_ub_cont_fv_pwl --probdata=',...
            filename,' --subprob=',num2str(numGens+2),' mip=gurobi lo=3']);
        irgdx('ubout.gdx');
        
        flag_rel = (u_i == 1);
        
        ub_cost(t,k) = tot_cost;
        ub_time(t,k) = time; 
        % set the power level of generators
        z_prev(~flag_rel) = 0;
        z_prev_idx(~flag_rel) = nan;
%         z_prev
        % DO BLAH BLAH BLHA
        z_prev_idx(flag_rel) = round((z_l(flag_rel) - ...
            q_min([u_i;0;0] == 1))./(Pstep(flag_rel)) + 1);
        
        % force z idx to be within bounds
        z_prev_idx(flag_rel)=min(z_prev_idx(flag_rel),z_idx_hi(flag_rel));
        z_prev_idx(flag_rel)=max(z_prev_idx(flag_rel),z_idx_lo(flag_rel));
        
        tmpP = Pow(:,flag_rel);
        z_prev(flag_rel) = tmpP(sub2ind(size(tmpP),...
            z_prev_idx(flag_rel)',1:length(z_prev_idx(flag_rel))));
        assert(all((z_prev(flag_rel)>=z_lo([flag_rel;false;false]))...
            & (z_prev(flag_rel)<=z_hi([flag_rel;false;false]))));
        
%         z_lo = max(
%         % for the buy and sell generators
%         % for buy generator
%         u_sol_mat_ub(t,numGens+1,k) = uGen(states(i),t,1);
%         % for sell generator
%         u_sol_mat_ub(t,numGens+2,k) = uGen(states(i),t,2);
        % run LP to determine generation levels
%         z_lo = 
    end
end
tot_run_time = toc;

% save(['max_Lg_r_d_',filename]);
