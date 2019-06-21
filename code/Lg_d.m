function [cost,z_sol,V,uGen] = Lg_d(q,c,c_bar,h_bar,Lu,Ld,ps,lambda,T,...
    numDemands)
% function [cost,z_sol,V,uGen] = Lg_d(q,c,c_bar,h_bar,Lu,Ld,ps,lambda,T)
% all inputs are in column form (if they are vectors)

% number of states = Lu + Ld
% first Lu correspond to min-up time and generator being on
% Lu+1 to Lu+Ld correspond to min-down time (gen being off)
% T = 168; % time periods
V = inf(Lu+Ld,T+1); % CHANGED RECENTLY (2014/01/25) FROM NAN TO INF
V(:,T+1) = 0; % terminal condition
uGen = nan(Lu+Ld,T); % generator on/off decisions
z_sol = zeros(size(lambda));
LP_cost = inf(size(lambda));
LP_cost_vec = inf(numDemands,1);
% numDemands = length(ps);

% Aeq = [1, 0, -q'; 0, 1, -c'; 0, 0, ones(1,length(q))];
% beq = [0; 0; 1];
% lb = zeros(length(q)+2,1);

for t = T:-1:1
    for j = 1:numDemands
        [z_sol(t,j),LP_cost_vec(j)] = lp_power(lambda(t,j),q,c);
    end
    LP_cost(t) = ps'*LP_cost_vec;
    
    uGen(1:(Lu-1),t) = 1; % generator must be on
    V(1:(Lu-1),t) = LP_cost(t) + c_bar + V(2:Lu,t+1);
    if (V(Lu+1,t+1) <= (LP_cost(t) + c_bar + V(Lu,t+1)))
        uGen(Lu,t) = 0; % turn off generator
        V(Lu,t) = V(Lu+1,t+1);
    else
        uGen(Lu,t) = 1; % turn on generator
        V(Lu,t) = LP_cost(t) + c_bar + V(Lu,t+1);
    end
    uGen((Lu+1):(Lu+Ld-1),t) = 0; % generator must be turned off
    V((Lu+1):(Lu+Ld-1),t) = V((Lu+2):(Lu+Ld),t+1);
    if (V(Lu+Ld,t+1) <= (LP_cost(t) + c_bar + h_bar + V(1,t+1)))
        uGen(Lu+Ld,t) = 0; % turn off generator
        V(Lu+Ld,t) = V(Lu+Ld,t+1);
    else
        uGen(Lu+Ld,t) = 1; % turn on generator
        V(Lu+Ld,t) = LP_cost(t) + c_bar + h_bar + V(1,t+1);
    end
end

tol = 0.0001
cost = V(Lu+Ld,1);
state = Lu+Ld;
for t = 1:T
    if abs(state - (Lu+Ld)) < tol
        % generator can either be on or off!
        % solve LP and DP update here
        if (V(state,t+1) <= (LP_cost(t) + c_bar + h_bar + V(1,t+1)))
            z_sol(t,:) = 0;
        else
            state = 1;
        end
    elseif (state < Lu)
        % generator must be on
        % and next state is l+1
        state = state + 1;
    elseif abs(state - Lu) < tol
        % generator can be either on or off!
        % solve LP and DP update here
        if (V(state+1,t+1) <= (LP_cost(t) + c_bar + V(state,t+1)))
            state = state + 1;
            z_sol(t,:) = 0;
        end
    elseif (state < Lu+Ld) && (state > Lu)
        % generator must be off
        % and next state is l+1
        state = state + 1;
        z_sol(t,:) = 0;
    else
        error('something is wrong with DP updates');
    end
end
