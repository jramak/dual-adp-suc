function [cost,z_sol,V,uGen] = Lg(q,c,c_bar,h_bar,Lu,Ld,lambda,T)
% all inputs are in column form (if they are vectors)

% number of states = Lu + Ld
% first Lu correspond to min-up time and generator being on
% Lu+1 to Lu+Ld correspond to min-down time (gen being off)
% T = 168; % time periods
V = inf(Lu+Ld,T+1); 
V(:,T+1) = 0; % terminal condition
uGen = nan(Lu+Ld,T); % generator on/off decisions
z_sol = zeros(size(lambda));
LP_cost = inf(size(lambda));

for t = T:-1:1
    [z_sol(t),LP_cost(t)] = lp_power(lambda(t),q,c);
    
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

cost = V(Lu+Ld,1);
state = Lu+Ld;
myeps = 0.01;
for t = 1:T
    if (abs(state - (Lu+Ld)) <= myeps)
        % generator can either be on or off!
        % solve LP and DP update here
        if (V(state,t+1) <= (LP_cost(t) + c_bar + h_bar + V(1,t+1)))
            z_sol(t) = 0;
        else
            state = 1;
        end
    elseif (state < Lu)
        % generator must be on
        % and next state is l+1
        state = state + 1;
    elseif (abs(state - Lu) <= myeps)
        % generator can be either on or off!
        % solve LP and DP update here
        if (V(state+1,t+1) <= (LP_cost(t) + c_bar + V(state,t+1)))
            state = state + 1;
            z_sol(t) = 0;
        end
    elseif (state < (Lu+Ld)) && (state > Lu)
        % generator must be off
        % and next state is l+1
        state = state + 1;
        z_sol(t) = 0;
    else
        error('something is wrong with DP updates');
    end
end