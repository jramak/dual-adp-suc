function lb_main(filename,maxIters,rho0,stepVec)

poolobj = parpool('local',8);

[~,~,~,~,~,~,~] = max_Lg_r(filename,maxIters,rho0,stepVec);

[~,~,~,~,~,~,~,~] = max_Lg_r_d(filename,maxIters,rho0,stepVec);

delete(poolobj);

irgdx(filename);
Lu = Lu;
numGens = length(Lu) - 2; % number of generators not including buy/sell

[status,result]=system(['gams per_info_mip --probdata=',...
            filename,' --subprob=',num2str(numGens+2),' mip=cplex o=/dev/null lo=0']);
        
% modelstatus = 8 means integer solution found
% make sure less than 0.01% optimality gap for all runs
irgdx('lb_per_info.gdx');
assert(all(abs(modelstatus - 8) <= 0.01));
assert(all(optcr <= 0.0002));
