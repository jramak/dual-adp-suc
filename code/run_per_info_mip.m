function run_per_info_mip(filename, num_scenarios)

irgdx(filename);
Lu = Lu;
numGens = length(Lu) - 2; % number of generators not including buy/sell

parfor i = 1:num_scenarios
  [status,result]=system(['gams per_info_mip_single --probdata=',...
                          filename,' --subprob=',num2str(numGens+2),...
                          ' --s=', num2str(i),...
                          ' mip=cplex o=/dev/null lo=0']);
  % modelstatus = 8 means integer solution found
  % make sure less than 0.01% optimality gap for all runs
  irgdx(['lb_per_info_', num2str(i), '.gdx']);
  assert(length(modelstatus) == 100);
  assert(all((abs(modelstatus - 8) <= 0.01) | (abs(modelstatus - 1) <= 0.01)));
  assert(all(optcr <= 0.1));
end


