tol = 0.0001
% 15 gen problem
mat_res = dir('./gen15/lb_bounds/lb_*');
num = length(mat_res);
m = [4, 6, 8];
sig = [15, 20, 25];

lb_bounds_15 = nan(num, 9);
solve_times_15 = nan(num, 6);
for idx_i = 1:num
  filename = mat_res(idx_i).name;
  cd(['gen15/lb_bounds/',filename]);
  m_v = str2num(filename(9));
  sig_v = str2num(filename(13:14));
  m_idx = find(abs(m - m_v) < tol);
  sig_idx = find(abs(sig - sig_v) < tol);
  lin_idx = 3 * (m_idx - 1) + sig_idx;
  lb_bounds_15(lin_idx, 1) = m_v * 0.1;
  lb_bounds_15(lin_idx, 2) = sig_v * 0.01;
  lb_bounds_15(lin_idx, 3) = 15;
  solve_times_15(lin_idx, 1) = m_v * 0.1;
  solve_times_15(lin_idx, 2) = sig_v * 0.01;
  solve_times_15(lin_idx, 3) = 15;

  load_tmp = dir('max_Lg_r_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_15(lin_idx, 4) = best_cost / 1e6;
    solve_times_15(lin_idx, 4) = tot_time / 60;
  end

  load_tmp = dir('max_Lg_r_d_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_15(lin_idx, 5) = best_cost / 1e6;
    solve_times_15(lin_idx, 5) = tot_time / 60;
  end

  load_tmp = dir('per_info_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_15(lin_idx, 6) = mean(obj_lb) / 1e6;
    lb_bounds_15(lin_idx, 7) = 1.96 * std(obj_lb) / (1e6 * sqrt(length(obj_lb)));
    solve_times_15(lin_idx, 6) = tot_time / 60;
  end
  cd ..
  cd ..
  cd ..

  ub_filename = ['ub_fv_mip_uc15_', num2str(m_v), 'md_', num2str(sig_v), 'sig'];
  cd(['gen15/ub_bounds/',ub_filename]);

  load_tmp = dir('ub_fv_*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    ub_bounds = sum(ub_cost);
    lb_bounds_15(lin_idx, 8) = mean(ub_bounds) / 1e6;
    lb_bounds_15(lin_idx, 9) = 1.96 * std(ub_bounds) / (1e6 * sqrt(length(ub_bounds)));
  end
  cd ..
  cd ..
  cd ..
end

% 30 gen problem
mat_res = dir('./gen30/lb_bounds/lb_*');
num = length(mat_res);

lb_bounds_30 = nan(num, 9);
solve_times_30 = nan(num, 6);
for idx_i = 1:num
  filename = mat_res(idx_i).name;
  cd(['gen30/lb_bounds/',filename]);
  m_v = str2num(filename(9));
  sig_v = str2num(filename(13:14));
  m_idx = find(abs(m - m_v) < tol);
  sig_idx = find(abs(sig - sig_v) < tol);
  lin_idx = 3 * (m_idx - 1) + sig_idx;
  lb_bounds_30(lin_idx, 1) = m_v * 0.1;
  lb_bounds_30(lin_idx, 2) = sig_v * 0.01;
  lb_bounds_30(lin_idx, 3) = 30;
  solve_times_30(lin_idx, 1) = m_v * 0.1;
  solve_times_30(lin_idx, 2) = sig_v * 0.01;
  solve_times_30(lin_idx, 3) = 30;

  load_tmp = dir('max_Lg_r_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_30(lin_idx, 4) = best_cost / 1e6;
    solve_times_30(lin_idx, 4) = tot_time / 60;
  end

  load_tmp = dir('max_Lg_r_d_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_30(lin_idx, 5) = best_cost / 1e6;
    solve_times_30(lin_idx, 5) = tot_time / 60;
  end

  load_tmp = dir('per_info_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_30(lin_idx, 6) = mean(obj_lb) / 1e6;
    lb_bounds_30(lin_idx, 7) = 1.96 * std(obj_lb) / (1e6 * sqrt(length(obj_lb)));
    solve_times_30(lin_idx, 6) = tot_time / 60;
  end
  cd ..
  cd ..
  cd ..

  ub_filename = ['ub_fv_mip_uc30_', num2str(m_v), 'md_', num2str(sig_v), 'sig'];
  cd(['gen30/ub_bounds/',ub_filename]);

  load_tmp = dir('ub_fv_*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    ub_bounds = sum(ub_cost);
    lb_bounds_30(lin_idx, 8) = mean(ub_bounds) / 1e6;
    lb_bounds_30(lin_idx, 9) = 1.96 * std(ub_bounds) / (1e6 * sqrt(length(ub_bounds)));
  end
  cd ..
  cd ..
  cd ..
end

% 50 gen problem
mat_res = dir('./gen50/lb_bounds/lb_*');
num = length(mat_res);

lb_bounds_50 = nan(num, 9);
solve_times_50 = nan(num, 6);
for idx_i = 1:num
  filename = mat_res(idx_i).name;
  cd(['gen50/lb_bounds/',filename]);
  m_v = str2num(filename(9));
  sig_v = str2num(filename(13:14));
  m_idx = find(abs(m - m_v) < tol);
  sig_idx = find(abs(sig - sig_v) < tol);
  lin_idx = 3 * (m_idx - 1) + sig_idx;
  lb_bounds_50(lin_idx, 1) = m_v * 0.1;
  lb_bounds_50(lin_idx, 2) = sig_v * 0.01;
  lb_bounds_50(lin_idx, 3) = 50;
  solve_times_50(lin_idx, 1) = m_v * 0.1;
  solve_times_50(lin_idx, 2) = sig_v * 0.01;
  solve_times_50(lin_idx, 3) = 50;

  load_tmp = dir('max_Lg_r_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_50(lin_idx, 4) = best_cost / 1e6;
    solve_times_50(lin_idx, 4) = tot_time / 60;
  end

  load_tmp = dir('max_Lg_r_d_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_50(lin_idx, 5) = best_cost / 1e6;
    solve_times_50(lin_idx, 5) = tot_time / 60;
  end

  load_tmp = dir('per_info_uc*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    lb_bounds_50(lin_idx, 6) = mean(obj_lb) / 1e6;
    lb_bounds_50(lin_idx, 7) = 1.96 * std(obj_lb) / (1e6 * sqrt(length(obj_lb)));
    solve_times_50(lin_idx, 6) = tot_time / 60;
  end
  cd ..
  cd ..
  cd ..

  ub_filename = ['ub_fv_mip_uc50_', num2str(m_v), 'md_', num2str(sig_v), 'sig'];
  cd(['gen50/ub_bounds/',ub_filename]);

  load_tmp = dir('ub_fv_*.mat');
  if ~(abs(length(load_tmp)) < tol)
    assert(abs(length(load_tmp) - 1) < tol);
    load(load_tmp.name);
    ub_bounds = sum(ub_cost);
    lb_bounds_50(lin_idx, 8) = mean(ub_bounds) / 1e6;
    lb_bounds_50(lin_idx, 9) = 1.96 * std(ub_bounds) / (1e6 * sqrt(length(ub_bounds)));
  end
  cd ..
  cd ..
  cd ..
end

lb_bounds = [lb_bounds_15; lb_bounds_30; lb_bounds_50]
solve_times = [solve_times_15; solve_times_30; solve_times_50]

input.data = lb_bounds;
input.dataFormat = {'%.1f',1,'%.2f',1,'%.0f',1,'%.2f',6};
input.tableColLabels = {'c','sig','# Gen','Lg','Lg_d','mu_PInfo','HW_PInfo','mean','HW'};
input.dataNanString = '-';
latex = latexTable(input);

input.data = solve_times;
input.dataFormat = {'%.1f',1,'%.2f',1,'%.0f',1,'%.1f',3};
input.tableColLabels = {'c','sig','# Gen','Lg','Lg_d','PInfo'};
input.dataNanString = '-';
latex = latexTable(input);
