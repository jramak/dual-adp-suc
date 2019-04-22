function scaledDemand = read_PJM_xls(maxDemand)

[NUMR,TXT,RAW] = xlsread('2013-day-ahead-demand-bid.xls','PJM');

% only use summer months: June, July, August
NUMR = NUMR(152:243, :);

day_labelR = inf(length(NUMR),1);
day_labelR(1:5) = 3:7; % data starts on a Tuesday!
idx = 6:length(day_labelR);
day_labelR(idx) = mod(idx - 6,7) + 1;

% remove bad values! (NaNs)
NUM = NUMR(~isnan(sum(NUMR,2)),:);
day_label = day_labelR(~isnan(sum(NUMR,2)));

% compute average demands
T = 168;
demandAvg = inf(T,1);
demandStd = inf(T,1);
tmp_idx = 1;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==1,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==1,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==2,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==2,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==3,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==3,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==4,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==4,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==5,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==5,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==6,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==6,:),1)';
tmp_idx = tmp_idx + 24;
demandAvg(tmp_idx:tmp_idx+23) = mean(NUM(day_label==7,:),1)';
demandStd(tmp_idx:tmp_idx+23) = std(NUM(day_label==7,:),1)';

assert(~any(isnan(demandAvg)));
assert(~any(isinf(demandAvg)));
assert(~any(isnan(demandStd)));
assert(~any(isinf(demandStd)));

% maxDemand = 4000;

scaledDemand = demandAvg*maxDemand/max(demandAvg);

% plot(1:T,scaledDemand);
