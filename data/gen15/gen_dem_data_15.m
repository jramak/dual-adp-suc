function gen_dem_data_15(maxDemand,percent_sig)

numOffers = 10;
numGens = 15;

s1.name = 'GenOfferPrice';
s1.form = 'full';
r_struct = rgdx('uc_all',s1);
r_price = r_struct.val;
price = r_price((numOffers+1):(numOffers+numGens),1:numOffers);

s2.name = 'GenOfferQty';
s2.form = 'full';
r_struct = rgdx('uc_all',s2);
r_qty = r_struct.val;
qty = r_qty((numOffers+1):(numOffers+numGens),1:numOffers);

s3.name = 'ColdStartCost';
s3.form = 'full';
r_struct = rgdx('uc_all',s3);
r_h_bar = r_struct.val;
% added spot buying and selling unit with 0 startup cost
h_bar = [r_h_bar((numOffers+1):(numOffers+numGens));0;0];

s4.name = 'NoLoadCost';
s4.form = 'full';
r_struct = rgdx('uc_all',s4);
r_c_bar = r_struct.val;
% added spot buying and selling unit with 0 transaction cost
c_bar = [r_c_bar((numOffers+1):(numOffers+numGens));0;0];

s5.name = 'MinDownTime';
s5.form = 'full';
r_struct = rgdx('uc_all',s5);
r_minDownTime = r_struct.val;
% added min down time of 1 for both buying and selling units
Ld = [r_minDownTime((numOffers+1):(numOffers+numGens));1;1];

s6.name = 'MinUpTime';
s6.form = 'full';
r_struct = rgdx('uc_all',s6);
r_minUpTime = r_struct.val;
% added min up time of 1 for both buying and selling units
Lu = [r_minUpTime((numOffers+1):(numOffers+numGens));1;1];

s7.name = 'xmax';
s7.form = 'full';
r_struct = rgdx('uc_all',s7);
r_qmax = r_struct.val;
qmax = r_qmax((numOffers+1):(numOffers+numGens));

s8.name = 'xmin';
s8.form = 'full';
r_struct = rgdx('uc_all',s8);
r_qmin = r_struct.val;
qmin = r_qmin((numOffers+1):(numOffers+numGens));

T = 168;
% d = 2000*ones(T,1);
% maxDemand = 2000;
% percent_sig = 0.15;
sigma = maxDemand*percent_sig; 
quadPts = 10;
[x,c] = lgwt(quadPts,maxDemand-4*sigma,maxDemand+4*sigma);

d = read_PJM_xls(maxDemand);

x_s = (d./maxDemand)*x';
c_s = (d./maxDemand)*c';
p_s = inf(size(c_s));
p_cs = inf(size(c_s)); % corrected/normalized p_s

d_check = inf(size(d));
d_check_cs = inf(size(d));
for i = 1:T
%     d_check(i) = sum(d(i).*mynormpdf(x_s(i,:),d(i),d(i)*0.05).*c_s(i,:));
    p_s(i,:) = mynormpdf(x_s(i,:),d(i),d(i)*percent_sig).*c_s(i,:);
    p_cs(i,:) = p_s(i,:)./sum(p_s(i,:)); % normalize just in case
%     d_check(i) = sum(d(i).*c_s(i,:));
    d_check(i) = sum(d(i).*p_s(i,:));
    d_check_cs(i) = sum(d(i).*p_cs(i,:));
end % not completely the same as d but that's, we can check later

% sum(p_s,2) % all the probabilities are the same. Only the demands change!
% sum(p_cs,2)
% sum(abs(d_check-d))/maxDemand
% sum(abs(d_check_cs-d))

% use normalized probability distribution
p_s = p_cs(1,:);

% generate demand 100 scenarios:
scenarios = 100;
D_idx_lb = gendist(p_s,T,scenarios);
D_idx_ub = gendist(p_s,T,scenarios);
D_s_lb = inf(size(D_idx_lb)); % demand scenarios for lower bound
D_s_ub = inf(size(D_idx_ub)); % demand scenarios for upper bound

for i = 1:scenarios
    for t = 1:T
       D_s_lb(t,i) = x_s(t,D_idx_lb(t,i));
       D_s_ub(t,i) = x_s(t,D_idx_ub(t,i));
    end
end

numPtsT = inf(numGens,1);

% check convexity (that price slopes are increasing for each generator)
for i = 1:numGens
    price_tmp = price(i,price(i,:)~=0);
    numPtsT(i) = length(price_tmp);
    assert(all(diff(price_tmp)>=0));
end

% generate points for piecewise linear function (q,c) for each generator
q_tmp = inf(numGens,numOffers+1);
c_tmp = inf(numGens,numOffers+1);

for i = 1:numGens
    q_tmp(i,1) = 0;
    c_tmp(i,1) = 0;
    for j = 2:(numPtsT(i)+1)
        q_tmp(i,j) = q_tmp(i,j-1) + qty(i,j-1);
        c_tmp(i,j) = c_tmp(i,j-1) + qty(i,j-1)*price(i,j-1);
    end
end

q = inf(numGens+2,numOffers+1);
c = inf(numGens+2,numOffers+1);
numPts = inf(numGens+2,1);
for i = 1:numGens
    tmp_idx = find(logical(-1*diff(q_tmp(i,:) <= qmin(i))));
    q(i,1:length(tmp_idx:(numOffers+1))) = q_tmp(i,tmp_idx:(numOffers+1));
    c(i,1:length(tmp_idx:(numOffers+1))) = c_tmp(i,tmp_idx:(numOffers+1));
    c(i,1) = (qmin(i) - q(i,1))*price(i,tmp_idx) + c(i,1);
    q(i,1) = qmin(i);
    numPts(i) = numPtsT(i) - tmp_idx + 2;
end
% spot market buying unit:
% q(numGens+1,1) = 0; c(numGens+1,1) = 0; 
% q(numGens+1,2) = 500; c(numGens+1,2) = 50000; % $100 slope per MW
% q(numGens+1,3) = maxDemand; c(numGens+1,3) = 50000+(maxDemand-500)*200; % $200 slope per MW
% numPts(numGens+1) = 3;
q(numGens+1,1) = 0; c(numGens+1,1) = 0; 
q(numGens+1,2) = 100; c(numGens+1,2) = 20000; % $200 slope per MW
q(numGens+1,3) = maxDemand; c(numGens+1,3) = 20000+(maxDemand-100)*850; % $850 slope per MW
numPts(numGens+1) = 3;
% spot market selling unit
% q(numGens+2,1) = -1000500; c(numGens+2,1) = -10001000; % $10 slope per MW
% q(numGens+2,2) = -500; c(numGens+2,2) = -1000; % $20 slope per MW
% q(numGens+2,3) = 0; c(numGens+2,3) = 0;

% q(numGens+2,1) = 0; c(numGens+2,1) = 0;
% q(numGens+2,2) = 500; c(numGens+2,2) = 1000; % -$20 slope per MW
% q(numGens+2,3) = 500+sum(qmax); c(numGens+2,3) = 1000+sum(qmax)*10; % -$10 slope per MW
% numPts(numGens+2) = 3;

q(numGens+2,1) = 0; c(numGens+2,1) = 0;
q(numGens+2,2) = sum(qmax); c(numGens+2,2) = sum(qmax)*10; % $10 slope per MW
                                                           % penalty for
                                                           % selling
numPts(numGens+2) = 2;

% for writing max and min generator levels
q_min = [qmin;q(numGens+1,1);q(numGens+2,1)];
q_max = [qmax;q(numGens+1,numPts(numGens+1));q(numGens+2,numPts(numGens+2))];

% for ramp up / down constraints
s9.name = 'RampDownRate';
s9.form = 'full';
r_struct = rgdx('uc_all',s9);
r_RampDownRate = r_struct.val;
% added min down time of 1 for both buying and selling units
Rd = [60.*r_RampDownRate((numOffers+1):(numOffers+numGens));...
    60*maxDemand;60*sum(qmax)];

s10.name = 'RampUpRate';
s10.form = 'full';
r_struct = rgdx('uc_all',s10);
r_RampUpRate = r_struct.val;
% added min up time of 1 for both buying and selling units
Ru = [60.*r_RampUpRate((numOffers+1):(numOffers+numGens));...
    60*maxDemand;60*sum(qmax)];

% % plotting piecewise linear functions for each generator
% for i = 1:numGens
%     figure;
%     plot(q(i,1:numPts),c(i,1:numPts));
% end

% q(gen#,demand pts), c(gen#,cost pts)
z_inc = nan(size(q_min));
nPts = 100; 
Pow = nan(nPts,length(q_min));
Feval = nan(nPts,length(q_min));
for i = 1:length(z_inc)
   [Pow(:,i),Feval(:,i)] = eval_pwl(q(i,1:numPts(i))',c(i,1:numPts(i))',q_min(i),q_max(i),nPts); 
   z_inc(i) = Pow(2,i) - Pow(1,i);
end

% x_s = 168 by 10 matrix representing for each time period possible demands
% D_s_lb (ub) = 168 by 100 matrix for lower (upper) bound runs
%             = each column represents one sample path

iwgdx(['uc15_',num2str(maxDemand),'md_',num2str(percent_sig*100),...
    'sig.gdx'],'q',q,'c',c,'numPts',numPts,'h_bar',h_bar,'c_bar',c_bar,...
    'Ld',Ld,'Lu',Lu,'d',d,'D_s_lb',D_s_lb,'D_s_idx_lb',D_idx_lb,...
    'D_s_ub',D_s_ub,'D_s_idx_ub',D_idx_ub,'x_s',x_s,'p_s',p_s,...
    'q_min',q_min,'q_max',q_max,'maxDemand',maxDemand,...
    'percent_sig',percent_sig,'Rd',Rd,'Ru',Ru,'nPts',nPts,...
    'z_inc',z_inc,'Pow',Pow,'Feval',Feval);