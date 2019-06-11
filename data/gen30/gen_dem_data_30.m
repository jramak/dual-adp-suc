function gen_dem_data_30(maxDemand_mult,maxDemand,percent_sig)
maxDemand = maxDemand * maxDemand_mult;

numOffers = 10;
numGens = 30;
maxGens = 1011;

s7.name = 'xmax';
s7.form = 'full';
r_struct = rgdx('uc_all',s7);
r_qmax = r_struct.val;
qmax = r_qmax((numOffers+1):(numOffers+maxGens));

s8.name = 'xmin';
s8.form = 'full';
r_struct = rgdx('uc_all',s8);
r_qmin = r_struct.val;
qmin = r_qmin((numOffers+1):(numOffers+maxGens));

s1.name = 'GenOfferPrice';
s1.form = 'full';
r_struct = rgdx('uc_all',s1);
r_price = r_struct.val;
price = r_price((numOffers+1):(numOffers+maxGens),1:numOffers);

s2.name = 'GenOfferQty';
s2.form = 'full';
r_struct = rgdx('uc_all',s2);
r_qty = r_struct.val;
qty = r_qty((numOffers+1):(numOffers+maxGens),1:numOffers);

s5.name = 'MinDownTime';
s5.form = 'full';
r_struct = rgdx('uc_all',s5);
r_minDownTime = r_struct.val;
% added min down time of 1 for both buying and selling units
Ld = r_minDownTime((numOffers+1):(numOffers+maxGens));

s6.name = 'MinUpTime';
s6.form = 'full';
r_struct = rgdx('uc_all',s6);
r_minUpTime = r_struct.val;
% added min up time of 1 for both buying and selling units
Lu = r_minUpTime((numOffers+1):(numOffers+maxGens));

isaninteger = @(x)isfinite(x) & x==floor(x);

flag_price = ((~all(price == 0, 2)) & (~all(qty == 0, 2)) & ...
             (qmax > qmin) & (sum(price == 0, 2) > 1) & ...
             (sum(qty == 0, 2) > 1) & (Ld > 0) & (Lu > 0) & ...
             (isaninteger(Ld)) & (isaninteger(Lu)) & (sum(qty, 2) > qmin));

sum(flag_price)

% randomly select 2 peak generator, 5 mid range generators, 23 base generators
idx = find(qmax > 800 & flag_price);
length(idx)
s_idx1 = datasample(RandStream('mt19937ar','Seed',4),idx,2,'Replace',false);
idx = find(qmax > 200 & qmax <= 800 & flag_price);
length(idx)
s_idx2 = datasample(RandStream('mt19937ar','Seed',5),idx,5,'Replace',false);
idx = find(qmax <= 200 & flag_price);
length(idx)
s_idx3 = datasample(RandStream('mt19937ar','Seed',6),idx,23,'Replace',false);

s1.name = 'GenOfferPrice';
s1.form = 'full';
r_struct = rgdx('uc_all',s1);
r_price = r_struct.val;
price = r_price((numOffers+[s_idx1; s_idx2; s_idx3])',1:numOffers);

s2.name = 'GenOfferQty';
s2.form = 'full';
r_struct = rgdx('uc_all',s2);
r_qty = r_struct.val;
qty = r_qty((numOffers+[s_idx1; s_idx2; s_idx3])',1:numOffers);

s3.name = 'ColdStartCost';
s3.form = 'full';
r_struct = rgdx('uc_all',s3);
r_h_bar = r_struct.val;
% added spot buying and selling unit with 0 startup cost
h_bar = [r_h_bar((numOffers+[s_idx1; s_idx2; s_idx3])');0;0];

s4.name = 'NoLoadCost';
s4.form = 'full';
r_struct = rgdx('uc_all',s4);
r_c_bar = r_struct.val;
% added spot buying and selling unit with 0 transaction cost
c_bar = [r_c_bar((numOffers+[s_idx1; s_idx2; s_idx3])');0;0];

s5.name = 'MinDownTime';
s5.form = 'full';
r_struct = rgdx('uc_all',s5);
r_minDownTime = r_struct.val;
% added min down time of 1 for both buying and selling units
Ld = [r_minDownTime((numOffers+[s_idx1; s_idx2; s_idx3])');1;1];

s6.name = 'MinUpTime';
s6.form = 'full';
r_struct = rgdx('uc_all',s6);
r_minUpTime = r_struct.val;
% added min up time of 1 for both buying and selling units
Lu = [r_minUpTime((numOffers+[s_idx1; s_idx2; s_idx3])');1;1];

s7.name = 'xmax';
s7.form = 'full';
r_struct = rgdx('uc_all',s7);
r_qmax = r_struct.val;
qmax = r_qmax((numOffers+[s_idx1; s_idx2; s_idx3])');

s8.name = 'xmin';
s8.form = 'full';
r_struct = rgdx('uc_all',s8);
r_qmin = r_struct.val;
qmin = r_qmin((numOffers+[s_idx1; s_idx2; s_idx3])');

price(qty == 0) = 0;

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

quadPts = length(p_s); % number of quadrature points in distribution

% generate demand 100 scenarios:
scenarios = 500;
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

% ensure max generation is captured in PWL points
for i = 1:numGens
    qmax(i) = min(q(i, numPts(i)), qmax(i));
    % get rid of repeated points
    for j = 1:numPts(i)
        if (qmax(i) == q(i,j))
            numPts(i) = j;
            break;
        end
    end
end

% spot market buying unit:
% q(numGens+1,1) = 0; c(numGens+1,1) = 0;
% q(numGens+1,2) = 500; c(numGens+1,2) = 50000; % $100 slope per MW
% q(numGens+1,3) = maxDemand; c(numGens+1,3) = 50000+(maxDemand-500)*200; % $200 slope per MW
% numPts(numGens+1) = 3;
q(numGens+1,1) = 0; c(numGens+1,1) = 0;
q(numGens+1,2) = 100; c(numGens+1,2) = 20000; % $200 slope per MW
q(numGens+1,3) = maxDemand+4*sigma; c(numGens+1,3) = 20000+(maxDemand+4*sigma-100)*850; % $850 slope per MW
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
Rd = [60.*r_RampDownRate((numOffers+[s_idx1; s_idx2; s_idx3])');...
    60*maxDemand;60*sum(qmax)];

s10.name = 'RampUpRate';
s10.form = 'full';
r_struct = rgdx('uc_all',s10);
r_RampUpRate = r_struct.val;
% added min up time of 1 for both buying and selling units
Ru = [60.*r_RampUpRate((numOffers+[s_idx1; s_idx2; s_idx3])');...
    60*maxDemand;60*sum(qmax)];

% % plotting piecewise linear functions for each generator
% for i = 1:numGens
%     figure;
%     plot(q(i,1:numPts),c(i,1:numPts));
% end

% q(gen#,demand pts), c(gen#,cost pts)
% z_inc = nan(size(q_min));
nPts = 50;
Pow = nan(nPts,length(q_min));
Feval = nan(nPts,length(q_min));
for i = 1:(numGens+2)
   [Pow(:,i),Feval(:,i)] = eval_pwl(q(i,1:numPts(i))',c(i,1:numPts(i))',q_min(i),q_max(i),nPts);
%     z_inc(i) = Pow(2,i) - Pow(1,i);
end

Pstep = (Pow(2,1:numGens) - Pow(1,1:numGens))';

% come up with fewer breakpoints for use in MIP in
% perfect information lower bound and expection upper bound
% iterate through and drop breakpoint in pwl function that are not
% needed
assert(all(numPts >= 2));
assert(max(numPts)>=2);
maxPts = 2+2*(max(numPts)-2);
qn = inf(length(q_min),maxPts);
cn = inf(length(q_min),maxPts);
qn(:,1) = q(:,1);
cn(:,1) = c(:,1);
numPts_n = 1.*ones(length(q_min),1); % initialize
delta_y = Feval(2,:) - Feval(1,:);
myeps = 0.00001;
for i = 1:numGens
    for j = 3:nPts
        if ~(abs(Feval(j,i) - Feval(j-1,i) - delta_y(i)) <= myeps)
            % need to use the previous point
            numPts_n(i) = numPts_n(i) + 1;
            qn(i,numPts_n(i)) = Pow(j-1,i);
            cn(i,numPts_n(i)) = Feval(j-1,i);
            delta_y(i) = Feval(j,i) - Feval(j-1,i);
        end
    end
    numPts_n(i) = numPts_n(i) + 1;
    qn(i,numPts_n(i)) = Pow(nPts,i);
    cn(i,numPts_n(i)) = Feval(nPts,i);
end
[~,max_pwl_pts] = size(q);
qn(numGens+1,1:(min(maxPts,max_pwl_pts))) = q(numGens+1,...
    1:(min(maxPts,max_pwl_pts)));
qn(numGens+2,1:(min(maxPts,max_pwl_pts))) = q(numGens+2,...
    1:(min(maxPts,max_pwl_pts)));
cn(numGens+1,1:(min(maxPts,max_pwl_pts))) = c(numGens+1,...
    1:(min(maxPts,max_pwl_pts)));
cn(numGens+2,1:(min(maxPts,max_pwl_pts))) = c(numGens+2,...
    1:(min(maxPts,max_pwl_pts)));
numPts_n(numGens+1) = numPts(numGens+1);
numPts_n(numGens+2) = numPts(numGens+2);

% x_s = 168 by 10 matrix representing for each time period possible demands
% D_s_lb (ub) = 168 by 100 matrix for lower (upper) bound runs
%             = each column represents one sample path

filename = ['uc30_',num2str(10*maxDemand_mult),'md_',...
            num2str(percent_sig*100),'sig.gdx'];

% initial state calculation
init_T = 72;
irgdx(['init_expdem_',filename], 'u_s', 'z_s');

numGens = length(Lu) - 2; % number of generators not including buy/sell

%z_init = [z_s(1:numGens, init_T); 0; 0];
z_init = zeros(numGens+2, 1);
u_init = [u_s(1:numGens, init_T); 0; 0];
for i = 1:numGens
    [~, tmp_idx] = min(abs(Pow(:, i) - z_s(i, init_T)));
    z_init(i) = Pow(tmp_idx, i);
end
gen_state = ones(numGens+2, 1) .* (Lu + Ld); % initial start state
stay_on = zeros(numGens+2, 1);
stay_off = zeros(numGens+2, 1);

for i = 1:numGens
    for t = 1:init_T
        if (gen_state(i) < Lu(i))
            gen_state(i) = gen_state(i) + 1;
        elseif ((gen_state(i) == Lu(i)) && (u_s(i, t) == 0))
            gen_state(i) = gen_state(i) + 1;
        elseif ((Lu(i) < gen_state(i)) && (gen_state(i) < (Lu(i) + Ld(i))))
            gen_state(i) = gen_state(i) + 1;
        elseif ((gen_state(i) == (Lu(i) + Ld(i))) && (u_s(i, t) == 1))
            gen_state(i) = 1;
        end
    end
end

for i = 1:numGens
    if gen_state(i) < Lu(i)
        stay_on(i) = Lu(i) - gen_state(i);
    end
    if (gen_state(i) > Lu(i)) && (gen_state(i) < (Lu(i) + Ld(i)))
        stay_off(i) = (Lu(i) + Ld(i)) - gen_state(i);
    end
end

iwgdx(filename,'q',q,'c',c,'numPts',numPts,'h_bar',h_bar,'c_bar',c_bar,...
    'Ld',Ld,'Lu',Lu,'d',d,'D_s_lb',D_s_lb,'D_s_idx_lb',D_idx_lb,...
    'D_s_ub',D_s_ub,'D_s_idx_ub',D_idx_ub,'x',x,'x_s',x_s,'p_s',p_s,...
    'q_min',q_min,'q_max',q_max,'maxDemand',maxDemand,...
    'percent_sig',percent_sig,'Rd',Rd,'Ru',Ru,'nPts',nPts,...
    'Pow',Pow,'Feval',Feval,'numGens',numGens,'T',T,'quadPts',quadPts,...
    'Pstep',Pstep,'qn',qn,'cn',cn,'numPts_n',numPts_n,'gen_state',gen_state,...
    'z_init',z_init,'u_init',u_init,'stay_on',stay_on,'stay_off',stay_off);
