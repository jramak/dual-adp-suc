$GDXin %probdata%
Sets
         i       generators              / 1 * %subprob% /
         t       time periods            / 1 * 168 /
         k       max num piecewise pts for buy sell gen  / 1 * 22 /
*         kk      max piecewise pts for fut reg gen       / 1 * 100 /
         s       sample paths            / 1 * 100 /
         ;

Scalar
         numGens       number of generators without buy and sell generators ;

numGens = %subprob% - 2 ;

alias(t,r) ;

Set      dynk(i,k) ;

Parameter
         numPts_n(i)     number of valid cost pts for each generator ;
$load numPts_n

dynk(i,k) = YES$(ord(k) le numPts_n(i));

*display k ;

Variables
         u(i,t)          equals 1 if generator i is on and 0 otherwise
         v(i,t)          turn on variable
         y(i,t)          piecewise linear cost
         z(i,t)          real power from generator i in time period t
*         g(i,t,kk)       variables for piecewise linear cost
         g(i,t,k)       variables for piecewise linear cost
         zc              objective (total errors)
         ;

Binary Variables u, v ;
Positive Variable y, z, g ;

*Parameter
*         q(i,k)        x coordinates quantity of piecewise function ;
*$load q
*
*Parameter
*         c(i,k)        y coordinates cost of piecewise function ;
*$load c

Parameter
         qn(i,k)        x coordinates quantity of piecewise function ;
$load qn

Parameter
         cn(i,k)        y coordinates cost of piecewise function ;
$load cn

Parameter
         c_bar(i)      minimum generator turn on cost ;
$load c_bar

Parameter
         h_bar(i)      generator start up cost ;
$load h_bar

*Parameter
*         Pow(kk,i)     possible power levels for each generator ;
*$load Pow
*
*Parameter
*         Feval(kk,i)   possible cost function vals for generator ;
*$load Feval

Parameter
         Lu(i)   generator minimum up time ;
$load Lu

Parameter
         Ld(i)   generator minimum down time ;
$load Ld

Parameter
         Ru(i)   generator ramp up rate ;
$load Ru

Parameter
         Rd(i)   generator ramp down rate ;
$load Rd

Parameter
         q_min(i)   generator min power level ;
$load q_min

Parameter
         q_max(i)   generator max power level ;
$load q_max

Parameter
         d(t)   demand over time t ;
$load d

Parameter
         D_s_lb(t,s)   demand sample paths ;
$load D_s_lb

$GDXin

Equations
         cost                  objective function
         demEq(t)              demand satisfaction constraint
         turnOn1(i)            initialize turn on variable v
         turnOn2(i,t)          turn variable constraints
         turnOnEq1(i,t)        first set of turn on inequalities
         turnOnEq2(i,t)        second set of turn on inequalities
         turnOffEq1(i,t)       first set of turn off inequalities
         turnOffEq2(i,t)       second set of turn off inequalities
         rampUpEq(i,t)         ramp up constraints
         rampDownEq(i,t)       ramp down constraints

         PWLEq1(i,t)           1st set of PWL eq
         PWLEq2(i,t)           2nd set of PWL eq
         PWLEq3(i,t)           3rd set of PWL eq

*         PWLEq1b(i,t)          1st set of PWL eq for buy gen
*         PWLEq2b(i,t)          2nd set of PWL eq for buy gen
*         PWLEq3b(i,t)          3rd set of PWL eq for buy gen
*
*         PWLEq1s(i,t)          1st set of PWL eq for sell gen
*         PWLEq2s(i,t)          2nd set of PWL eq for sell gen
*         PWLEq3s(i,t)          3rd set of PWL eq for sell gen
         ;

cost ..                zc=e=sum((i,t),y(i,t)+c_bar(i)*u(i,t)+h_bar(i)*v(i,t));
demEq(t) ..            sum(i$(ord(i) ne %subprob%),z(i,t))-z('%subprob%',t)=e=d(t);
turnOn1(i) ..          v(i,'1')=e=u(i,'1');
turnOn2(i,t)$(ord(t) gt 1) .. v(i,t)=g=u(i,t)-u(i,t-1);
turnOnEq1(i,t)$(ord(t) le Lu(i)) .. sum(r$(ord(r) ge 1 and ord(r) le ord(t)),v(i,r))=l=u(i,t);
turnOnEq2(i,t)$(ord(t) gt Lu(i)) .. sum(r$(ord(r) ge ord(t)-Lu(i)+1 and ord(r) le ord(t)),v(i,r))=l=u(i,t);
turnOffEq1(i,t)$(ord(t) le Ld(i)) .. sum(r$(ord(r) ge 1 and ord(r) le ord(t)),v(i,r))=l=1-u(i,t-Ld(i));
turnOffEq2(i,t)$(ord(t) gt Ld(i)) .. sum(r$(ord(r) ge ord(t)-Ld(i)+1 and ord(r) le ord(t)),v(i,r))=l=1-u(i,t-Ld(i));
rampUpEq(i,t) .. z(i,t) =l= z(i,t-1) + Ru(i) + v(i,t)*q_min(i);
rampDownEq(i,t) .. z(i,t-1) - Rd(i) - (1-u(i,t))*q_min(i) =l= z(i,t);

*PWLEq1(i,t)$(ord(i) le numGens) .. sum(kk,g(i,t,kk)) =e= u(i,t);
*PWLEq2(i,t)$(ord(i) le numGens) .. z(i,t) =e= sum(kk,Pow(kk,i)*g(i,t,kk));
*PWLEq3(i,t)$(ord(i) le numGens) .. y(i,t) =e= sum(kk,Feval(kk,i)*g(i,t,kk));

PWLEq1(i,t) .. sum(dynk(i,k),g(i,t,k)) =e= u(i,t);
PWLEq2(i,t) .. z(i,t) =e= sum(dynk(i,k),qn(i,k)*g(i,t,k));
PWLEq3(i,t) .. y(i,t) =e= sum(dynk(i,k),cn(i,k)*g(i,t,k));

** for buy generator
*PWLEq1b(i,t)$(ord(i) eq (numGens+1)) .. sum(dynk(i,k),gb(i,t,k)) =e= u(i,t);
*PWLEq2b(i,t)$(ord(i) eq (numGens+1)) .. z(i,t) =e= sum(dynk(i,k),q(i,k)*gb(i,t,k));
*PWLEq3b(i,t)$(ord(i) eq (numGens+1)) .. y(i,t) =e= sum(dynk(i,k),c(i,k)*gb(i,t,k));
*
** for sell generator
*PWLEq1s(i,t)$(ord(i) eq (numGens+2)) .. sum(dynk(i,k),gs(i,t,k)) =e= u(i,t);
*PWLEq2s(i,t)$(ord(i) eq (numGens+2)) .. z(i,t) =e= sum(dynk(i,k),q(i,k)*gs(i,t,k));
*PWLEq3s(i,t)$(ord(i) eq (numGens+2)) .. y(i,t) =e= sum(dynk(i,k),c(i,k)*gs(i,t,k));


Model UC /all/ ;

* 0.01% optimality gap and 1 hr time limit
* 0.5% optimality gap and 1 hr time limit
UC.optcr = 0.005 ;
UC.reslim = 3600 ;
UC.threads = 1 ;

* remove according to Steve because upper bound of 0 is stored
*g.up(i,t,k)$(ord(k) gt numPts(i)) = 0 ;

Parameter cost_gen(i) cost incurred for each generator for each sample path ;

Parameter cost_tot cost incurred for each sample path ;

Parameter objEst_lb, objVal_ub ;

Parameter time keep track of total time ;

*Parameter time_s(s) keep track of total time per sample path ;

Parameters optca, optcr, modelstatus ;

Parameters       u_s(i,t)
                 v_s(i,t)
                 y_s(i,t)
                 z_s(i,t)
                 g_s(i,t,k)
                 zc_s
                 ;

time = timeelapsed ;

*start solve
d(t) = D_s_lb(t,'%scenario%') ;

*time_s(s) = timeelapsed ;

Solve UC using mip minimizing zc ;

*time_s(s) = timeelapsed - time_s(s) ;

cost_gen(i) = sum(t,y.l(i,t)+c_bar(i)*u.l(i,t)+h_bar(i)*v.l(i,t));

u_s(i,t) = u.l(i,t) ;
v_s(i,t) = v.l(i,t) ;
y_s(i,t) = y.l(i,t) ;
z_s(i,t) = z.l(i,t) ;
g_s(i,t,k) = g.l(i,t,k) ;
zc_s = zc.l ;

objEst_lb = UC.objest ;
objVal_ub = UC.objval ;

optca = abs(UC.objest - UC.objval) ;
optcr = optca / max(abs(UC.objest),abs(UC.objval)) ;
modelstatus = UC.Modelstat ;

time = timeelapsed - time ;

cost_tot = sum(i,cost_gen(i)) ;

execute_unloadIdx 'lb_per_info_%scenario%.gdx', objEst_lb, objVal_ub, cost_tot,
cost_gen, time, optca, optcr, modelstatus, u_s, v_s, y_s, z_s,
g_s, zc_s ;
