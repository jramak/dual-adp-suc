$GDXin %probdata%
Sets
         i       generators              / 1 * %subprob% /
         t       time periods            / 1 * 168 /
         k       max num piecewise pts for buy sell gen  / 1 * 22 /
*         kk      max piecewise pts for fut reg gen       / 1 * 100 /
         s       sample paths            / 1 * 50 /
         ;

Scalar
         numGens       number of generators without buy and sell generators ;

numGens = %subprob% - 2 ;

alias(t,t1,t2) ;
*alias(t,t2) ;

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

*Parameter d_save(t) copy of d(t) ;
*d_save(t) = d(t) ;

*Parameter
*         D_s_ub(t,s)   demand sample paths ;
*$load D_s_ub

$GDXin

Parameter time keep track of total time ;

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
turnOnEq1(i,t)$(ord(t) le Lu(i)) .. sum(t1$(ord(t1) ge 1 and ord(t1) le ord(t)),v(i,t1))=l=u(i,t);
turnOnEq2(i,t)$(ord(t) gt Lu(i)) .. sum(t1$(ord(t1) ge ord(t)-Lu(i)+1 and ord(t1) le ord(t)),v(i,t1))=l=u(i,t);
turnOffEq1(i,t)$(ord(t) le Ld(i)) .. sum(t1$(ord(t1) ge 1 and ord(t1) le ord(t)),v(i,t1))=l=1-u(i,t-Ld(i));
turnOffEq2(i,t)$(ord(t) gt Ld(i)) .. sum(t1$(ord(t1) ge ord(t)-Ld(i)+1 and ord(t1) le ord(t)),v(i,t1))=l=1-u(i,t-Ld(i));
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
*PWLEq2b(i,t)$(ord(i) eq (numGens+1)) .. z(i,t) =e= sum(dynk(i,k),qn(i,k)*gb(i,t,k));
*PWLEq3b(i,t)$(ord(i) eq (numGens+1)) .. y(i,t) =e= sum(dynk(i,k),cn(i,k)*gb(i,t,k));
*
** for sell generator
*PWLEq1s(i,t)$(ord(i) eq (numGens+2)) .. sum(dynk(i,k),gs(i,t,k)) =e= u(i,t);
*PWLEq2s(i,t)$(ord(i) eq (numGens+2)) .. z(i,t) =e= sum(dynk(i,k),qn(i,k)*gs(i,t,k));
*PWLEq3s(i,t)$(ord(i) eq (numGens+2)) .. y(i,t) =e= sum(dynk(i,k),cn(i,k)*gs(i,t,k));

Model UC /all/ ;

* 0.01% optimality gap and 1 hr time limit
*UC.optcr = 0 ;
* 0.5% optimality gap and 1 hr time limit
UC.optcr = 0.005 ;
UC.reslim = 3600 ;

* remove according to Steve because upper bound of 0 is stored
*g.up(i,t,k)$(ord(k) gt numPts(i)) = 0 ;

*Parameter cost_gen(i,s) cost incurred for each generator for each sample path ;

*Parameter cost_tot(s) cost incurred for each sample path ;

*Parameter time_s(s) keep track of total time per sample path ;
*Parameter time_s_t(s,t) keep track of total time per sample path per time period;

*Parameters optca(t,s), optcr(t,s), modelstatus(t,s) ;

*Parameters       u_s(i,t,s)
*                 v_s(i,t,s)
*                 y_s(i,t,s)
*                 z_s(i,t,s)
*                 g_s(i,t,k,s)
*                 zc_s(s)
*                 tmp_idx(t2)
*                 tmp_cost(t2)
*                 tmp_slope(t2)
*                 tmp_z(t2)
*                 ;

* force buy and sell generators to be turned on
u.fx(i,t)$(ord(i) gt numGens) = 1 ;

*cost_gen(i,s) = 0;

time = timeelapsed ;

Solve UC using mip minimizing zc ;

*loop(s,
*         time_s(s) = timeelapsed ;
*
*         loop(t2,
*                 u.fx(i,t2) = u.l(i,t2) ;
*                 u.fx(i,t2)$(ord(i) gt numGens) = 1;
*
*                 d(t2) = D_s_ub(t2,s) ;
*
*                 Solve UC using mip minimizing zc ;
*
*                 z.fx(i,t2) = z.l(i,t2) ;
*
*                 cost_gen(i,s)=cost_gen(i,s)+y.l(i,t2)+c_bar(i)*u.l(i,t2)+h_bar(i)*v.l(i,t2) ;
*
*                 u_s(i,t2,s) = u.l(i,t2) ;
*                 v_s(i,t2,s) = v.l(i,t2) ;
*                 y_s(i,t2,s) = y.l(i,t2) ;
*                 z_s(i,t2,s) = z.l(i,t2) ;
*                 g_s(i,t2,k,s) = g.l(i,t2,k) ;
*                 zc_s(s) = zc.l ;
*
*                 optca(t2,s) = abs(UC.objest - UC.objval) ;
*                 optcr(t2,s) = optca(t2,s) / max(abs(UC.objest),abs(UC.objval)) ;
*                 modelstatus(t2,s) = UC.Modelstat ;
*         );
*
*         u.lo(i,t) = 0 ;
*         u.up(i,t) = 1 ;
*         z.lo(i,t) = 0 ;
*         z.up(i,t) = inf ;
*
*         d(t) = d_save(t) ;
*
*         time_s(s) = timeelapsed - time_s(s) ;
*);

*cost_tot(s) = sum(i,cost_gen(i,s)) ;

time = timeelapsed - time ;

scalars optca, optcr, modelstatus ;
optca = abs(UC.objest - UC.objval) ;
optcr = optca / max(abs(UC.objest),abs(UC.objval)) ;
modelstatus = UC.Modelstat ;

Parameter cost_gen(i) cost incurred for each generator ;

Parameter cost_tot cost incurred over all generators ;

cost_gen(i) = sum(t,y.l(i,t)+c_bar(i)*u.l(i,t)+h_bar(i)*v.l(i,t));
cost_tot = sum(i,cost_gen(i)) ;

Parameters       u_s(i,t)
                 v_s(i,t)
                 y_s(i,t)
                 z_s(i,t)
                 g_s(i,t,k)
                 zc_s
                 ;

u_s(i,t) = u.l(i,t) ;
v_s(i,t) = v.l(i,t) ;
y_s(i,t) = y.l(i,t) ;
z_s(i,t) = z.l(i,t) ;
g_s(i,t,k) = g.l(i,t,k) ;
zc_s = zc.l ;

*execute_unloadIdx 'ub_expdem_mip.gdx', cost_tot, cost_gen, time, time_s, optca, optcr, modelstatus, u_s, v_s, y_s, z_s, g_s, zc_s ;
execute_unloadIdx 'init_expdem_mip.gdx', cost_tot, cost_gen, time, optca, optcr, modelstatus, u_s, v_s, y_s, z_s, g_s, zc_s ;
