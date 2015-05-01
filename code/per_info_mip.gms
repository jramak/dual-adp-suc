$GDXin %probdata%
Sets
         i       generators              / 1 * %subprob% /
         t       time periods            / 1 * 168 /
         k      max num piecewise pts   / 1 * 11 /
         s       sample paths            / 1 * 500 /
         ;

Scalar
         numGens       number of generators without buy and sell generators ;

numGens = %subprob% - 2 ;

alias(t,r) ;

Set      dynk(i,k) ;

Parameter
         numPts(i)     number of valid cost pts for each generator ;
$load numPts

dynk(i,k) = YES$(ord(k) le numPts(i));

*display k ;

Variables
         u(i,t)          equals 1 if generator i is on and 0 otherwise
         v(i,t)          turn on variable
         y(i,t)          piecewise linear cost
         z(i,t)          real power from generator i in time period t
         g(i,t,k)        variables for piecewise linear cost
         zc              objective (total errors)
         ;

Binary Variables u, v ;
Positive Variable y, z, g ;

Parameter
         q(i,k)        x coordinates quantity of piecewise function ;
$load q

Parameter
         c(i,k)        y coordinates cost of piecewise function ;
$load c

Parameter
         c_bar(i)      minimum generator turn on cost ;
$load c_bar

Parameter
         h_bar(i)      generator start up cost ;
$load h_bar

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
*PWLEq1(i,t) ..         sum(dynk(i,k),(g(gi,t,k)$dynk(i,k))) =e= u(i,t);
*PWLEq2(i,t) ..         z(i,t) =e= sum(dynk(i,k),(q(i,k)$dynk(i,k))*(g(i,t,k)$dynk(i,k)));
*PWLEq3(i,t) ..         y(i,t) =e= sum(dynk(i,k),(c(i,k)$dynk(i,k))*(g(i,t,k)$dynk(i,k)));
PWLEq1(i,t) ..         sum(dynk(i,k),g(i,t,k)) =e= u(i,t);
PWLEq2(i,t) ..         z(i,t) =e= sum(dynk(i,k),q(i,k)*g(i,t,k));
PWLEq3(i,t) ..         y(i,t) =e= sum(dynk(i,k),c(i,k)*g(i,t,k));

Model UC /all/ ;

* 0.01% optimality gap and 1 hr time limit
UC.optcr = 0.0001 ;
UC.reslim = 3600 ;

* remove according to Steve because upper bound of 0 is stored
*g.up(i,t,k)$(ord(k) gt numPts(i)) = 0 ;

Parameter cost_gen(i,s) cost incurred for each generator for each sample path ;

Parameter cost_tot(s) cost incurred for each sample path ;

Parameter time keep track of total time ;

Parameter time_s(s) keep track of total time per sample path ;

Parameters optca(s), optcr(s), modelstatus(s) ;

Parameters       u_s(i,t,s)
                 v_s(i,t,s)
                 y_s(i,t,s)
                 z_s(i,t,s)
                 g_s(i,t,k,s)
                 zc_s(s)
                 ;

time = timeelapsed ;
loop(s,
         d(t) = D_s_lb(t,s) ;

         time_s(s) = timeelapsed ;

         Solve UC using mip minimizing zc ;

         time_s(s) = timeelapsed - time_s(s) ;

         cost_gen(i,s) = sum(t,y.l(i,t)+c_bar(i)*u.l(i,t)+h_bar(i)*v.l(i,t));

         u_s(i,t,s) = u.l(i,t) ;
         v_s(i,t,s) = v.l(i,t) ;
         y_s(i,t,s) = y.l(i,t) ;
         z_s(i,t,s) = z.l(i,t) ;
         g_s(i,t,k,s) = g.l(i,t,k) ;
         zc_s(s) = zc.l ;

         optca(s) = abs(UC.objest - UC.objval) ;
         optcr(s) = optca(s) / max(abs(UC.objest),abs(UC.objval)) ;
         modelstatus(s) = UC.Modelstat ;
);
time = timeelapsed - time ;

cost_tot(s) = sum(i,cost_gen(i,s)) ;

execute_unloadIdx 'lb_per_info.gdx', cost_tot, cost_gen, time, time_s, optca, optcr, modelstatus, u_s, v_s, y_s, z_s, g_s, zc_s ;
