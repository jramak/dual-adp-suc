$GDXin %data%
Sets
         i       generators              / 1 * 11 /
         t       time periods            / 1 * 168 /
         k       max num piecewise pts   / 1 * 11 /
         s       sample paths            / 1 * 1 /
         ;

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
*         z(i,t)          real power from generator i in time period t
         zi(i,t)         integer representing generation level for gen i time t
*         zb(t)           production level for buy generator
*         zs(t)           production level for sell generator
         g(i,t,k)        variables for piecewise linear cost
         zc              objective (total errors)
         ;

Binary Variables u, v ;
*Positive Variable y, z, g ;
Positive Variables y, g ;
Integer Variable zi ;

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
         D_s(t,s)   demand sample paths ;
$load D_s

Parameter
         z_inc(i)   increment on generator power level ;
$load z_inc

Parameter
         nPts    number of discretization points for power generation ;
$load nPts

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
         rampUpEqG(i,t)        ramp up constraints for generators
         rampDownEqG(i,t)      ramp down constraints for generators
         PWLEq1(i,t)           1st set of PWL eq
         PWLEq2G(i,t)          2nd set of PWL eq for generators
         PWLEq2C(i,t)          2nd set of PWL eq for buy or sell gen
         PWLEq3(i,t)           3rd set of PWL eq
         ;

cost ..                zc=e=sum((i,t),y(i,t)+c_bar(i)*u(i,t)+h_bar(i)*v(i,t));
*demEq(t) ..            sum(i$(ord(i) ne 11),z(i,t))-z('11',t)=e=d(t);
demEq(t) ..            sum(i$(ord(i) lt 10),u(i,t)*q_min(i) + z_inc(i)*zi(i,t))+zi('10',t)-zi('11',t)=e=d(t);
turnOn1(i) ..          v(i,'1')=e=u(i,'1');
turnOn2(i,t)$(ord(t) gt 1) .. v(i,t)=g=u(i,t)-u(i,t-1);
turnOnEq1(i,t)$(ord(t) le Lu(i)) .. sum(r$(ord(r) ge 1 and ord(r) le ord(t)),v(i,r))=l=u(i,t);
turnOnEq2(i,t)$(ord(t) gt Lu(i)) .. sum(r$(ord(r) ge ord(t)-Lu(i)+1 and ord(r) le ord(t)),v(i,r))=l=u(i,t);
turnOffEq1(i,t)$(ord(t) le Ld(i)) .. sum(r$(ord(r) ge 1 and ord(r) le ord(t)),v(i,r))=l=1-u(i,t-Ld(i));
turnOffEq2(i,t)$(ord(t) gt Ld(i)) .. sum(r$(ord(r) ge ord(t)-Ld(i)+1 and ord(r) le ord(t)),v(i,r))=l=1-u(i,t-Ld(i));
rampUpEqG(i,t)$(ord(i) lt 10) .. u(i,t)*q_min(i) + z_inc(i)*zi(i,t) =l= u(i,t-1)*q_min(i) + z_inc(i)*zi(i,t-1) + Ru(i) + v(i,t)*q_min(i);
rampDownEqG(i,t)$(ord(i) lt 10) .. u(i,t-1)*q_min(i) + z_inc(i)*zi(i,t-1) - Rd(i) - (1-u(i,t))*q_min(i) =l= u(i,t)*q_min(i) + z_inc(i)*zi(i,t);
PWLEq1(i,t) ..         sum(dynk(i,k),g(i,t,k)) =e= u(i,t);
PWLEq2G(i,t)$(ord(i) lt 10) ..         u(i,t)*q_min(i) + z_inc(i)*zi(i,t) =e= sum(dynk(i,k),q(i,k)*g(i,t,k));
PWLEq2C(i,t)$(ord(i) ge 10) ..         zi(i,t) =e= sum(dynk(i,k),q(i,k)*g(i,t,k));
PWLEq3(i,t) ..         y(i,t) =e= sum(dynk(i,k),c(i,k)*g(i,t,k));

*Model UC /cost,demEq,turnOn1,turnOn2,turnOnEq1,turnOnEq2,turnOffEq1,turnOffEq2,rampUpEqG,rampUpEqG,rampUpEqC,rampUpEqC,PWLEq1,PWLEq2G,PWLEq2C,PWLEq3/ ;
Model UC / all / ;

UC.optcr = 0.0001 ;
*UC.optcr = 0 ;
*UC.prioropt = 2 ;

* remove according to Steve because upper bound of 0 is stored
*g.up(i,t,k)$(ord(k) gt numPts(i)) = 0 ;

Parameter cost_gen(i,s) cost incurred for each generator for each sample path ;

Parameter cost_tot(s) cost incurred for each sample path ;

Parameter time keep track of total time ;

Parameter time_s(s) keep track of total time per sample path ;

Parameters       u_s(i,t,s)
                 v_s(i,t,s)
                 y_s(i,t,s)
                 z_s(i,t,s)
                 g_s(i,t,k,s)
                 zc_s(s)
                 ;

zi.lo(i,t) = 0;
zi.up(i,t) = nPts - 1;
zi.prior('10',t) = inf;
zi.prior('11',t) = inf;
zi.lo('10',t) = q_min('10');
zi.up('10',t) = q_max('10');
zi.lo('11',t) = q_min('11');
zi.up('11',t) = q_max('11');

time = timeelapsed ;
loop(s,
         d(t) = D_s(t,s) ;

         time_s(s) = timeelapsed ;

         Solve UC using mip minimizing zc ;

         time_s(s) = timeelapsed - time_s(s) ;

         cost_gen(i,s) = sum(t,y.l(i,t)+c_bar(i)*u.l(i,t)+h_bar(i)*v.l(i,t));

         u_s(i,t,s) = u.l(i,t) ;
         v_s(i,t,s) = v.l(i,t) ;
         y_s(i,t,s) = y.l(i,t) ;
         z_s(i,t,s) = zi.l(i,t) ;
         g_s(i,t,k,s) = g.l(i,t,k) ;
         zc_s(s) = zc.l ;
);
time = timeelapsed - time ;

cost_tot(s) = sum(i,cost_gen(i,s)) ;
