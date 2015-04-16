Sets
         i       generators                              / 1 * %subprob% /
         k       max num piecewise pts                   / 1 * 11 /
         kk      max piecewise pts for fut val func      / 1 * 100 /
         ;

Scalar
         numGens       number of generators without buy and sell generators ;

numGens = %subprob% - 2 ;

Variables
         u(i)            generator on off decisions
         v(i)            turn on variable
         y(i)            piecewise linear cost
         yy(i)           piecewise linear cost for future value function
         z(i)            real power from generator i
         g(i,kk)         variables for piecewise linear cost
         gb(k)           pwl variable for buy generator
         gs(k)           pwl variable for sell generator
         zc              objective
         ;

Positive Variable y, z, g, gb, gs ;
Binary Variable u, v ;

$GDXin %probdata%
Set      dynk(i,k) ;

Parameter
         numPts(i)     number of valid cost pts for each generator ;
$load numPts

dynk(i,k) = YES$(ord(k) le numPts(i));

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
         Pow(kk,i)     possible power levels for each generator ;
$load Pow

Parameter
         Feval(kk,i)   possible cost function vals for generator ;
$load Feval

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

$GDXin

$GDXIN ubdata.gdx
Parameter
         u_prev(i)     on off generator decisions  ;
$load u_prev

Parameter
         z_prev(i)     previous power level to enforce ramping constraints  ;
$load z_prev

Parameter
         Voff(i)         future cost incurred when generator off  ;
$onUNDF
$load Voff

Parameter
         Von(i,kk)       future cost incurred when generator on ;
$onUNDF
$load Von

Parameter
         u_l(i)     on off generator decisions  ;
$load u_l

Parameter
         v_l(i)     turn on generator decisions  ;
$load v_l

Scalar
         d          observed demand that needs to be satisfied ;
$load d

$GDXIN

Scalar time keep track of total time ;

Equations
         cost                  objective function
         demEq                 demand satisfaction constraint
         turnOn(i)             equation to set turn on variable
         rampUpEq(i)           ramp up constraint
         rampDownEq(i)         ramp down constraint

         PWLEq1(i)             1st set of PWL eq
         PWLEq2(i)             2nd set of PWL eq
         PWLEq3(i)             3rd set of PWL eq
         PWLEq4(i)             4th set of PWL eq

         PWLEq1b(i)            1st set of PWL eq for buy gen
         PWLEq2b(i)            2nd set of PWL eq for buy gen
         PWLEq3b(i)            3rd set of PWL eq for buy gen

         PWLEq1s(i)            1st set of PWL eq for sell gen
         PWLEq2s(i)            2nd set of PWL eq for sell gen
         PWLEq3s(i)            3rd set of PWL eq for sell gen
         ;

cost ..              zc=e=sum(i,y(i)+c_bar(i)*u(i)+h_bar(i)*v(i))+sum(i$(ord(i) le numGens),yy(i)+(1-u(i))*Voff(i)) ;
demEq ..             sum(i$(ord(i) ne %subprob%),z(i))-z('%subprob%')=e=d;
turnOn(i) ..         v(i)=g=u(i)-u_prev(i) ;
rampUpEq(i)$(ord(i) le numGens) .. z(i) =l= z_prev(i) + Ru(i) + v(i)*q_min(i);
rampDownEq(i)$(ord(i) le numGens) .. z_prev(i) - Rd(i) - (1-u(i))*q_min(i) =l= z(i);

PWLEq1(i)$(ord(i) le numGens) .. sum(kk,g(i,kk)) =e= u(i);
PWLEq2(i)$(ord(i) le numGens) .. z(i) =e= sum(kk,Pow(kk,i)*g(i,kk));
PWLEq3(i)$(ord(i) le numGens) .. y(i) =e= sum(kk,Feval(kk,i)*g(i,kk));
PWLEq4(i)$(ord(i) le numGens) .. yy(i) =e= sum(kk,Von(i,kk)*g(i,kk));

* for buy generator
PWLEq1b(i)$(ord(i) eq (numGens+1)) .. sum(dynk(i,k),gb(k)) =e= u(i);
PWLEq2b(i)$(ord(i) eq (numGens+1)) .. z(i) =e= sum(dynk(i,k),q(i,k)*gb(k));
PWLEq3b(i)$(ord(i) eq (numGens+1)) .. y(i) =e= sum(dynk(i,k),c(i,k)*gb(k));

* for sell generator
PWLEq1s(i)$(ord(i) eq (numGens+2)) .. sum(dynk(i,k),gs(k)) =e= u(i);
PWLEq2s(i)$(ord(i) eq (numGens+2)) .. z(i) =e= sum(dynk(i,k),q(i,k)*gs(k));
PWLEq3s(i)$(ord(i) eq (numGens+2)) .. y(i) =e= sum(dynk(i,k),c(i,k)*gs(k));


u.fx(i)$(ord(i) le numGens) = u_l(i) ;
u.fx(i)$(ord(i) gt numGens) = 1 ;
v.fx(i)$(ord(i) le numGens) = v_l(i) ;
v.fx(i)$(ord(i) gt numGens) = 1 ;

Model UC /all/ ;

UC.optca = 0 ;
UC.optcr = 0 ;

time = timeelapsed;

Solve UC using mip minimizing zc ;

time = timeelapsed - time ;

scalars optca, optcr, modelstatus ;
optca = abs(UC.objest - UC.objval) ;
optcr = optca / max(abs(UC.objest),abs(UC.objval)) ;
modelstatus = UC.Modelstat ;

parameters z_lo(i), z_hi(i) ;
z_lo(i) = z_prev(i) - Rd(i) - (1-u.l(i))*q_min(i) ;
z_hi(i) = z_prev(i) + Ru(i) + v.l(i)*q_min(i) ;

scalars c_bar_cost, h_bar_cost, pw_cost, pw_gen_cost, tot_cost ;
c_bar_cost = sum(i,c_bar(i)*u.l(i));
h_bar_cost = sum(i,h_bar(i)*v.l(i));
pw_cost = sum(i,y.l(i));
pw_gen_cost = sum(i$(ord(i) le numGens),y.l(i));
tot_cost = pw_cost + c_bar_cost + h_bar_cost ;

parameter z_l(i) ;
z_l(i) = z.l(i) ;

execute_unloadIdx 'ubout.gdx', tot_cost, z_l, z_lo, z_hi, optca, time, modelstatus ;
