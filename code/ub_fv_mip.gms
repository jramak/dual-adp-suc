Sets
         i       generators                              / 1 * %subprob% /
         k       max num piecewise pts                   / 1 * 11 /
         kk      max piecewise pts for fut val func      / 1 * %nPts% /
         j       for each demand instance                / 1 * 10 /
         ;

Scalar
         numGens       number of generators without buy and sell generators ;

numGens = %subprob% - 2 ;

Variables
         u(i)            generator on off decisions
         v(i)            turn on variable
         y(i,j)          piecewise linear cost
         yy(i,j)         piecewise linear cost for future value function
         z(i,j)          real power from generator i
         g(i,kk,j)       variables for piecewise linear cost
         gb(k,j)         pwl variable for buy generator
         gs(k,j)         pwl variable for sell generator
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
         p_s(j)        probability of the demand instances ;
$load p_s

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
         dj(j)      observed demand that needs to be satisfied ;
$load dj

Parameter
         Voff(i)         future cost incurred when generator off  ;
$onUNDF
$load Voff

Parameter
         Von(i,kk)       future cost incurred when generator on ;
$onUNDF
$load Von

Parameter u_fx_lo(i) ;
$load u_fx_lo

Parameter u_fx_hi(i) ;
$load u_fx_hi

$GDXIN

Scalar time keep track of total time ;

Equations
         cost                    objective function
         demEq(j)                demand satisfaction constraint
         turnOn(i)               equation to set turn on variable
         rampUpEq(i,j)           ramp up constraint
         rampDownEq(i,j)         ramp down constraint

         PWLEq1(i,j)             1st set of PWL eq
         PWLEq2(i,j)             2nd set of PWL eq
         PWLEq3(i,j)             3rd set of PWL eq
         PWLEq4(i,j)             4th set of PWL eq

         PWLEq1b(i,j)            1st set of PWL eq for buy gen
         PWLEq2b(i,j)            2nd set of PWL eq for buy gen
         PWLEq3b(i,j)            3rd set of PWL eq for buy gen

         PWLEq1s(i,j)            1st set of PWL eq for sell gen
         PWLEq2s(i,j)            2nd set of PWL eq for sell gen
         PWLEq3s(i,j)            3rd set of PWL eq for sell gen
         ;

cost ..              zc=e=sum(j,p_s(j)*(sum(i,y(i,j)+c_bar(i)*u(i)+h_bar(i)*v(i))+sum(i$(ord(i) le numGens),yy(i,j)+(1-u(i))*Voff(i)))) ;
demEq(j) ..          sum(i$(ord(i) ne %subprob%),z(i,j))-z('%subprob%',j)=e=dj(j);
turnOn(i) ..         v(i)=g=u(i)-u_prev(i) ;
rampUpEq(i,j)$(ord(i) le numGens) .. z(i,j) =l= z_prev(i) + Ru(i) + v(i)*q_min(i);
rampDownEq(i,j)$(ord(i) le numGens) .. z_prev(i) - Rd(i) - (1-u(i))*q_min(i) =l= z(i,j);

PWLEq1(i,j)$(ord(i) le numGens) .. sum(kk,g(i,kk,j)) =e= u(i);
PWLEq2(i,j)$(ord(i) le numGens) .. z(i,j) =e= sum(kk,Pow(kk,i)*g(i,kk,j));
PWLEq3(i,j)$(ord(i) le numGens) .. y(i,j) =e= sum(kk,Feval(kk,i)*g(i,kk,j));
PWLEq4(i,j)$(ord(i) le numGens) .. yy(i,j) =e= sum(kk,Von(i,kk)*g(i,kk,j));

* for buy generator
PWLEq1b(i,j)$(ord(i) eq (numGens+1)) .. sum(dynk(i,k),gb(k,j)) =e= u(i);
PWLEq2b(i,j)$(ord(i) eq (numGens+1)) .. z(i,j) =e= sum(dynk(i,k),q(i,k)*gb(k,j));
PWLEq3b(i,j)$(ord(i) eq (numGens+1)) .. y(i,j) =e= sum(dynk(i,k),c(i,k)*gb(k,j));

* for sell generator
PWLEq1s(i,j)$(ord(i) eq (numGens+2)) .. sum(dynk(i,k),gs(k,j)) =e= u(i);
PWLEq2s(i,j)$(ord(i) eq (numGens+2)) .. z(i,j) =e= sum(dynk(i,k),q(i,k)*gs(k,j));
PWLEq3s(i,j)$(ord(i) eq (numGens+2)) .. y(i,j) =e= sum(dynk(i,k),c(i,k)*gs(k,j));


* force buy and sell generators to be turned on
u.lo(i)$(ord(i) le numGens) = u_fx_lo(i) ;
u.up(i)$(ord(i) le numGens) = u_fx_hi(i) ;
u.fx(i)$(ord(i) gt numGens) = 1 ;
v.fx(i)$(ord(i) gt numGens) = 1 ;

Model UC /all/ ;

*UC.optca = 0 ;
*UC.optcr = 0 ;
* 0.5% optimality gap
UC.optcr = 0.005 ;
UC.threads = 4 ;

time = timeelapsed ;

Solve UC using mip minimizing zc ;

time = timeelapsed - time ;

scalars optca, optcr, modelstatus ;
optca = abs(UC.objest - UC.objval) ;
optcr = optca / max(abs(UC.objest),abs(UC.objval)) ;
modelstatus = UC.Modelstat ;

scalars c_bar_cost, h_bar_cost ;

parameter pw_cost(j), pw_gen_cost(j), tot_cost(j) ;
c_bar_cost = sum(i,c_bar(i)*u.l(i));
h_bar_cost = sum(i,h_bar(i)*v.l(i));
pw_cost(j) = sum(i,y.l(i,j));
pw_gen_cost(j) = sum(i$(ord(i) le numGens),y.l(i,j));
tot_cost(j) = pw_cost(j) + c_bar_cost + h_bar_cost ;

parameter u_l(i), z_l(i,j) ;
u_l(i) = u.l(i) ;
z_l(i,j) = z.l(i,j) ;

execute_unloadIdx 'ubout.gdx', tot_cost, u_l, z_l, optca, optcr, time, modelstatus ;
