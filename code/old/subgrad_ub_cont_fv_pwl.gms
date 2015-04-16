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
         y(i)            piecewise linear cost
         yy(i)           piecewise linear cost for future value function
         z(i)            real power from generator i
         g(i,kk)         variables for piecewise linear cost
         gb(k)           pwl variable for buy generator
         gs(k)           pwl variable for sell generator
*         gg(i,kk)        variables for piecewise linear future value function
         zc              objective
         ;

*Positive Variable y, z, g, gg ;
Positive Variable y, z, g, gb, gs ;
Binary Variable u ;

$GDXin %probdata%
Set      dynk(i,k) ;
*Set dynkb(k) ;
*Set dynks(k) ;

Parameter
         numPts(i)     number of valid cost pts for each generator ;
$load numPts

dynk(i,k) = YES$(ord(k) le numPts(i));
*dynkb(k) = YES$(ord(k) le numPts(i$(ord(i) eq (numGens+1))));
*dynks(k) = YES$(ord(k) le numPts(i$(ord(i) eq (numGens+2))));
*display dynkb, dynks ;
*$exit

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

$GDXin

$GDXIN ubdata.gdx
Parameter
         z_lo(i)    lowest possible power level for each generator ;
$load z_lo

Parameter
         z_hi(i)    highest possible power level for each generator ;
$load z_hi

Parameter
         u_i(i)     on off generator decisions  ;
$load u_i

Parameter
         v_i(i)     turn on generator decisions  ;
$load v_i

Scalar
         d          observed demand that needs to be satisfied ;
$load d

*Parameter
*         fv_i(i)    generators which require future value function optimized  ;
*$load fv_i

*Scalar
*         sumV       amount of non-z related future value function ;
*$load sumV

*Parameter
*         Plevels(i,kk)   discretized power levels ;
*$load Plevels

*Parameter
*         fV(i,kk)   discretized points on the value function ;
*$onUNDF
*$load fV

Parameter
         Voff(i)         future cost incurred when generator off
$onUNDF
$load Voff

Parameter
         Von(i,kk)       future cost incurred when generator on
$onUNDF
$load Von

Parameter
         states(i)  current state of each generator ;
$load states

Parameter
         z_idx_lo(i)   valid lo idx for power level for future value function ;
$load z_idx_lo

Parameter
         z_idx_hi(i)   valid hi idx for power level for future value function ;
$load z_idx_hi

Set      dynikk(i,kk) ;

dynikk(i,kk) = YES$((ord(kk) ge z_idx_lo(i)) and (ord(kk) le z_idx_hi(i)));

$GDXIN

*display dynikk, Von, Voff, q ;

Scalar time keep track of total time ;

Equations
         cost                  objective function
         demEq                 demand satisfaction constraint
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
*         PWLEq11(i)            1st set of PWL eq for future value function
*         PWLEq22(i)            2nd set of PWL eq for future value function
*         PWLEq33(i)            3rd set of PWL eq for future value function
         ;

*cost ..              zc=e=sum(i,y(i)+c_bar(i)*u(i)+h_bar(i)*v_i(i)) ;
cost ..              zc=e=sum(i,y(i)+c_bar(i)*u(i)+h_bar(i)*v_i(i))+sum(i$(ord(i) le numGens),yy(i)+(1-u(i))*Voff(i)) ;
demEq ..             sum(i$(ord(i) ne %subprob%),z(i))-z('%subprob%')=e=d;
*PWLEq1(i) ..         sum(dynk(i,k),g(i,k)) =e= u(i);
*PWLEq2(i) ..         z(i) =e= sum(dynk(i,k),q(i,k)*g(i,k));
*PWLEq3(i) ..         y(i) =e= sum(dynk(i,k),c(i,k)*g(i,k));
*PWLEq11(i)$(ord(i) le numGens) .. sum(dynikk(i,kk),gg(i,kk)) =e= u(i) ;
*PWLEq22(i)$(ord(i) le numGens) .. z(i) =e= sum(dynikk(i,kk),Plevels(i,kk)*gg(i,kk));
*PWLEq33(i)$(ord(i) le numGens) .. yy(i) =e= sum(dynikk(i,kk),Von(i,kk)*gg(i,kk));

PWLEq1(i)$(ord(i) le numGens) .. sum(dynikk(i,kk),g(i,kk)) =e= u(i);
PWLEq2(i)$(ord(i) le numGens) .. z(i) =e= sum(dynikk(i,kk),Pow(kk,i)*g(i,kk));
PWLEq3(i)$(ord(i) le numGens) .. y(i) =e= sum(dynikk(i,kk),Feval(kk,i)*g(i,kk));
PWLEq4(i)$(ord(i) le numGens) .. yy(i) =e= sum(dynikk(i,kk),Von(i,kk)*g(i,kk));

* for buy generator
PWLEq1b(i)$(ord(i) eq (numGens+1)) .. sum(dynk(i,k),gb(k)) =e= u(i);
PWLEq2b(i)$(ord(i) eq (numGens+1)) .. z(i) =e= sum(dynk(i,k),q(i,k)*gb(k));
PWLEq3b(i)$(ord(i) eq (numGens+1)) .. y(i) =e= sum(dynk(i,k),c(i,k)*gb(k));

* for sell generator
PWLEq1s(i)$(ord(i) eq (numGens+2)) .. sum(dynk(i,k),gs(k)) =e= u(i);
PWLEq2s(i)$(ord(i) eq (numGens+2)) .. z(i) =e= sum(dynk(i,k),q(i,k)*gs(k));
PWLEq3s(i)$(ord(i) eq (numGens+2)) .. y(i) =e= sum(dynk(i,k),c(i,k)*gs(k));


* force buy / sell generators to not have a cost
*yy.fx(i)$(ord(i) gt numGens) = 0  ;
z.lo(i) = z_lo(i) ;
z.up(i) = z_hi(i) ;

* fix on / off decisions to be the ones from subgradient method
u.fx(i) = u_i(i) ;

Model UC /all/ ;

*UC.optcr = 0.005 ;
UC.optcr = 0 ;

time = timeelapsed;

Solve UC using mip minimizing zc ;

time = timeelapsed - time ;

scalar c_bar_cost, h_bar_cost, pw_cost, pw_gen_cost, tot_cost, fval_cost ;
*pw_buy_cost, pw_sell_cost ;

c_bar_cost = sum(i,c_bar(i)*u.l(i));
h_bar_cost = sum(i,h_bar(i)*v_i(i));
pw_cost = sum(i,y.l(i));
pw_gen_cost = sum(i$(ord(i) le numGens),y.l(i));
tot_cost = pw_cost + c_bar_cost + h_bar_cost ;
*pw_sell_cost = y.l('%subprob%');

parameter z_l(i) ;
z_l(i) = z.l(i) ;

*display z_l, u.l ;
*display z_l, u.l, yy.l, gg.l ;

*scalar tcomp, texec, telapsed ;
*
*tcomp = TimeComp ;
*texec = TimeExec ;
*telapsed = TimeElapsed ;
*
*display tcomp, texec, telapsed ;

*display z.l, zc.l, y.l, g.l, c_bar_cost, h_bar_cost, pw_cost, pw_gen_cost, tot_cost ;
*display z.l, zc.l, y.l, g.l, yy.l, gg.l, c_bar_cost, h_bar_cost, pw_cost, pw_gen_cost, tot_cost ;
execute_unloadIdx 'ubout.gdx', tot_cost, z_l, time ;
