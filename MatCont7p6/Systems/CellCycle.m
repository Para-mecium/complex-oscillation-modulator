function out = CellCycle
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
dydt=[par_k1-(par_kp2+par_kpp2*kmrgd(2))*kmrgd(1);
(par_kp3+par_kpp3*kmrgd(3))*(1-kmrgd(2))/(par_J3+1-kmrgd(2))-par_k4*par_m*kmrgd(2)*kmrgd(1)/(par_J4+kmrgd(2));
par_kp5+par_kpp5*(par_m*kmrgd(1))^par_n/(par_J5^par_n+(par_m*kmrgd(1))^par_n)-par_k6*kmrgd(3);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(CellCycle);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
jac=[ - par_kp2 - kmrgd(2)*par_kpp2 , -kmrgd(1)*par_kpp2 , 0 ; -(kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4) , (kmrgd(1)*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - ((par_kp3 + kmrgd(3)*par_kpp3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 - (kmrgd(1)*par_m*par_k4)/(kmrgd(2) + par_J4) - (par_kp3 + kmrgd(3)*par_kpp3)/(par_J3 - kmrgd(2) + 1) , -(par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1) ; (par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , -par_k6 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
jacp=[ 1 , -kmrgd(1) , -kmrgd(1)*kmrgd(2) , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -(kmrgd(2) - 1)/(par_J3 - kmrgd(2) + 1) , -(kmrgd(3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1) , -(kmrgd(1)*kmrgd(2)*par_m)/(kmrgd(2) + par_J4) , ((par_kp3 + kmrgd(3)*par_kpp3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 , (kmrgd(1)*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 , 0 , 0 , 0 , 0 , 0 , -(kmrgd(1)*kmrgd(2)*par_k4)/(kmrgd(2) + par_J4) ; 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , (kmrgd(1)*par_m)^par_n/((kmrgd(1)*par_m)^par_n + par_J5^par_n) , -kmrgd(3) , -(par_J5^(par_n - 1)*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n)/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , (par_kpp5*log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^par_n)/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_kpp5*(kmrgd(1)*par_m)^par_n*(log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^par_n + par_J5^par_n*log(par_J5)))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , (kmrgd(1)*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (kmrgd(1)*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
hess1=[ 0 , -par_kpp2 , 0 ; 0 , (kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (par_m*par_k4)/(kmrgd(2) + par_J4) , 0 ; (2*par_m^2*par_n^2*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 - (2*par_m^2*par_n^2*par_kpp5*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 + (par_m^2*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_m^2*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , 0 ];
hess2=[ -par_kpp2 , 0 , 0 ; (kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (par_m*par_k4)/(kmrgd(2) + par_J4) , (2*kmrgd(1)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (2*(par_kp3 + kmrgd(3)*par_kpp3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^3 - (2*(par_kp3 + kmrgd(3)*par_kpp3))/(par_J3 - kmrgd(2) + 1)^2 - (2*kmrgd(1)*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 , - par_kpp3/(par_J3 - kmrgd(2) + 1) - (par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 ; 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 ; 0 , - par_kpp3/(par_J3 - kmrgd(2) + 1) - (par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 , 0 ; 0 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
hessp1=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp2=[ -1 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp3=[ -kmrgd(2) , -kmrgd(1) , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp4=[ 0 , 0 , 0 ; 0 , - (kmrgd(2) - 1)/(par_J3 - kmrgd(2) + 1)^2 - 1/(par_J3 - kmrgd(2) + 1) , 0 ; 0 , 0 , 0 ];
hessp5=[ 0 , 0 , 0 ; 0 , - kmrgd(3)/(par_J3 - kmrgd(2) + 1) - (kmrgd(3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 , -(kmrgd(2) - 1)/(par_J3 - kmrgd(2) + 1) ; 0 , 0 , 0 ];
hessp6=[ 0 , 0 , 0 ; -(kmrgd(2)*par_m)/(kmrgd(2) + par_J4) , (kmrgd(1)*kmrgd(2)*par_m)/(kmrgd(2) + par_J4)^2 - (kmrgd(1)*par_m)/(kmrgd(2) + par_J4) , 0 ; 0 , 0 , 0 ];
hessp7=[ 0 , 0 , 0 ; 0 , (par_kp3 + kmrgd(3)*par_kpp3)/(par_J3 - kmrgd(2) + 1)^2 + (2*(par_kp3 + kmrgd(3)*par_kpp3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^3 , (par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^2 ; 0 , 0 , 0 ];
hessp8=[ 0 , 0 , 0 ; (kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 , (kmrgd(1)*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (2*kmrgd(1)*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 , 0 ; 0 , 0 , 0 ];
hessp9=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp10=[ 0 , 0 , 0 ; 0 , 0 , 0 ; (par_m*par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_m*par_n*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , 0 ];
hessp11=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , -1 ];
hessp12=[ 0 , 0 , 0 ; 0 , 0 , 0 ; (2*par_J5^(par_n - 1)*par_m*par_n^2*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 - (par_J5^(par_n - 1)*par_m*par_n^2*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , 0 ];
hessp13=[ 0 , 0 , 0 ; 0 , 0 , 0 ; (par_m*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_m*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 + (par_m*par_n*par_kpp5*log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1)*(log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^par_n + par_J5^par_n*log(par_J5)))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 - (2*par_m*par_n*par_kpp5*log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 + (2*par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1)*(log(kmrgd(1)*par_m)*(kmrgd(1)*par_m)^par_n + par_J5^par_n*log(par_J5)))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 , 0 , 0 ];
hessp14=[ 0 , 0 , 0 ; -(kmrgd(2)*par_k4)/(kmrgd(2) + par_J4) , (kmrgd(1)*kmrgd(2)*par_k4)/(kmrgd(2) + par_J4)^2 - (kmrgd(1)*par_k4)/(kmrgd(2) + par_J4) , 0 ; (par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 - (2*kmrgd(1)*par_m*par_n^2*par_kpp5*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 + (2*kmrgd(1)*par_m*par_n^2*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 + (kmrgd(1)*par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) - (kmrgd(1)*par_m*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
hessp(:,:,10) =hessp10;
hessp(:,:,11) =hessp11;
hessp(:,:,12) =hessp12;
hessp(:,:,13) =hessp13;
hessp(:,:,14) =hessp14;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
tens31=[ 0 , 0 , 0 ; 0 , 0 , 0 ; (6*par_m^3*par_n^3*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1)*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 - (2*par_m^3*par_n^2*par_kpp5*(2*par_n - 2)*(kmrgd(1)*par_m)^(2*par_n - 3))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 - (2*par_m^3*par_n^2*par_kpp5*(kmrgd(1)*par_m)^(par_n - 1)*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 + (par_m^3*par_n*par_kpp5*(kmrgd(1)*par_m)^(par_n - 3)*(par_n - 1)*(par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n) + (2*par_m^3*par_n^2*par_kpp5*(2*par_n - 2)*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(2*par_n - 3))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 - (6*par_m^3*par_n^3*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1)*(kmrgd(1)*par_m)^(2*par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^4 + (2*par_m^3*par_n^2*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 1)*(kmrgd(1)*par_m)^(par_n - 2)*(par_n - 1))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^3 - (par_m^3*par_n*par_kpp5*(kmrgd(1)*par_m)^par_n*(kmrgd(1)*par_m)^(par_n - 3)*(par_n - 1)*(par_n - 2))/((kmrgd(1)*par_m)^par_n + par_J5^par_n)^2 , 0 , 0 ];
tens32=[ 0 , 0 , 0 ; 0 , (2*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (2*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 ; 0 , (2*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (2*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 , 0 ; 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 ; (2*par_m*par_k4)/(kmrgd(2) + par_J4)^2 - (2*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 , (6*kmrgd(1)*kmrgd(2)*par_m*par_k4)/(kmrgd(2) + par_J4)^4 - (6*(par_kp3 + kmrgd(3)*par_kpp3)*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^4 - (6*kmrgd(1)*par_m*par_k4)/(kmrgd(2) + par_J4)^3 - (6*(par_kp3 + kmrgd(3)*par_kpp3))/(par_J3 - kmrgd(2) + 1)^3 , - (2*par_kpp3)/(par_J3 - kmrgd(2) + 1)^2 - (2*par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^3 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , - (2*par_kpp3)/(par_J3 - kmrgd(2) + 1)^2 - (2*par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^3 , 0 ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , - (2*par_kpp3)/(par_J3 - kmrgd(2) + 1)^2 - (2*par_kpp3*(kmrgd(2) - 1))/(par_J3 - kmrgd(2) + 1)^3 , 0 ; 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_k1,par_kp2,par_kpp2,par_kp3,par_kpp3,par_k4,par_J3,par_J4,par_kp5,par_kpp5,par_k6,par_J5,par_n,par_m)
