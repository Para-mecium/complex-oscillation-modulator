function out = MorrisLecar
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = @der4;
out{9} = @der5;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
minf=(1+tanh((kmrgd(1)-par_v1)/par_v2))/2;
winf=(1+tanh((kmrgd(1)-par_v3)/par_v4))/2;
tau=1/cosh((kmrgd(1)-par_v3)/(2*par_v4));
dydt=[(-par_gCa*minf*(kmrgd(1)-par_VCa)-par_gK*kmrgd(2)*(kmrgd(1)-par_VK)-par_gL*(kmrgd(1)-par_VL)+par_Iapp)/par_C;
par_phi*(winf-kmrgd(2))/tau;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(MorrisLecar);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
jac=[ -(par_gL + kmrgd(2)*par_gK + par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) - (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2))/par_C , -(par_gK*(kmrgd(1) - par_VK))/par_C ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) , -par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)) ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
jacp=[ (par_gL*(kmrgd(1) - par_VL) - par_Iapp + par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) + kmrgd(2)*par_gK*(kmrgd(1) - par_VK))/par_C^2 , (kmrgd(2)*par_gK)/par_C , -(kmrgd(2)*(kmrgd(1) - par_VK))/par_C , (par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2))/par_C , -((kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2))/par_C , par_gL/par_C , -(kmrgd(1) - par_VL)/par_C , -(par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2) , -(par_gCa*(kmrgd(1) - par_v1)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2^2) , 0 , 0 , 0 , 1/par_C ; 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) , (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4^2) , cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2) , 0 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
hess1=[ ((par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2 - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2)/par_C , -par_gK/par_C ; (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) + (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^2 , -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) ];
hess2=[ -par_gK/par_C , 0 ; -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
hessp1=[ (par_gL + kmrgd(2)*par_gK + par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) - (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2))/par_C^2 , (par_gK*(kmrgd(1) - par_VK))/par_C^2 ; 0 , 0 ];
hessp2=[ 0 , par_gK/par_C ; 0 , 0 ];
hessp3=[ -kmrgd(2)/par_C , -(kmrgd(1) - par_VK)/par_C ; 0 , 0 ];
hessp4=[ -(par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2) , 0 ; 0 , 0 ];
hessp5=[ -(tanh((kmrgd(1) - par_v1)/par_v2)/2 - ((kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2) + 1/2)/par_C , 0 ; 0 , 0 ];
hessp6=[ 0 , 0 ; 0 , 0 ];
hessp7=[ -1/par_C , 0 ; 0 , 0 ];
hessp8=[ -((par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2) - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2)/par_C , 0 ; 0 , 0 ];
hessp9=[ -((par_gCa*(kmrgd(1) - par_v1)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2^2) + (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2^2) - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_v1)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^3)/par_C , 0 ; 0 , 0 ];
hessp10=[ 0 , 0 ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^2 , (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) ];
hessp11=[ 0 , 0 ; (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^3) + (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^3) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^3 , (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3))/(2*par_v4^2) ];
hessp12=[ 0 , 0 ; (sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) - (cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) , -cosh((kmrgd(1) - par_v3)/(2*par_v4)) ];
hessp13=[ 0 , 0 ; 0 , 0 ];
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
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
tens31=[ ((par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^3 - (3*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2 + (2*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^2*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^3)/par_C , 0 ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(8*par_v4^3) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^3 - (3*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(8*par_v4^3) + (3*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^3) - (2*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^2*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^3 , -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) ];
tens32=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) , 0 ];
tens33=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
tens41=[ ((4*par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^3 + (8*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^2*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^3 - (8*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^4 - (4*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^3*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^4)/par_C , 0 ; (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(16*par_v4^4) - (2*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^4 - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(4*par_v4^4) + (8*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^4 + (4*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^3*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^4 - (4*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^2*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^4 + (3*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^4) , -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(8*par_v4^3) ];
tens42=[ 0 , 0 ; -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(8*par_v4^3) , 0 ];
tens43=[ 0 , 0 ; -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(8*par_v4^3) , 0 ];
tens44=[ 0 , 0 ; 0 , 0 ];
tens45=[ 0 , 0 ; -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(8*par_v4^3) , 0 ];
tens46=[ 0 , 0 ; 0 , 0 ];
tens47=[ 0 , 0 ; 0 , 0 ];
tens48=[ 0 , 0 ; 0 , 0 ];
tens4(:,:,1,1,1) =tens41;
tens4(:,:,1,1,2) =tens42;
tens4(:,:,1,2,1) =tens43;
tens4(:,:,1,2,2) =tens44;
tens4(:,:,2,1,1) =tens45;
tens4(:,:,2,1,2) =tens46;
tens4(:,:,2,2,1) =tens47;
tens4(:,:,2,2,2) =tens48;
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp)
tens51=[ ((8*par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^3)/par_v2^5 - (20*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^3*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^4 - (40*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^4 + (44*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^2*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^5 + (8*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^4*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^5)/par_C , 0 ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(32*par_v4^5) - (8*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^3)/par_v4^5 - (5*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/(2*par_v4^5) - (5*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(32*par_v4^5) + (5*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(4*par_v4^5) - (5*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^2*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^5 - (8*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^4*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^5 + (20*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^5 + (10*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^3*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^5 - (44*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^2*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^5 , -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(16*par_v4^4) ];
tens52=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(16*par_v4^4) , 0 ];
tens53=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(16*par_v4^4) , 0 ];
tens54=[ 0 , 0 ; 0 , 0 ];
tens55=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(16*par_v4^4) , 0 ];
tens56=[ 0 , 0 ; 0 , 0 ];
tens57=[ 0 , 0 ; 0 , 0 ];
tens58=[ 0 , 0 ; 0 , 0 ];
tens59=[ 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(16*par_v4^4) , 0 ];
tens510=[ 0 , 0 ; 0 , 0 ];
tens511=[ 0 , 0 ; 0 , 0 ];
tens512=[ 0 , 0 ; 0 , 0 ];
tens513=[ 0 , 0 ; 0 , 0 ];
tens514=[ 0 , 0 ; 0 , 0 ];
tens515=[ 0 , 0 ; 0 , 0 ];
tens516=[ 0 , 0 ; 0 , 0 ];
tens5(:,:,1,1,1,1) =tens51;
tens5(:,:,1,1,1,2) =tens52;
tens5(:,:,1,1,2,1) =tens53;
tens5(:,:,1,1,2,2) =tens54;
tens5(:,:,1,2,1,1) =tens55;
tens5(:,:,1,2,1,2) =tens56;
tens5(:,:,1,2,2,1) =tens57;
tens5(:,:,1,2,2,2) =tens58;
tens5(:,:,2,1,1,1) =tens59;
tens5(:,:,2,1,1,2) =tens510;
tens5(:,:,2,1,2,1) =tens511;
tens5(:,:,2,1,2,2) =tens512;
tens5(:,:,2,2,1,1) =tens513;
tens5(:,:,2,2,1,2) =tens514;
tens5(:,:,2,2,2,1) =tens515;
tens5(:,:,2,2,2,2) =tens516;
