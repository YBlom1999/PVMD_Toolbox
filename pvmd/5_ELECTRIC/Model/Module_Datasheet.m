function [V_,day,night,param] = Module_Datasheet(Acell,J_,T_,Sf,numCells,I,DSValues,InputElectric,ButterFly,JphGenPro)
% Module_Datasheet Calculates the IV curves of each cell of the PV module
% for every moment in time by using datasheet values
%
% This function uses the datasheet values to reconstruct the IV curves for
% every point in time for each cell
%
% Parameters
% ----------
% Acell : double
%   The area of each cell
% J_ : double
%   the absorbed current density in the cell
% T_ : double
%   The temperature of the cell
% Sf : double
%   The shading factor
% numCells : double
%   the number of cells
% I : double
%   The current axis of the IV curves
% DSValues : cell
%   The datasheet values of the cell
% ButterFly : Boolean
%   It indicates if the module has a butterfly connection
% JphGenPro : double
%   The photogenerated current from the GenPro simulation
%
% Returns
% -------
% V_ : double
%   The Voltage axis of the IV curves
% night : logic
%   Indicator of which time steps are night
% param : double
%   The equivalent circuit parameters
%
% Developed by Youri Blom

%% Reading parameters

if ButterFly
    Voc_STC=2*DSValues(1)/numCells;
    Isc_STC=DSValues(2)/2;
    Vmpp_STC=2*DSValues(3)/numCells;
    Impp_STC=DSValues(4)/2;
else
    Voc_STC=DSValues(1)/numCells;
    Isc_STC=DSValues(2);
    Vmpp_STC=DSValues(3)/numCells;
    Impp_STC=DSValues(4);

end
Kp=DSValues(5)*1e-2;
Kv=DSValues(6)*1e-2;
Ki=DSValues(7)*1e-2;

%Constants
k=1.3806e-23;
q=1.6022e-19;


%% Construct temperature dependence
T=(-35:10:95)+273.15;
T_STC = 298.15;

Iph = zeros(1,length(T));
J0 = zeros(1,length(T));
n = zeros(1,length(T));
Rsh = zeros(1,length(T));
Rs = zeros(1,length(T));

Jsc_T = zeros(1,length(T));
Pmpp_T = zeros(1,length(T));
Voc_T = zeros(1,length(T));

n_0 = 1.3;
Rsh_0 = 1e2;
Rs_0 = 1e-5;
for T_i = 1:length(T)
    factor_Pmpp = (1+Kp*(T(T_i)-T_STC));
    factor_Isc = (1+Ki*(T(T_i)-T_STC));
    factor_Voc = (1+Kv*(T(T_i)-T_STC));
    factor_mpp = sqrt(factor_Isc*factor_Voc/factor_Pmpp);

    Jsc_new = Isc_STC/Acell*factor_Isc;
    Jmpp_new = Impp_STC/Acell*factor_Isc/factor_mpp;
    Voc_new = Voc_STC*factor_Voc;
    Vmpp_new = Vmpp_STC*factor_Voc/factor_mpp;

    Vth = k*T(T_i)/q;

    Jph_0 = Jsc_new;
    J0_0 = Jsc_new/(exp(Voc_new/n0/Vth)-1);

    x0 = [Jph_0,J0_0,n_0,Rsh_0,Rs_0];
    eq_fmin = @(x) OptimizeParameters(x,Jsc_new,Jmpp_new,Voc_new,Vmpp_new,Vth);
    options = optimset('Display','none');
    options.TolFun = 1e-8;
    x = fminsearchbnd(eq_fmin,x0,[0,0,0,0,0],[],options);

    Iph(T_i) = x(1);
    J0(T_i) = x(2);
    n(T_i) = x(3);
    Rsh(T_i) = max(x(4),1e1);
    Rs(T_i) = x(5);

    %Calculate Isc, Pmpp, and Voc
    Voltage = 0:0.001:Voc_new*1.2;
    z=(Rs(T_i)*J0(T_i)/(n*Vth*(1+Rs(T_i)/Rsh(T_i))))*exp((Rs(T_i)*(Iph(T_i)+J0(T_i))+Voltage)./(n*Vth*(1+Rs(T_i)/Rsh(T_i))));
    Current=(Iph(T_i)+J0(T_i)-Voltage/(Rsh(T_i)))/(1+Rs(T_i)/Rsh(T_i))-lambertw(z).*(n*Vth)/Rs(T_i);
    Power = Voltage.*Current;
    
    Jsc_T(T_i) = Current(1);
    Pmpp_T(T_i) = max(Power);
    [~,Ind_oc] = min(abs(Current));
    Voc_T(T_i) = Voltage(Ind_oc);

end
R_sh_T = polyfit(T,log(Rsh),4);
R_s_T = polyfit(T,log(Rs),4);
J0_T = polyfit(T,log(J0),4);
n_T = polyfit(T,n,2);
J_sc_T = polyfit(T,Iph,1);

%Calculate the actual temperature dependencies

Ki_actual = (Jsc_T(end)-Jsc_T(1))/(Isc_STC/Acell)/(T(end)-T(1));
Kp_actual = (Pmpp_T(end)-Pmpp_T(1))/(Vmpp_STC*Impp_STC/Acell)/(T(end)-T(1));
Kv_actual = (Voc_T(end)-Voc_T(1))/(Voc_STC)/(T(end)-T(1));

disp(append('actual Ki is ',num2str(100*Ki_actual),'%'))
disp(append('actual Kp is ',num2str(100*Kp_actual),'%'))
disp(append('actual Kv is ',num2str(100*Kv_actual),'%'))

%% Construct irradiance dependence
G=[50,100,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500];
Iabs = JphGenPro*G/1000;
T = 298.15;

Iph = zeros(1,length(G));
J0 = zeros(1,length(G));
n = zeros(1,length(G));
Rsh = zeros(1,length(G));
Rs = zeros(1,length(G));

for G_i = 1:length(G)
    Jsc_new = Isc_STC/Acell*G(G_i)/1000;
    Jmpp_new = Impp_STC/Acell*G(G_i)/1000;
    Voc_new = Voc_STC+n_0*k*T/q*log(G(G_i)/1000);
    Vmpp_new = Vmpp_STC+n*k*T/q*log(G(G_i)/1000)*Vmpp_STC/Voc_STC;

    Vth = k*T/q;

    Jph_0 = Jsc_new;
    J0_0 = Jsc_new/(exp(Voc_new/n0/Vth)-1);
    

    x0 = [Jph_0,J0_0,n_0,Rsh_0,Rs_0];
    eq_fmin = @(x) OptimizeParameters(x,Jsc_new,Jmpp_new,Voc_new,Vmpp_new,Vth);
    options = optimset('Display','none');
    options.TolFun = 1e-8;
    x = fminsearchbnd(eq_fmin,x0,[0,0,0,0,0],[],options);

    Iph(G_i) = x(1);
    J0(G_i) = x(2);
    n(G_i) = x(3);
    Rsh(G_i) = max(x(4),1e1);
    Rs(G_i) = x(5);

end
J_sc_J = polyfit(Iabs,Iph,2);
R_sh_J = polyfit(Iabs,log(Rsh),4);
R_s_J = polyfit(Iabs,log(Rs),4);
J0_J = polyfit(Iabs,log(J0),4);
n_J = polyfit(Iabs,n,2);

%% Running current mapping

%Reduce current based on shading
J_=J_*(1-Sf*0.01);

%Sum over all cells
j=sum(J_,2);
%timesteps with zero current (night)
night=find(j==0);
%timesteps with current greater than zero (day)
day=find(j>0);
%delete timesteps (rows) with zero current
J_(night,:)=[];
T_(night,:)=[];

Tmin=min(min(T_));

%Turn matrix into a vector
J_=reshape(J_',numel(J_),1);
T_=reshape(T_',numel(T_),1);
%Mapp cells into 0.4 A/m2 and 0.3 degC bins
J_bin=0.4;%0.4;
T_bin=0.3;
pointer=zeros(length(J_),2);
pointer(:,1)=ceil((J_-J_bin)/J_bin);
pointer(:,2)=ceil((T_-Tmin)/T_bin);
%Identify unique conditions bins
cond=unique(pointer,'rows');
%Removes conditions with J_ < J_bin A/m2
index=find(cond(:,1)>1,1);
cond=cond(index:end,:);



%% Caculating IV curves

%IV curves for each condition with resolution of I
sim=zeros(length(cond(:,1)),length(I));
%Voltage steps
V=[-15:-1 0:5e-3:1.5];

export_param=zeros(length(cond(:,1)),5); %added by Youri

%Calculate IV-curves for all bins
for i=1:length(cond(:,1))
    J=cond(i,1)*J_bin-0.5*J_bin;
    T=cond(i,2)*T_bin-0.5*T_bin+Tmin;
    %prepareing parameters for one diode model solved by lambert-wfunction
    Iph=Acell*polyval(J_sc_T,T).*polyval(J_sc_J,J)./polyval(J_sc_T,298.15);
    Vth=k*T/q;
    %Correction for cell area
    %due to fitting JV (instead of IV) from ASA
    %Rs & Rsh would be m2*Ohm and I0 would be in A/m2 without corrrection
    Rsh=max((1/Acell)*exp(polyval(R_sh_T,T)).*exp(polyval(R_sh_J,J))./exp(polyval(R_sh_T,298.15)),0.001);
    Rs=(1/Acell)*exp(polyval(R_s_T,T)).*exp(polyval(R_s_J,J))./exp(polyval(R_s_T,298.15));
    I0=Acell*exp(polyval(J0_T,T)).*exp(polyval(J0_J,J))./exp(polyval(J0_T,298.15)) ;
    n=polyval(n_T,T).*polyval(n_J,J)./polyval(n_T,298.15);
    %one diode model solved by lambert-wfunction
    z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
    I_=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;

    % Apply different model for negative currents
    if ~isnan(InputElectric.reverseParam)
        Be = InputElectric.reverseParam(1);
        phi_t = InputElectric.reverseParam(2);
        V_b = InputElectric.reverseParam(3);
        c = InputElectric.reverseParam(4);
        ind_neg = V < 0;
        I_n = Iph-V(ind_neg)/Rsh+c*V(ind_neg).^2;
        K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V(ind_neg)))));
        I_(ind_neg) = I_n./(1-K_e);
    end
    
    %Interpolation to comon currents
    ind = isinf(I_).*I_<0;
    I_(ind) = -1e5*V(ind);

    sim(i,:)=interp1(I_,V,I,'linear','extrap');

    export_param(i,:)=[Iph,Rs,Rsh,n,I0]; %added by Youri
end

%% Remap the bins in the time slots
V_=zeros(length(day)*numCells,length(I));
param = zeros(length(day)*numCells,5); %added by Youri
for i=1:length(cond(:,1))
    index=find(pointer(:,1)==cond(i,1) & pointer(:,2)==cond(i,2));
    V_(index,:)=repmat(sim(i,:),length(index),1);
    param(index,:) = repmat(export_param(i,:),length(index),1); %added by Youri
end


end

function [error] = OptimizeParameters(x,Isc,Impp,Voc,Vmpp,Vth)
% OptimzeParameters predicts the error of an IV curve with respect to the external parameters.
% x contains the parameters of one-diode model, that can be used to make
% the IV curve
%
% Parameters
% ----------
% x : array
%   The values of the one-diode model
% n: double
%   The specified ideality factor
% Isc: double
%   The desired short-circuit current
% Impp: double
%   The desired maximum power point current
% Voc: double
%   The desired open circuit voltage
% Vmpp: double
%   The desired maximum power point voltage
% Vth: double
%   The thermal voltage
%
% Returns
% --------
% error: double
%   The error that quantifies the quality of the fit.

Jph = x(1);
J0 = x(2);
n = x(3);
Rsh = max(x(4),1e1);
Rs = x(5);

%Error 1 is the difference between the current source and the diode + shunt
%current at short circuit + output current
err1 = (Jph-J0*(exp(Isc*Rs/(n*Vth))-1)-Isc*Rs/Rsh-Isc)^2;

%Error 2 is the difference between the current source and the diode + shunt
%current at open circuit condition
err2 = (Jph-J0*(exp(Voc/(n*Vth))-1)-Voc/Rsh)^2;

%Error 1 is the difference between the current source and the diode + shunt
%current at maximum power point + output current
err3 = (Jph-J0*(exp((Impp*Rs+Vmpp)/(n*Vth))-1)-(Impp*Rs+Vmpp)/Rsh-Impp)^2;

%Error 4 is the derivative of the power (this should be 0 at MPP).
exp_term = (Impp*Rs+Vmpp)/(n*Vth);
der_cur = -1*(1/Rsh + J0/(Vth*n)*exp(exp_term))/(1+Rs/Rsh+J0*Rs/(Vth*n)*exp(exp_term));
err4 = (der_cur*Vmpp+Impp)^2;



error = err1+err2+err3+err4;
end

function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
  x0u(k:n) = [];
end

% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end

% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end

% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);

% undo the variable transformations into the original space
x = xtransform(xu,params);

% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);

% Use a nested function as the OutputFcn wrapper
  function stop = outfun_wrapper(x,varargin);
    % we need to transform x first
    xtrans = xtransform(x,params);
    
    % then call the user supplied OutputFcn
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    
  end

end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function

% transform
xtrans = xtransform(x,params);

% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});

end % sub function intrafun end

% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end