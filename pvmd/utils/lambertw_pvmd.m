function w = lambertw_pvmd(branch, x)
% Lambert_W Lambert's W function.
%    W = lambertw(Z) solves W*exp(W) = Z.
%    W = lambertw(K,Z) is the K-th branch of this multi-valued function.
 
%    References:
%    [1] Robert M. Corless, G. H. Gonnet, D. E. G. Hare,
%    D. J. Jeffrey, and D. E. Knuth, "On the Lambert W Function",
%    Advances in Computational Mathematics, volume 5, 1996, pp. 329-359.
 
%    [2] Corless, Jeffrey, Knuth, "A Sequence of Series
%    for The Lambert W Function", ISSAC '97
%    Copyright Lateef Adewale Kareem 2022.
if nargin < 2
    x = branch;  branch = 0;
end
% Effective starting guess
v = inf*ones(size(x));
if numel(branch) == 1
    if branch == 0
       w = ones(size(x));  % Start above -1
    else  
       w = -2*ones(size(x));  % Start below -1
    end
    if(branch == 0 || branch  == -1)
        w = Haley(w, v, x, 0);
    else
        w = Haley(w, v, x, branch);
    end
else
    w = ones(size(x));
    w(branch ~= 0) = -2;
    indx = branch == 0 | branch == -1;
    w(indx) = Haley(w(indx), v(indx), x(indx), 0);
    indx = ~(branch == 0 & branch ==-1);
    w(indx) = Haley(w(indx), v(indx), x(indx), branch(indx));
end
w(x==0)=0;
w(x==-1/exp(1))=-1;
end
function w = Haley(w, v, x, branch)
    % Haley's method

%     fig = figure;
%     hold on; box on;
%     w_range = 0.001:0.001:1;
%     f_range = w_range + log(w_range)-log(x);
%     plot(w_range,f_range);

    %indices below 1e-2 take the output W = x
    ind = x>10^-1;

    w_sel = w(ind);
    x_sel = x(ind);
    v_sel = v(ind);

    while any(abs(w_sel - v_sel)./abs(w_sel) > 1.e-15)
       v_sel = w_sel;
       f = w_sel + log(w_sel) - log(x_sel) - 2*pi*1i*branch;
       fp = 1 + 1./w_sel;
       fpp = -1./(w_sel.*w_sel);
       % Iterate to make this quantity zero
       w_sel = w_sel - f ./ (fp - f .* fpp ./ (2 * fp));

       %if w is negative log(w) becomes complex, therefore a different
       %iteration method is selected
       ind_neg = w_sel < 0;
       w_sel(ind_neg) = v_sel(ind_neg)-0.1*f(ind_neg);
    end
    w(ind) = w_sel;
    w(~ind) = 0;
    for n = 1:10
        w(~ind) = w(~ind) + (-n)^(n-1)/factorial(n)*x(~ind).^n;
    end
end
