function [x,resnorm,residual,exitflag,output,lambda,jacobian] = fitter_CARS_resid_interp(x0,lb,ub,PrecondBandWidth_Data,sqrtflag)
% This is an auto generated MATLAB file from Optimization Tool.

% Start with the default options
options = optimset;
% Modify options setting
options = optimset(options,'Display', 'iter-detailed');
options = optimset(options,'Algorithm', 'trust-region-reflective');
options = optimset(options,'MaxIter', 30);
%options = optimset(options,'LevenbergMarquardt', 'off');
options = optimset(options,'PrecondBandWidth', PrecondBandWidth_Data);

if (sqrtflag == 'y')
        [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@CARS_resid_sqrt_interp,x0,lb,ub,options);
else
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@CARS_resid_interp_nrb_multiply,x0,lb,ub,options);
end