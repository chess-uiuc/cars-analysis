function [x,resnorm,residual,exitflag,output,lambda,jacobian] = autogencode_3D_interp(x0,lb,ub,PrecondBandWidth_Data)
% This is an auto generated MATLAB file from Optimization Tool.

% Start with the default options
options = optimset;
% Modify options setting
options = optimset(options,'Display', 'iter-detailed');
options = optimset(options,'Algorithm', 'trust-region-reflective');
options = optimset(options,'MaxIter', 30);
%options = optimset(options,'LevenbergMarquardt', 'off');
options = optimset(options,'PrecondBandWidth', PrecondBandWidth_Data);


[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@CARS_resid_3Dlib_interp,x0,lb,ub,options);
