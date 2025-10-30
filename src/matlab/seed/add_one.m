function y = add_one(x)
%ADD_ONE  Minimal demo function for CI smoke test.
%   y = ADD_ONE(x) returns x + 1.
%
%   Example:
%     add_one(41)  % -> 42
%
%   This lives in src/ so MATLAB's path handling in CI can use project root.
    y = x + 1;
end
