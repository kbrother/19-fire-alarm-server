function S = tt_matrix2cellstr(M)
%TT_MATRIX2CELLSTR Convert a matrix to a cell array of strings.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.


fmt = get(0,'FormatSpacing');
format compact
S = evalc('disp(M)');
if isempty(S)
    S = {''};
    return;
end
set(0,'FormatSpacing',fmt)
S = textscan(S,'%s','delimiter','\n','whitespace','');
S = S{1};
end
