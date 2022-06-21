function progresss(i, i_max, string_in)
%PROGRESSS - Simple iteration counter utility.
% Prints updating iteration counter on a single line.
%
% Syntax:  progresss(i,i_max,string_in)
%
% Inputs:
%    i - Current iteration index.
%    i_max - Maximum iteration index.
%    string_in - String to precede counter.
%
% Example:
%    progresss(1, 10, 'Iterating... ')
%      Prints sting_in followed by counter.
%    progresss(5, 10, 'Iterating... ')
%      Erases counter and prints new counter.
%    progresss(10, 10, 'Itreating... ')
%      Prints final counter and newline.
%    progresss(11, 10, 'Iterating... ') 
%      Erases newline, prints 'Done' and newline.

% Author: David Abramian
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2020; Last revision: 13-May-2021


digits_str = num2str(i_max);
n_digits = length(digits_str);

string = ['%', num2str(n_digits), 'i/%', num2str(n_digits), 'i'];

if i == 1  % first iter
    fprintf([string_in, string], i, i_max)
    return
elseif i < i_max  % standard iter
    fprintf([repmat('\b',1,2*n_digits+1), string], i, i_max)
    return
elseif i == i_max  % last iter
    fprintf([repmat('\b',1,2*n_digits+1), string, '\n'], i, i_max)
    return
elseif i == i_max + 1  % optional last iter
    fprintf('\b Done \n')
end

