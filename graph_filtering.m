function coeff = graph_filtering(L, l_max, kernel, signals, cheb_ord)
%GRAPH_FILTERING Perform graph filtering with Chebyshev polynomials.
%Performs graph filtering of a set of signals with a specified spectral
%kernel. The spectral kernel is approximated with Chebyshev polynomials,
%which are then applied directly on the graph Laplacian. This avoids the
%costly eigendecomposition of the graph Laplacian matrix.
%
% Syntax:  graph_filtering(f_fmri,f_mask,G,tau)
%
% Inputs:
%    L - Laplacian matrix of graph on which data reside.
%    l_max - Highest eigenvalue of graph Laplacian.
%    kernel - Function handle to spectral kernel of filter. Kernel function
%      has support in the [0,2] range.
%    signals - Array containing signals defined on the vertices of the
%      graph. Signals are arranged in the columns of the array.
%    cheb_ord - Order of Chebyshev polynomials used to approximate the
%      spectral kernel. Higher values provide a better approximation of the
%      filter kernel.
%
% Examples:
%    coeff = graph_filtering(G.L, G.l_max, @(x) exp(-x*tau), ...
%      signals, cheb_ord)
%
% See also: DSS_FILTER_FMRI

% Author: David Abramian, adapted from Martin Larsson, 2017
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% check input
assert(size(L,2) == size(L,1), 'A is not square.');
assert(size(signals,1) == size(L,1),...
    'Signal size does not match number of vertices');

arange = [0 l_max];

% approximate kernel with Chebyshev polynomial
c = sgwt_cheby_coeff(kernel, cheb_ord, cheb_ord+1, arange);

% decompose signal and return low-pass component
coeff = sgwt_cheby_op(signals, L, c, arange);

