function [dp,lambda] = SDT(Hit,FA,Miss,CR)
% [dp,lambda] = SDT(Hit,FA,Miss,CR)
%   Gaussian Equal Variance Signal Detection Theory
% 
% Inputs: 
% Hit    the number of hits
% FA     the number of false alarms
% Miss   the number of misses
% CR     the number of correct rejections
%
% Outputs:
% dp     dprime
% lambda sensory criterion
%

nsignal = Hit+Miss;
nnoise  = FA+CR;

pHit    = Hit/nsignal;
pMiss   = Miss/nsignal;
pFA     = FA/nnoise;
pCR     = CR/nnoise;

lambda  = icdf( 'norm', pCR,   0, 1 );
l_star  = icdf( 'norm', pMiss, 0, 1 );
dp      = lambda - l_star;
end