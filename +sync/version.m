function ver = version()
% v1 2019/05/15
% -Two aligner classes now exist under the "BaseAligner" superclass
%   -Original "Aligner" class renamed to "PatternAligner"
%   -New "SimpleAligner" class handles pre-matched pulse trains
% -Remove "align()" method; alignment is now 
%
% v2 2019/09/02
% -Fix bug where sync.Network didn't cache aligners correctly

ver = 2;
end