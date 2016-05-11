function [CF] = fitness_function_ULA_GA(excite,er)
%% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: Beam Pattern Synthesis in Sensor Arrays Using Optimisation
%   Algorithms
%
%   A driver .m file for the implementation of a Greedy algorithm
%   for optimising a ULA with N elements
%
%   Version: 0.2 - 22/04/2015
%
%   This function holds the cost) functions which
%   the genetic algorithms attempt to minimise
%   The function has two inputs as follows:
%
%   Inputs:
%   1) excite = excitation vector applied to the array
%   2) erBrokenOriginal - the original eR vector which represents the
%      original ULA with broken elements.
%
%   As well as a single output:
%   1) CF - The objective (cost) function value
%-------------------------------------------------------------------------%
%% /Fitness Function

% seperate excitation vector for algorithm
excite = [(excite(1) + excite(2)*1j)...
    (excite(3) + excite(4)*1j)...
    (excite(5) + excite(6)*1j)...
    (excite(7) + excite(8)*1j)...
    (excite(9) + excite(10)*1j)...
    (excite(11) + excite(12)*1j) ...
    (excite(13) + excite(14)*1j)...
    (excite(15) + excite(16)*1j)...
    (excite(17) + excite(18)*1j) ...
    (excite(19) + excite(20)*1j)...
    (excite(21) + excite(22)*1j)...
    (excite(23) + excite(24)*1j)];

er = er*excite'; % create resultant beampattern

% locating ratios between peaks of the current beam pattern
mainBeam = max(abs(er(190:210))); % set bounds in which Main Beam can be found
maxSidelobe1 = max(abs(er(1:175))); % set bounds in which max SL1 can be found + don't care region
maxSidelobe2 = max(abs(er(225:end))); % set bounds in which max SL2 can be found + don't care region
maxSidelobe = max([maxSidelobe1 maxSidelobe2]); % choose maximum sidelobe
beamRatio = maxSidelobe/mainBeam; % calculate corresponding Beam Ratio

%-- weights below are not used, but may be if cost function is more complex --%
% W1 = 1; % beamRatio weight
% W2 = 5; % maxSidelobe weight
% W3 = 1; % mainBeam weight

CF = beamRatio; % the Objective/Cost/Fitness Function

end

