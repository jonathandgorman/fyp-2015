function [CF] = fitness_function_GA_ULA_Inactive_Element(excite,erBroken)
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: Beam Pattern Synthesis in Sensor Arrays Using Optimisation
%   Algorithms
%
%   This function holds the cost functions which
%   the genetic algorithm minimises.
%   This function is also capable of finding the peaks of the a beampattern,
%   which are used as part of the
%   cost function. The function has two inputs as follows:
%
%   Inputs:
%   1) excite = excitation vector applied to the array
%   2) erBroken - the original eR vector which represents the
%      original ULA with broken elements.
%
%   As well as a single output:
%   1) CF - The cost function value
%-------------------------------------------------------------------------%
%% Objective/Cost/Fitness Function
excite = [(excite(1) + excite(2)*1j)...
    (excite(3) + excite(4)*1j)...
    (excite(5) + excite(6)*1j)...
    (excite(7) + excite(8)*1j)...
    (excite(9) + excite(10)*1j)...
    (0 + 0*1j) ...
    (excite(13) + excite(14)*1j)...
    (excite(15) + excite(16)*1j)...
    (excite(17) + excite(18)*1j) ...
    (excite(19) + excite(20)*1j)...
    (excite(21) + excite(22)*1j)...
    (excite(23) + excite(24)*1j)];

er = erBroken*excite'; % create beampattern

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

CF = beamRatio; % the objective (cost) function

end

