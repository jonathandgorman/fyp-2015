function [mainBeam, maxSidelobe,beamRatio] = peakFinderULA(er)
%-------------------------------------------------------------------------%
%   ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: The Application of Optimisation ALgorithms in Antenna
%   Beampattern Synthesis
%
% A function which is used to determine a number of parameters of an array
% pattern, which may be used as optimisation criteria. The function takes
% an array pattern (er) as an input and outputs a number of paramaters
% defining the array pattern, including:
%   
%   - Main Beam 
%   - Maximum Sidelobe 
%   - Beam Ratio 
%
%-------------------------------------------------------------------------%

mainBeam = max(abs(er(190:210))); % set bounds in which Main Beam can be found
maxSidelobe1 = max(abs(er(1:175))); % set bounds in which max SL1 can be found + don't care region
maxSidelobe2 = max(abs(er(225:end))); % set bounds in which max SL2 can be found + don't care region
maxSidelobe = max([maxSidelobe1 maxSidelobe2]); % choose maximum sidelobe
beamRatio = maxSidelobe/mainBeam; % calculate corresponding Beam Ratio

end

