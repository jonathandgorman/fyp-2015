function [er] = erGenMatULA(d,phi,N,nAngle,thetaArray)
%-------------------------------------------------------------------------%
% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
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
%  A function that utilises MATLAB's repmat (repeat matrix) function in
%  order to build the matrices necessary for generating an array pattern
%  for a Uniform Linear Array.
%
%-------------------------------------------------------------------------%

elementComponent = repmat(((1:N)*1j),length(thetaArray),1);% generates the matrix that determines the number of elements
thetaComponent = repmat((2*pi*d*cos(thetaArray(1:nAngle)')),1,length(1:N));% generates the matrix that determines the number of angles
er = (exp(elementComponent.*thetaComponent + phi));% generates the matrix that is the product of the above ( + phi)

end


