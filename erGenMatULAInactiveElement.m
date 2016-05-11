function [er] = erGenMatULAInactiveElement(d,phi,N,nAngle,thetaArray)
%-------------------------------------------------------------------------%
% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: The Application of Optimisation ALgorithms in Antenna
%   Beampattern Synthesis
%
%  A function that utilises MATLAB's repmat (repeat matrix) function in
%  order to build the matrices necessary for generating an array pattern
%  for a ULA with an inactive element number 6
%
%-------------------------------------------------------------------------%

% generates the matrix that determines the number of elements
elementComponent = repmat(((1:N)*1j),length(thetaArray),1);

% generates the matrix that determines the number of angles
thetaComponent = repmat((2*pi*d*cos(thetaArray(1:nAngle)')),1,length(1:N));

% generates the matrix that is the product of the above ( + phi)
exMatrix = (exp(elementComponent.*thetaComponent + phi));

% a replica of teh matrix above, but with a broken element
exMatrix2 = [exMatrix(:,(1:5)) zeros(length(exMatrix),1) exMatrix(:,(7:N))] ; % One prominent element missing

% the matrix to be outputted
er = exMatrix2;

end


