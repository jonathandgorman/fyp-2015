function [er] = phasedArray(d,phi,N,nAngle)
%-------------------------------------------------------------------------%
% a function which calculates the array pattern for a phased array, given
% a number of array parameters as inputs. These inputs are outlined as
% follows:
%
% d = spacing (in terms of wavelength)
% phi = progressive phase diffence in radians
% N = number of elements in the array
% nangle = number angles between 0 and 2*pi defined
%
% a phased array is an array of antennas in which the relative phases of 
% the respective signals feeding the antennas are varied in such a way that 
% the effective radiation pattern of the array is reinforced in a desired 
% direction and suppressed in undesired directions.
%
%   University College Dublin (UCD)
%   School of Electronic, Electrical & Computer Engineering
%
%   Author: Jonathan Gorman
%   Student Number: 10310781
%
%-------------------------------------------------------------------------%
%% initialisation of local variables

thetaArray = []; % array containing theta values
er = zeros(1,nAngle); % initialises array for Er

for m = (1:(nAngle)) % loops through the number of angles
    
    theta = ((m-1)*pi)/(nAngle); % angles of theta
    thetaArray = [thetaArray theta];
    psi = 2*pi*d*cos(theta) + phi; % value of psi
    
    if (psi ~= 0) % if psi is not equal to zero
        
        er(m) = sin(N*psi/2)/sin(psi/2);
        
    else % otherwise if psi is equal to zero
        
        er(m) = 1;
        
    end
    
end % end for loop



end

