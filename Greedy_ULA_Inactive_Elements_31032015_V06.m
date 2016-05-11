%% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%-------------------------------------------------------------------------%
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: Beam Pattern Synthesis in Sensor Arrays Using Optimisation
%   Algorithms
%
%   A driver .m file for the implementation of a Greedy algorithm
%   for optimising a ULA with N elements, where 1 element is inactive
%
%   Version: 0.3 - 22/04/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%-------------------------------------------------------------------------%
%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for k = 1:20 % loop for repeated output file generation
    
    %% Preamble
    format long % longer cmd window format
    clear all % clears all initialised variables
    close all % closes all open windows
    clc % clears the command window
    
    % Simulation start output
    disp('Simulation started - ULA With an Inactive Element Beam Pattern Synthesis Using a Greedy Algorithm')
    disp(date)
    disp(datestr(now, 'HH:MM:SS'))
    disp('-------------------------------------------------------------------')
    
    %% Declaration of array paramaters
    d = 0.5;  % spacing (in terms of wavelength)
    phi = 0; % progressive phase difference in radians
    N = 12;  % number of elements in the array
    nAngle = 400; % number of angles between 0 and pi
    
    R = 20; % relative distance between sidelobe level and mainbeam in dB, (e.g., R = 30)
    exciteInitial = 1.0*[1+0j (1+0j) (1+0j) (1+0j) (1+0j) (0+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j)]; % Initial excitations with element set as inactive
    BRArray = []; % array to hold changes in the beam ratio
    stdDevArray = []; % array to hold changes in the beam ratio
    thetaArray = zeros(1,nAngle); % array to hold values of theta
    aArray = []; % array to hold dolph-chebyshev weights
    
    seed = rand()*1e9; % generate random seed number
    rng(seed) % seed random number generator
    
    % for loop which generate theta and dolph-chebyshev values
    for m = 1:nAngle % cycle through angles
        theta = ((m-1)*pi)/nAngle; % angles of theta
        thetaArray(m) = theta;% storing the array elements
        
        
        %----------------- Code for producing Dolph-Chebyshev pattern ------------------%
        %     % checking whether to use doplh.m or dolph2.m fucntion to calculate
        %     % weights for dolph-chebysehev method
        %     if ((mod(N,2) ~= 0)&&(d <=0.49)) % if odd numebr of elements and d < 0.5
        %         [a,dph] = dolph2(d,rad2deg(theta),N,R); % dolph2.m function called, returns weights and beamwidth
        %     else
        %         [a,dph] = dolph(d,rad2deg(theta),N,R); % dolph.m function called, returns weights and beamwidth
        %     end
        %     aArray = [aArray;a]; % Array containing Dolph Chebyshev weights
        %-------------------------------------------------------------------------------%
        
    end
    
    % generate phased weights for ULA
    [er] = erGenMatULAInactiveElement(d,phi,N,nAngle,thetaArray);
    
    % subject broken array to initial excitation
    erInitial = er*exciteInitial';
    
    % find peaks of Er for array
    [ML, maxSL, BR] = peakFinderULAInactiveElement(erInitial)
    
    % store original beampattern parameters
    origBR = BR; % store original beam ratio
    origMaxSL = maxSL; % store original max sidelobe
    origML = ML; % store original main lobe
    
    % set initial parameters as initial optimum
    optBR= BR; % store initial beam ratio as optimum
    optMaxSL = maxSL; % store initial max side lobe as optimum
    optML= ML; % store initial main lobe as optimum
    
    exciteOptimum = exciteInitial; % setting initial excitation current vector as the optimum
    
    %% Greedy Algorithm - ULA with an Inactive Element
    % algorithm stopping criteria variables - N.B. Should set variables
    maxIters = 1e7; % max number of iterations 
    maxTime = 10; % max time limit in seconds - for time based experiments
    minBR = 0.075; % minimum target BR - for perfromance based experiments
   
    tic % start timer
    for n = 1:maxIters
        
        %-- 0nly one of the two switch statements below should be active at any time--%
        % switch statement that controls stdDev variable according to
        % simulation time
        switch  logical(true) % switch according to the statement that is true
            case (toc < (maxTime/8))
                stdDev = 0.5; % standard deviation 1
            case ((toc > (maxTime/8)) && (toc < (maxTime/7)))
                stdDev = 0.5; % standard deviation 2
            case ((toc > (maxTime/7)) && (toc < (maxTime/6)))
                stdDev = 0.25; % standard deviation 3
            case ((toc > (maxTime/6)) && (toc < (maxTime/5)))
                stdDev = 0.25; % standard deviation 4
            case ((toc > (maxTime/5)) && (toc < (maxTime/4)))
                stdDev = 0.1; % standard deviation 5
            case ((toc > (maxTime/4)) && (toc < (maxTime/3)))
                stdDev = 0.05; % standard deviation 6
            case ((toc > (maxTime/3)) && (toc < (maxTime/2)))
                stdDev = 0.01; % standard deviation 7
            otherwise
                stdDev = 0.005; % standard deviation 8
        end
        
        stdDevArray = [stdDevArray stdDev];
        
        % switch statement that controls stdDev variable according to
        % current iteration value
        %         switch  logical(true) % switch according to the statement that is true
        %             case (n < (maxIters/8))
        %                 stdDev = 0.5 % standard deviation 1
        %             case ((n > (maxIters/8)) && (n < (maxIters/7)))
        %                 stdDev = 0.5 % standard deviation 2
        %             case ((n > (maxIters/7)) && (n < (maxIters/6)))
        %                 stdDev = 0.25 % standard deviation 3
        %             case ((n > (maxIters/6)) && (n < (maxIters/5)))
        %                 stdDev = 0.25 % standard deviation 4
        %             case ((n > (maxIters/5)) && (n < (maxIters/4)))
        %                 stdDev = 0.1 % standard deviation 5
        %             case ((n > (maxIters/4)) && (n < (maxIters/3)))
        %                 stdDev = 0.05 % standard deviation 6
        %             case ((n > (maxIters/3)) && (n < (maxIters/2)))
        %                 stdDev = 0.01 % standard deviation 7
        %             otherwise
        %                 stdDev = 0.05 % standard deviation 8
        %         end
        
        gauss = stdDev*(randn(1,N)+randn(1,N)*1j); % create a new gaussian random perturbation
        gauss(6) = 0; % setting 6th element equal to zero
        exciteStore = exciteOptimum; % optimum excitation stored
        exciteGauss = exciteOptimum + gauss; % gaussian random perturbation added to optimum excitation found
        
        %------------------------------ UNCOMMENT FOR BOUNDED EXPERIMENTS --------------------------------%
        exciteGauss((imag(exciteGauss) > 5)) = (real(exciteGauss((imag(exciteGauss) > 5)))+ 5j); % if imag component is greater than +5, set it to + 5
        exciteGauss((real(exciteGauss) > 5)) = (5 + imag(exciteGauss((real(exciteGauss) > 5)))*1j); % if real component is greater than +5, set it to + 5
        exciteGauss((imag(exciteGauss) < -5)) = (real(exciteGauss((imag(exciteGauss) < -5)))+ -5j); % if imag component is less than -5, set it to - 5
        exciteGauss((real(exciteGauss) < -5)) = (-5 + imag(exciteGauss((real(exciteGauss) < -5)))*1j); % if real component is less than -5, set it to - 5
        %------------------------------ UNCOMMENT FOR BOUNDED EXPERIMENTS --------------------------------%
        
        erGreedy = er*exciteGauss'; % resultant beampattern following addition of random perturbation
        [ML, maxSL, BR] = peakFinderULAInactiveElement(erGreedy); % array parameters calculated for new array pattern
        
        % current array parameters stored
        currBR= BR; % store current value of BR
        currMaxSL = maxSL; % store current value of maxSL
        currMB = ML; % store current value of MB
        
        if (currBR < optBR)% if beam ratio has been reduced
            
            optBR = currBR;% current BR value is set as optimum
            optMaxSL = currMaxSL;% current max SL value is set as optimum
            optML = currMB; % current MB value is set as optimum
            
            exciteOptimum = exciteGauss; % the optimum excitation current
            
            %          % update plot to show algorithm progress
            %          drawnow
            %          plot((20*log10(abs(erPhased*exciteRand')/max(abs(erPhased*exciteRand')))),'k')
            %          grid on
            %          xlabel('Angle \theta (degrees)','FontSize', 25)
            %          ylabel('Er Magnitude (dB)','FontSize', 25)
            %          legend('Greedy Optimised ULA')
            %          axis([0 400 -50 10])
            %
            
        end
        
        BRArray = [BRArray optBR]; % array storing the changes in the BR
        
        % check to see if maximum set time has elapsed
        if ((toc >= maxTime) || (optBR <= minBR))% if current time is greater than time limit
            iters = n; % store corresponding iteration count
            endTime = toc;  % store corresponding finishing time
            break % bvreak out of loop
        end       
        
    end % end of Greedy algorithm
    
    erGreedy = er*exciteOptimum'; % set optimised excitation current vector
    
    %% Greedy Algorithm ULA Simulation Plots
    % plot change in CF against iterations
    figure
    plot(BRArray,'k')
    grid on
    xlabel('Iterations (n)','FontSize', 25)
    ylabel('Cost Function','FontSize', 25)
    hLegend = legend('Change in Cost Function');
    set(hLegend,'FontSize',15);
    
    % change in standard deviation versus iterations
    figure
    plot(stdDevArray,'k')
    grid on
    ylim([0 0.5])
    xlabel('Iterations (n)','FontSize', 25)
    ylabel('Sigma','FontSize', 25)
    hLegend = legend('Change in the Value of Sigma');
    set(hLegend,'FontSize',15);
    
    % linear plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),abs(erInitial),'r--')
    hold on
    plot(rad2deg(thetaArray),abs(erGreedy),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Gain (dB)','FontSize', 25)
    hLegend = legend('Initial Beampattern','Greedy Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % decibel plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),(20*log10(abs(erGreedy)/max(abs(erGreedy)))),'b-.')
    hold on
    plot(rad2deg(thetaArray),(20*log10(abs(erInitial)/max(abs(erInitial')))),'r--')
    ylim([-40,5])
    grid on
    
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Gain (dB)','FontSize', 25)
    hLegend = legend('12 element ULA - element 6 inactive','Greedy Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    plot(rad2deg(thetaArray),(20*log10(abs(erGreedy)/max(abs(erGreedy)))),'k')
    
    % polar plot comparing the initial beampattern with the optimised beampattern
    figure
    polar(thetaArray,abs(erInitial)','r--')
    hold on
    polar(thetaArray,abs(erGreedy)','k')
    grid on
    hLegend = legend('Initial Beampattern','Greedy Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % plot of optimised excitation current distribution
    figure
    plot(exciteOptimum,'*')
    hold on
    plot(1,0,'rx')
    axis square
    grid on
    xlabel('Real Component','FontSize', 25)
    ylabel('Imaginary Component','FontSize', 25)
    hLegend = legend('Optimised Excitation Current Values','Initial Excitation Current Values');
    set(hLegend,'FontSize',15);
    
    %% Statistics and Results
    % display simulation statistics and results
    disp('The optimised excitation current vector:')
    disp(exciteOptimum')
    disp('The number of iterations performed:')
    disp(iters)
    disp('The change in the Main Lobe:')
    disp(optML - origML)
    disp('The reduction in the Max Side Lobe in decibels:')
    disp(((20*log10(abs(optMaxSL)/max(abs(erGreedy)))))+ abs((20*log10(abs(origMaxSL)/max(abs(erInitial))))))
    disp('The reduction in the Beam Ratio:')
    disp(optBR - origBR)
    disp('The number of changes in the beam ratio:')
    disp(numel(unique(BRArray)))
    disp('The norm of the optimised excitation current distribution:')
    disp(norm(exciteOptimum))
    disp('-------------------------------------------------------------------')
    
    %% End of simulation
    str1 = sprintf('./Greedy_ULA_Inactive_Elements_matfiles/BR_%d_Greedy_ULA_Inactive_Element_iters_%d_time_%d',optBR,iters,endTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'iters','endTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','erGreedy','erInitial','exciteOptimum','thetaArray','BRArray','stdDevArray')
    
    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display time
    % end of simulation
    
% end

