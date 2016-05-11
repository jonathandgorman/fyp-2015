%% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%-------------------------------------------------------------------------%
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: The Application of Optimisation Algorithms in Antenna
%   Beampattern Synthesis
%
%   A driver .m file for the implementation of a Greedy algorithm
%   for optimising a ULA with N elements
%
%   Version: 0.5 - 01/04/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%-------------------------------------------------------------------------%
%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for k = 1:20 % loop for repeated output file generation
    
    % Preamble
    format long % longer cmd window format
    clear all % clears all initialised variables
    close all % closes all open windows
    clc % clears the command window
    % Simulation start output
    disp('Simulation started - WSN Beam Pattern Synthesis Using a Greedy Algorithm')
    disp(date)
    disp(datestr(now, 'HH:MM:SS'))
    disp('-------------------------------------------------------------------')
    
    %%-- comment out the below code in order to create a random pattern --%
    seedInitial = 11; % set seed for random number generation of the network
    rng(seedInitial) % seed the random number generator
    
    %% Creating a circle and populating with random point sensor coordinates
    % Creating the circle and the random points
    R = 2.0; % radius of circle
    xCentre = 0.0; % x-coordinate of centre
    yCentre = 0.0; % y-coordinate of centre
    N = 32; % number of points
    r = R*sqrt(rand(N,1)); % random distance from origin
    t = 2*pi*rand(N,1); % random angle
    x = xCentre + r.*cos(t); % x-coordinate of point
    y = yCentre + r.*sin(t);% y-coordinate of point
    XRX = -2.5; % x-coordinate of the receiver
    YRX = 0; % y-coordinate of the receiver
    t2 = linspace(0,2*pi); % creates equally spaced points for 360 degrees
    X = xCentre + R*cos(t2); % creates X coordinates for circle
    Y = yCentre + R*sin(t2); % creates Y coordinates for circle
    
    % sensor plot for 32 WSN array
    plot(xCentre,yCentre,'r+',X,Y,'k',x,y,'b*',XRX,YRX,'m^') % plots centre, circle and point coordinates
    axis equal % axis equal
    xlabel('\it{d} (wavelengths)','FontSize',15)
    ylabel('\it{d} (wavelengths)','FontSize',15)
    hlegend = legend('Network Origin', 'Network Sensor','Network Boundry','Receiver');
    set(hlegend,'FontSize',15)
    grid on
    ylim([-(R+1.0) (R+1.0)])
    xlim([-(R+1.0) (R+1.0)])
    
    % randomising random number generator for remainder simulation
    seed = rand()*1e9; % randomise seed
    rng((seed + rand()*1e6)) % random number generator seeded
    
    %% Generate beampattern from sensor elements
    nAngle = 400; % number angles between 0 and 2*pi defined
    BRArray = []; % array to hold changes in the beam ratio
    stdDevArray = []; % array to hold changes in the beam ratio
    thetaArray = zeros(1,nAngle); % array to hold values of theta
    aArray = []; % array to hold dolph-chebyshev weights
    
    % generate angles of theta
    thetaArray = zeros(1, nAngle); % initialise array to hold angles
    for m = 1:nAngle % loop through number of angles
        
        theta = ((m)*2*pi)/nAngle; % angles of theta for the array
        thetaArray(m) = theta; % storing the array elements
        
    end
    
    % generate phasor for each element at each angle of theta
    [th,rad] = cart2pol(x,y); % change sensor coordinates to polar form
    phasorMatrix = zeros(length(thetaArray),N); % create matrix to hold phasors
    
    for k = (1:length(thetaArray)) % loop through angles
        for n = (1:N)
            phasor = exp(-1*1j*2*pi*rad(n)*cos(thetaArray(k) - th(n))); % caclulate phasor
            phasorMatrix(k,n) = phasor; % store phasor
        end
    end
    
    % generate excitation vector such that resultant beam is steered at 180 degrees
    exciteInitial = conj(phasorMatrix((400),:))*1.0;
    
    % calculate initial beampattern
    erInitial = phasorMatrix*exciteInitial';
    
    % find peaks of initial beampattern
    [ML, maxSL,BR] = peakFinder2D(erInitial);
    
    % store original beampattern parameters
    origBR = BR; % store original beam ratio
    origMaxSL = maxSL; % store original max sidelobe
    origML = ML; % store original main lobe
    
    % set initial parameters as initial optimum
    optBR= BR; % store initial beam ratio as optimum
    optMaxSL = maxSL; % store initial max side lobe as optimum
    optML= ML; % store initial main lobe as optimum
    
    exciteOptimum = exciteInitial; % setting initial excitation current vector as the optimum
    
    %% Greedy Algorithm - 2D WSN
    % algorithm stopping criteria variables
    maxIters = 1e7; % max number of iterations
    maxTime = 10; % max time in seconds - for time based testing
    minBR = 0.02; % minimum BR for performance based testing             
    tic % start timer
    
    for n = 1:maxIters
        
        %-- 0nly one of the two switch statements below should be active at any time--%
        % switch statement that controls stdDev variable according to
        % simulation time
        switch  logical(true) % switch according to the statement that is true
            case (toc < (maxTime/8))
                stdDev = 0.5; % standard deviation 1
            case ((toc > (maxTime/8)) && (toc < (maxTime/7)))
                stdDev = 0.4; % standard deviation 2
            case ((toc > (maxTime/7)) && (toc < (maxTime/6)))
                stdDev = 0.25; % standard deviation 3
            case ((toc > (maxTime/6)) && (toc < (maxTime/5)))
                stdDev = 0.2; % standard deviation 4
            case ((toc > (maxTime/5)) && (toc < (maxTime/4)))
                stdDev = 0.1; % standard deviation 5
            case ((toc > (maxTime/4)) && (toc < (maxTime/3)))
                stdDev = 0.05; % standard deviation 6
            case ((toc > (maxTime/3)) && (toc < (maxTime/2)))
                stdDev = 0.01; % standard deviation 7
            otherwise
                stdDev = 0.005; % standard deviation 8
        end
        
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
        
        stdDevArray = [stdDevArray stdDev]; % populate array to hold change in standard deviation variables
        
        gauss = stdDev*(randn(1,N) + randn(1,N)*1j); % create random excitation according to the standard deviation
        exciteStore = exciteOptimum; % optimum excitation stored
        exciteGauss = exciteOptimum + gauss; % gaussian random perturbation added to optimum excitation found
        
        %------------------------------ UNCOMMENT FOR BOUNDED EXPERIMENTS --------------------------------%
        %          exciteGauss((imag(exciteGauss) > 5)) = (real(exciteGauss((imag(exciteGauss) > 5)))+ 5j); % if imag component is greater than +5, set it to + 5
        %          exciteGauss((real(exciteGauss) > 5)) = (5 + imag(exciteGauss((real(exciteGauss) > 5)))*1j); % if real component is greater than +5, set it to + 5
        %          exciteGauss((imag(exciteGauss) < -5)) = (real(exciteGauss((imag(exciteGauss) < -5)))+ -5j); % if imag component is less than -5, set it to - 5
        %          exciteGauss((real(exciteGauss) < -5)) = (-5 + imag(exciteGauss((real(exciteGauss) < -5)))*1j); % if real component is less than -5, set it to - 5
        %------------------------------ UNCOMMENT FOR BOUNDED EXPERIMENTS --------------------------------%        
        
        erGreedy = phasorMatrix*exciteGauss'; % resultant beampattern following addition of random perturbation
        [ML, maxSL, BR] = peakFinder2D(erGreedy); % array parameters calculated for new array pattern
        
        % current array parameters stored
        currBR= BR; % store current value of BR
        currMaxSL = maxSL; % store current value of maxSL
        currMB = ML; % store current value of MB
        
        if (currBR < optBR)% if beam ratio has been reduced
            
            optBR= currBR;% current BR value is set as optimum
            optMaxSL = currMaxSL;% current max SL value is set as optimum
            optML= currMB; % current MB value is set as optimum
            
            exciteOptimum = exciteGauss; % the optimum excitation current
            
            %          % update plot to show algorithm progress
            %          drawnow
            %          plot((20*log10(abs(erPhased*exciteRand')/max(abs(erPhased*exciteRand')))),'k')
            %          grid on
            %          xlabel('Angle \theta (degrees)','FontSize', 25)
            %          ylabel('Er Magnitude (dB)','FontSize', 25)
            %          legend('Greedy Optimised ULA')
            %          axis([0 400 -50 10])
            
        end
        
        BRArray = [BRArray optBR]; % array storing the changes in the BR
        
        % check to see if maximum set time has elapsed
        if ((toc >= maxTime) || (optBR <= minBR)) % if current time is greater than time limit or BR less than specified minimum targer
            iters = n; % store corresponding iteration count
            endTime = toc;  % store corresponding finishing time
            break % break out of loop
        end
        
    end % end of greedy algorithm
    
    erGreedy = phasorMatrix*exciteOptimum'; % set optimised excitation current vector
    
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
    ylim([0 0.55])
    xlabel('Iterations (m)','FontSize', 25)
    ylabel('Sigma','FontSize', 25)
    hLegend = legend('Change in the Value of Sigma');
    set(hLegend,'FontSize',15);
    
    % linear plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),abs(erInitial),'k')
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
    plot(rad2deg(thetaArray),(20*log10(abs(erInitial)/max(abs(erInitial)))),'k')
    ylim([-40,5])
    xlim([0,360])
    hold on
    plot(rad2deg(thetaArray),(20*log10(abs(erGreedy)/max(abs(erGreedy)))),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Gain (dB)','FontSize', 25)
    hLegend = legend('Initial Beam Pattern of 32 element WSN array','Greedy Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
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
    % save to file
    str1 = sprintf('./Greedy_2D_matfiles/unboundedPerformanceBR_%d_Greedy_2D_iters_%d_time_%d',optBR,iters,endTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'iters','endTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','erGreedy','erInitial','exciteOptimum','thetaArray','BRArray','stdDevArray')
    
    disp('-------------------------------------------------------------------')
    toc % stops clock timer
    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display time
    % end of simulation
    
% end
