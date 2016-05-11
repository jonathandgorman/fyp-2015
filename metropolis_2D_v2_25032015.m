%% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%-------------------------------------------------------------------------%
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: Beam Pattern Synthesis in Sensor Arrays Using Optimisation
%   Algorithms
%
%   A driver .m file for the implementation of a Metropolis algorithm
%   for optimising a ULA with N elements
%
%   Version: 0.3 - 17/04/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%   3) added additional stopping criteria and current bounds for final thesis
%-------------------------------------------------------------------------%

%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for k = 1:20 % loop for repeated output file generation

    % Preamble
    format long % longer cmd window format
    clear all % clears all initialised variables
    close all % closes all open windows
    clc % clears the command window
    
    % Simulation start output
    disp('Simulation started - WSN Beam Pattern Synthesis Using a Metropolis Algorithm')
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
    xlabel('d (wavelengths)','FontSize',25)
    ylabel('d (wavelengths)','FontSize',25)
    legend('Network Origin', 'Network Sensor','Network Boundry','Receiver')
    grid on
    ylim([-(R+1.0) (R+1.0)])
    xlim([-(R+1.0) (R+1.0)])
    
    % randomising random number generator for remainder simulation
    seed = rand()*1e9; % randomise seed
    rng(seed + 1e1) % random number generator seeded
    
    %% Generate beampattern from sensor elements
    nAngle = 400; % number angles between 0 and 2*pi defined
    BRArray = []; % array to hold changes in the beam ratio
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
    
    % calculate initial phased matrix
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
    
    optBR = origBR; % CF with current excitations
    currBR = origBR; % CF with current excitations
    exciteCurrent = exciteInitial; % set current excitation vector as the initial
    exciteOptimum = exciteInitial; % set current optimum excitation vector as the initial
    
    %% Metropolis Algorithm
    % Creating temperature profile
    maxIters = 20000; % maximum number of iterations defined
    maxTime = 10; % max time in seconds
    T0 = 0.2; % initial max temperature
    Tfinal = 0.0001; % final minimum temperature
    Tarray = [T0 zeros(1,maxIters-1)]; % array of temperature valeus
    alpha = nthroot(Tfinal,maxIters); % value of alpha
    minBR = 0.2; % minimum BR for performance measurments
    
    for n = 2:maxIters
        Tarray(n) = Tarray(n-1)*alpha;
    end
    
    probArray = []; % initialise array holding probabilities of acceptance values
    
    tic % start timer
    for n = 1:maxIters
        
        gauss = Tarray(n)*(randn(1,N) + randn(1,N)*1j); % create random excitation according to the standard deviation
        exciteGauss = exciteCurrent + gauss; % add random excitation to current vector
        
        %         exciteGauss((imag(exciteGauss) > 5)) = (real(exciteGauss((imag(exciteGauss) > 5)))+ 5j); % if imag component is greater than +5, set it to + 5
        %         exciteGauss((real(exciteGauss) > 5)) = (5 + imag(exciteGauss((real(exciteGauss) > 5)))*1j); % if real component is greater than +5, set it to + 5
        %         exciteGauss((imag(exciteGauss) < -5)) = (real(exciteGauss((imag(exciteGauss) < -5)))+ -5j); % if imag component is less than -5, set it to - 5
        %         exciteGauss((real(exciteGauss) < -5)) = (-5 + imag(exciteGauss((real(exciteGauss) < -5)))*1j); % if real component is less than -5, set it to - 5
        %
        erMetropolis = phasorMatrix*exciteGauss';  % generate of beampattern with random excitation
        
        [ML, maxSL, BR] = peakFinder2D(erMetropolis); % array parameters calculated for new array pattern
        gaussBR = BR; % find corresponding cost function
        gaussMaxSL = maxSL; % find corresponding cost function
        gaussML = ML; % find corresponding cost function
        
        if (gaussBR < optBR) % if random CF is better than optimum CF found
            
            exciteOptimum = exciteGauss; % store optimum excitation vector
            optBR = gaussBR; % set random CF as new optimum CF
            optMaxSL = gaussMaxSL; % set random Max SL as new optimum maxSL
            optML = gaussML; % set random ML as new optimum ML
            
            %--- continually draws plot, commented out for efficiency ---%
            %         drawnow
            %         plot((20*log10(abs(erMetropolis)/max(abs(erMetropolis)))))
            %         axis([0 400 -45 0])
            
        end
        
        prob = exp((-(currBR - optBR)/Tarray(n))); % calculate out probability of acceptance
        probArray = [probArray prob]; % add value to array of probabilities of acceptance
        BRArray = [BRArray optBR]; % add optimum BR value to array
        
        if (gaussBR < currBR) % if random CF is better than the current excitation vector, accept it
            
            exciteCurrent = exciteGauss; %  set current excitation to random excitation
            currBR = gaussBR; % set current CF to Gauss CF
            
        elseif (prob > rand) % Otherwise, take it with probability of acceptance
            
            exciteCurrent = exciteGauss; %  set current excitation to random excitations
            currBR = gaussBR; % set current CF to Gauss CF
            
        else % otherwise
            
            exciteCurrent; % keep the current excitation vector
            
        end
        
        % check to see if maximum set time has elapsed
        if ((toc >= maxTime) || (optBR <= minBR)) % if current time is greater than time limit or minimum BR met
            iters = n; % store corresponding iteration count
            endTime = toc;  % store corresponding finishing time
            break % bvreak out of loop
        end
        
    end
    
    erMetropolis = phasorMatrix*exciteOptimum'; % set optimised excitation current vector
    
    %% Metropolis Algorithm ULA Simulation Plots
    % plot change in CF against iterations
    figure
    plot(BRArray,'k')
    grid on
    xlabel('Iterations (n)','FontSize', 25)
    ylabel('Cost Function','FontSize', 25)
    hLegend = legend('Change in Cost Function');
    set(hLegend,'FontSize',15);
    
    % change in probability of acceptance versus iterations
    figure
    plot(probArray,'k')
    grid on
    xlabel('Iterations (n)','FontSize', 25)
    ylabel('Probability of Acceptance','FontSize', 25)
    hLegend = legend('Probability of Acceptance');
    set(hLegend,'FontSize',15);
    
    % Temperature profile of Metropolis Algorithm
    figure
    plot(Tarray,'k')
    grid on
    xlabel('Iterations (n)','FontSize', 25)
    ylabel('Temperature','FontSize', 25)
    hLegend = legend('Temperature Profile');
    set(hLegend,'FontSize',15);
    
    % linear plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),abs(erInitial),'r--')
    hold on
    plot(rad2deg(thetaArray),abs(erMetropolis),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Gain (dB)','FontSize', 25)
    hLegend = legend('Initial Beampattern','Metropolis Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % decibel plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),(20*log10(abs(erInitial)/max(abs(erInitial)))),'r--')
    ylim([-40,5])
    hold on
    plot(rad2deg(thetaArray),(20*log10(abs(erMetropolis)/max(abs(erMetropolis)))),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Er Magnitude (dB)','FontSize', 25)
    hLegend = legend('Initial Beampattern','Metropolis Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % polar plot comparing the initial beampattern with the optimised beampattern
    figure
    polar(thetaArray,abs(erInitial)','r--')
    hold on
    polar(thetaArray,abs(erMetropolis)','k')
    grid on
    hLegend = legend('Initial Beampattern','Metropolis Optimised Beampattern');
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
    disp(((20*log10(abs(optMaxSL)/max(abs(erMetropolis)))))+ abs((20*log10(abs(origMaxSL)/max(abs(erInitial))))))
    disp('The reduction in the Beam Ratio:')
    disp(optBR - origBR)
    disp('The number of changes in the Beam Ratio:')
    disp(numel(unique(BRArray)))
    disp('The norm of the optimised excitation current distribution:')
    disp(norm(exciteOptimum))
    disp('-------------------------------------------------------------------')
    
    %% End of simulation
    str1 = sprintf('./Metropolis_2D_matfiles/unBoundedPerfromanceBR_%d_Metropolis_2D_iters_%d_time_%d',optBR,iters,endTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'iters','endTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','erMetropolis','erInitial','exciteOptimum','thetaArray','BRArray','probArray','Tarray')
    
    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display
    % end of simulation
    
% end

