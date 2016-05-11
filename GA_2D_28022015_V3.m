%% ME Electronic & Computer Engineering Final Year Project (EEEN40240)
%-------------------------------------------------------------------------%
%   University College Dublin (UCD)
%   School of Electrical, Electronic & Communications Engineering
%
%   Author: Jonathan D. Gorman
%   Project: Beam Pattern Synthesis in Sensor Arrays Using Optimisation
%   Algorithms
%
%   A driver .m file for the implementation of a Genetic algorithm
%   for optimising a ULA with N elements
%
%   Version: 0.3 - 01/04/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%-------------------------------------------------------------------------%

%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for n = 1:20 % loop for repeated output file generation
    
    % Preamble
    format long % longer cmd window format
    clear all % clears all initialised variables
    close all % closes all open windows
    clc % clears the command window
    
    % Simulation start output
    disp('Simulation started - WSN Array Beam Pattern Synthesis Using a Greedy Algorithm')
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
    rng(seed + 1e2) % random number generator seeded
    
    %% Generate contribution from array elements
    nAngle = 400; % number angles between 0 and 2*pi defined
    
    % generate angles of theta
    thetaArray = zeros(1, nAngle);
    for m = 1:nAngle
        
        theta = ((m)*2*pi)/nAngle; % angles of theta for the array
        thetaArray(m) = theta; % storing the array elements
        
    end
    
    % generate phasor for each element at each angle of theta
    [th,rad] = cart2pol(x,y); % change sensor coordinates to polar form
    phasorMatrix = zeros(length(thetaArray),N); % create matrix to hold phasors
    
    for k = (1:length(thetaArray)) % loop through angles
        for n = (1:N)
            phasor = exp(-1j*2*pi*rad(n)*cos(thetaArray(k) - th(n)));
            phasorMatrix(k,n) = phasor;
        end
    end
    
    exciteInitial = (conj(phasorMatrix((400),:)));% generate initial excitation vector
    erInitial = phasorMatrix*exciteInitial'; % calculate initial beampattern
    [ML, maxSL,BR] = peakFinder2D(erInitial);    % find peaks of initial beampattern
     
    % store original beampattern parameters
    origBR = BR; % store original beam ratio
    origMaxSL = maxSL; % store original max sidelobe
    origML = ML; % store original main lobe
    
    % set initial parameters as initial optimum
    optBR= BR; % store initial beam ratio as optimum
    optMaxSL = maxSL; % store initial max side lobe as optimum
    optML= ML; % store initial main lobe as optimum
    
    exciteOptimum = exciteInitial; % setting initial excitation current vector as the optimum
    
    %% Genetic Algorithm
    % algorithm parameters and stopping criteria
    genNumber = inf; % number of generations
    maxTime = 10; % max time limit - for time based testing
    minBR = 0.2; % minimum CF target - for perfromance based testing
    
    % splitting the excitation vector into real and imaginary parts so that
    % it may be manipulated by the genetic algorithm
    exciteImag = imag(exciteInitial);
    exciteReal = real(exciteInitial);
    exciteArraySplit = ...
        [exciteReal(1) exciteImag(1) ...
        exciteReal(2) exciteImag(2) ...
        exciteReal(3) exciteImag(3) ...
        exciteReal(4) exciteImag(4) ...
        exciteReal(5) exciteImag(5) ...
        exciteReal(6) exciteImag(6) ...
        exciteReal(7) exciteImag(7) ...
        exciteReal(8) exciteImag(8) ...
        exciteReal(9) exciteImag(9) ...
        exciteReal(10) exciteImag(10) ...
        exciteReal(11) exciteImag(11) ...
        exciteReal(12) exciteImag(12) ...
        exciteReal(13) exciteImag(13) ...
        exciteReal(14) exciteImag(14) ...
        exciteReal(15) exciteImag(15) ...
        exciteReal(16) exciteImag(16) ...
        exciteReal(17) exciteImag(17) ...
        exciteReal(18) exciteImag(18) ...
        exciteReal(19) exciteImag(19) ...
        exciteReal(20) exciteImag(20) ...
        exciteReal(21) exciteImag(21) ...
        exciteReal(22) exciteImag(22) ...
        exciteReal(23) exciteImag(23) ...
        exciteReal(24) exciteImag(24) ...
        exciteReal(25) exciteImag(25) ...
        exciteReal(26) exciteImag(26) ...
        exciteReal(27) exciteImag(27) ...
        exciteReal(28) exciteImag(28) ...
        exciteReal(29) exciteImag(29) ...
        exciteReal(30) exciteImag(30) ...
        exciteReal(31) exciteImag(31) ...
        exciteReal(32) exciteImag(32) ...
        ];
    
    ub = 5*[ones(1,64)]; % upper bound of excitations
    lb = -5*[ones(1,64)]; % lower bound of excitations
    
    FitnessFunction = @(exciteArraySplit) fitness_function_2D_GA(exciteArraySplit,phasorMatrix);% set initial cost function
    % GA options set - N.B: Add the following for dynamic plots 'PlotFcns',{@gaplotbestf,@gaplotstopping}
    options = gaoptimset('StallGenLimit',9999999999999,'Generations',genNumber,'timeLimit',maxTime,'FitnessLimit',minBR);
    tic
    [exciteOptimum,fitnessFunctionEval,exitFlag,output] = ga(FitnessFunction,64,[],[],[],[],[],[],[],options) % calling the GA algorithm
    endTime = toc;
    % print to screen options linked to simulannealbnd
    fprintf('The number of generations was : %d\n', output.generations);
    fprintf('The number of function evaluations was : %d\n', output.funccount);
    fprintf('The best function value found was : %g\n', fitnessFunctionEval);
    
    % store number of generations completed
    genNumber = output.generations;
    
    % recombining initial exciattion vector
    exciteOptimum = [(exciteOptimum(1) + exciteOptimum(2)*1j)...
        (exciteOptimum(3) + exciteOptimum(4)*1j)...
        (exciteOptimum(5) + exciteOptimum(6)*1j)...
        (exciteOptimum(7) + exciteOptimum(8)*1j)...
        (exciteOptimum(9) + exciteOptimum(10)*1j)...
        (exciteOptimum(11) + exciteOptimum(12)*1j) ...
        (exciteOptimum(13) + exciteOptimum(14)*1j)...
        (exciteOptimum(15) + exciteOptimum(16)*1j)...
        (exciteOptimum(17) + exciteOptimum(18)*1j) ...
        (exciteOptimum(19) + exciteOptimum(20)*1j)...
        (exciteOptimum(21) + exciteOptimum(22)*1j)...
        (exciteOptimum(23) + exciteOptimum(24)*1j)...
        (exciteOptimum(25) + exciteOptimum(26)*1j)...
        (exciteOptimum(27) + exciteOptimum(28)*1j)...
        (exciteOptimum(29) + exciteOptimum(30)*1j)...
        (exciteOptimum(31) + exciteOptimum(32)*1j)...
        (exciteOptimum(33) + exciteOptimum(34)*1j)...
        (exciteOptimum(35) + exciteOptimum(36)*1j)...
        (exciteOptimum(37) + exciteOptimum(38)*1j)...
        (exciteOptimum(39) + exciteOptimum(40)*1j)...
        (exciteOptimum(41) + exciteOptimum(42)*1j)...
        (exciteOptimum(43) + exciteOptimum(44)*1j)...
        (exciteOptimum(45) + exciteOptimum(46)*1j)...
        (exciteOptimum(47) + exciteOptimum(48)*1j)...
        (exciteOptimum(49) + exciteOptimum(50)*1j)...
        (exciteOptimum(51) + exciteOptimum(52)*1j)...
        (exciteOptimum(53) + exciteOptimum(54)*1j)...
        (exciteOptimum(55) + exciteOptimum(56)*1j)...
        (exciteOptimum(57) + exciteOptimum(58)*1j)...
        (exciteOptimum(59) + exciteOptimum(60)*1j)...
        (exciteOptimum(61) + exciteOptimum(62)*1j)...
        (exciteOptimum(63) + exciteOptimum(64)*1j)];
    
    % find peaks of Er for array
    [ML, maxSL, BR] = peakFinder2D(phasorMatrix*exciteOptimum');
    
    % store original beampattern parameters
    optBR = BR;
    optMaxSL = maxSL;
    optML = ML;
    
    %% Genetic Algorithm ULA Simulation Plots
    % linear plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),abs(erInitial),'r--')
    hold on
    plot(rad2deg(thetaArray),abs(phasorMatrix*exciteOptimum'),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Gain (dB)','FontSize', 25)
    hLegend = legend('Initial Beampattern','Genetic Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % decibel plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),(20*log10(abs(erInitial)/max(abs(erInitial)))),'r--')
    ylim([-40,5])
    hold on
    plot(rad2deg(thetaArray),(20*log10(abs(phasorMatrix*exciteOptimum')/max(abs(phasorMatrix*exciteOptimum')))),'k')
    grid on
    xlabel('Angle \theta (degrees)','FontSize', 25)
    ylabel('Er Magnitude (dB)','FontSize', 25)
    hLegend = legend('Initial Beampattern','Genetic Optimised Beampattern');
    set(hLegend,'FontSize',15);
    hold off
    
    % polar plot comparing the initial beampattern with the optimised beampattern
    figure
    polar(thetaArray,abs(erInitial)','r--')
    hold on
    polar(thetaArray,abs(phasorMatrix*exciteOptimum')','k')
    grid on
    hLegend = legend('Initial Beampattern','Genetic Optimised Beampattern');
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
    disp('The number of generations performed:')
    disp(genNumber)
    disp('The change in the Main Lobe:')
    disp(optML - origML)
    disp('The reduction in the Max Side Lobe in decibels:')
    disp(((20*log10(abs(optMaxSL)/max(abs(phasorMatrix*exciteOptimum')))))+ abs((20*log10(abs(origMaxSL)/max(abs(erInitial))))))
    disp('The reducution in the Beam Ratio:')
    disp(fitnessFunctionEval - origBR)
    disp('The norm of the optimised excitation current distribution:')
    disp(norm(exciteOptimum))
    disp('-------------------------------------------------------------------')
    
    %% End of simulation
    str1 = sprintf('./GA_2D_matfiles/BR_%d_GA_ULA_gens_%d_time_%d',optBR,genNumber,maxTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'genNumber','maxTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','phasorMatrix','exciteOptimum','thetaArray','erInitial','endTime')
    
    disp('-------------------------------------------------------------------')
    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display time
    
% end
