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
%   for optimising a ULA with N elements, where one element is inactive
%
%   Version: 0.1 - 31/03/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%-------------------------------------------------------------------------%

%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for k = 1:20 % loop for repeated output file generation
    
    format long % longer cmd window format
    clear all % clears all initialised variables
    close all % closes all open windows
    clc % clears the command window
    
    % Simulation start output
    disp('Simulation started - ULA With an Inactive Element Beam Pattern Synthesis Using a Genetic Algorithm')
    disp(date)
    disp(datestr(now, 'HH:MM:SS'))
    disp('-------------------------------------------------------------------')
    
    %% Declaration of array paramaters
    d = 0.5;  % spacing (in terms of wavelength)
    phi = 0; % progressive phase difference in radians
    N = 12;  % number of elements in the array
    nAngle = 400; % number of angles between 0 and pi
    
    % Dolph Chebyshev array parameters
    R = 20; % relative distance between sidelobe level and mainbeam in dB, (e.g., R = 30)
    thetaArray = zeros(1,nAngle); % array to hold values of theta
    aArray = []; % array to hold dolph-chebyshev weights
    
    seed = rand()*1e9; % generate random seed number
    rng(seed) % seed random number generator
    
    % Initial excitation vector
    excite = 1.0*[1 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0];
    exciteInitial = [excite(1) + excite(2)*1j ...
        excite(3) + excite(4)*1j ...
        excite(5) + excite(6)*1j ...
        excite(7) + excite(8)*1j ...
        excite(9) + excite(10)*1j ...
        excite(11) + excite(12)*1j ...
        excite(13) + excite(14)*1j ...
        excite(15) + excite(16)*1j ...
        excite(17) + excite(18)*1j ...
        excite(19) + excite(20)*1j ...
        excite(21) + excite(22)*1j ...
        excite(23) + excite(24)*1j ...
        ];
    
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
    
    % generate initial beampattern for array
    [er] = erGenMatULAInactiveElement(d,phi,N,nAngle,thetaArray);
    
    % subject broken array to initial excitation current vector
    erInitial = er*exciteInitial';
    
    % find peaks of Er for array
    [ML, maxSL, BR] = peakFinderULAInactiveElement(erInitial);
    
    % store original beampattern parameters
    origBR = BR; % store original beam ratio
    origMaxSL = maxSL; % store original max sidelobe
    origML = ML; % store original main lobe
    
    %% Genetic Algorithm
    % algorithm parameters and stopping criteria
    genNumber = inf; % number of generations
    maxTime = 20; % max time limit
    minBR = 0.18; % minimum BR for performance comparison
    
    % initial excitations and bounds on excitations
    exciteInitial = 1*[1 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0]; % Initial excitations
    ub = 5*[1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1]; % upper bound of excitations
    lb = -5*[1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1]; % lower bound of excitations
    
    % initialise fitness (cost) function
    FitnessFunction = @(exciteInitial) fitness_function_GA_ULA_Inactive_Element(exciteInitial,er);
    
    % GA options set - N.B: Add the following for dynamic plots 'PlotFcns',{@gaplotbestf,@gaplotstopping}
    options = gaoptimset('StallGenLimit',inf,'Generations',genNumber,'FitnessLimit',minBR,'TimeLimit',maxTime)
    
    tic % start timer
    % calling the GA algorithm
    [exciteOptimum,fitnessFunctionEval,exitFlag,output] = ga(FitnessFunction,24,[],[],[],[],lb,ub,[],options)
    endTime = toc; % stop timer
    
    % print to screen options linked to simulannealbnd
    fprintf('The number of generations was : %d\n', output.generations);
    fprintf('The number of function evaluations was : %d\n', output.funccount);
    fprintf('The best function value found was : %g\n', fitnessFunctionEval);
    
    % store number of generations completed
    genNumber = output.generations;
    
    % recombining initial excitation vector
    exciteInitial = [(exciteInitial(1) + exciteInitial(2)*1j)...
        (exciteInitial(3) + exciteInitial(4)*1j)...
        (exciteInitial(5) + exciteInitial(6)*1j)...
        (exciteInitial(7) + exciteInitial(8)*1j)...
        (exciteInitial(9) + exciteInitial(10)*1j)...
        (exciteInitial(11) + exciteInitial(12)*1j) ...
        (exciteInitial(13) + exciteInitial(14)*1j)...
        (exciteInitial(15) + exciteInitial(16)*1j)...
        (exciteInitial(17) + exciteInitial(18)*1j) ...
        (exciteInitial(19) + exciteInitial(20)*1j)...
        (exciteInitial(21) + exciteInitial(22)*1j)...
        (exciteInitial(23) + exciteInitial(24)*1j)];
    
    % recombining optimum exciattion vector
    exciteOptimum = [(exciteOptimum(1) + exciteOptimum(2)*1j) (exciteOptimum(3) + exciteOptimum(4)*1j) (exciteOptimum(5) + exciteOptimum(6)*1j) ...
        (exciteOptimum(7) + exciteOptimum(8)*1j) (exciteOptimum(9) + exciteOptimum(10)*1j) (exciteOptimum(11) + exciteOptimum(12)*1j) ...
        (exciteOptimum(13) + exciteOptimum(14)*1j) (exciteOptimum(15) + exciteOptimum(16)*1j) (exciteOptimum(17) + exciteOptimum(18)*1j) ...
        (exciteOptimum(19) + exciteOptimum(20)*1j) (exciteOptimum(21) + exciteOptimum(22)*1j) (exciteOptimum(23) + exciteOptimum(24)*1j)];
    
    
    % find peaks of Er for array
    [ML, maxSL, BR] = peakFinderULAInactiveElement(er*exciteOptimum');
    
    % store optimised beampattern parameters
    optBR = BR;
    optMaxSL = maxSL;
    optML = ML;
    
    %% Genetic Algorithm ULA Simulation Plots
    % linear plot comparing the initial beampattern with the optimised beampattern
    figure
    plot(rad2deg(thetaArray),abs(erInitial),'r--')
    hold on
    plot(rad2deg(thetaArray),abs(er*exciteOptimum'),'k')
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
    plot(rad2deg(thetaArray),(20*log10(abs(er*exciteOptimum')/max(abs(er*exciteOptimum')))),'k')
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
    polar(thetaArray,abs(er*exciteOptimum')','k')
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
    disp(((20*log10(abs(optMaxSL)/max(abs(er*exciteOptimum')))))+ abs((20*log10(abs(origMaxSL)/max(abs(erInitial))))))
    disp('The reducution in the Beam Ratio:')
    disp(fitnessFunctionEval - origBR)
    disp('The norm of the optimised excitation current distribution:')
    disp(norm(exciteOptimum))
    disp('-------------------------------------------------------------------')
    
    %% End of Simulation
    str1 = sprintf('./GA_ULA_Inactive_Elements_matfiles/BR_%d_GA_ULA_Inactive_Element_gens_%d_time_%d',optBR,genNumber,maxTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'genNumber','maxTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','erInitial','exciteOptimum','thetaArray','er','endTime')
    
    disp('-------------------------------------------------------------------')
    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display time
    
% end
