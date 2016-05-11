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
%   for optimising a ULA with N elements
%
%   Version: 0.2 - 22/04/2015
%
%   Version specific comments:
%   1) updated commentary
%   2) vectorised where possible to increase efficiency
%   3) Add bounded currents in a later version
%-------------------------------------------------------------------------%
%--- uncomment for loop and corresponding 'end' for numerous samples ---%
% for k = 1:20 % loop for repeated output file generation

    %% Preamble
    format long % longer number format
    clear all % clears all variables
    close all % closes all open windows
    clc % clears the command window

    % Simulation start output
    disp('Simulation started ...')
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
    exciteInitial = 1.0*[1+0j (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j) (1+0j)]; % Initial current excitation vector
    BRArray = []; % array to hold changes in the beam ratio
    thetaArray = zeros(1,nAngle); % array to hold values of theta
    aArray = []; % array to hold dolph-chebyshev weights
    probArray = []; % initialise array holding probabilities of acceptance values

    seed = rand()*1e9; % generate random seed number
    rng(seed) % seed random number generator

    % for loop which generate theta and dolph-chebyshev values
    for m = 1:nAngle % cycle through angles
        theta = ((m-1)*pi)/nAngle; % angles of theta
        thetaArray(m) = theta;% storing the array elements

        %----------------- Code for producing Dolph-Chebyshev pattern comment out if not required ------------------%
        %     % checking whether to use doplh.m or dolph2.m fucntion to calculate
        %     % weights for dolph-chebysehev method
        %     if ((mod(N,2) ~= 0)&&(d <=0.49)) % if odd numebr of elements and d < 0.5
        %         [a,dph] = dolph2(d,rad2deg(theta),N,R); % dolph2.m function called, returns weights and beamwidth
        %     else
        %         [a,dph] = dolph(d,rad2deg(theta),N,R); % dolph.m function called, returns weights and beamwidth
        %     end
        %     aArray = [aArray;a]; % Array containing Dolph Chebyshev weights
        %------------------------------------------------------------------------------------------------------------%

    end

    % generate initial phasor matrix
    [er] = erGenMatULA(d,phi,N,nAngle,thetaArray);

    % subject broken array to initial excitation current vector
    erInitial = er*exciteInitial';

    % find peaks of Er for array
    [ML, maxSL, BR] = peakFinderULA(erInitial)

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

    %% Metropolis Algorithm - ULA
    % Creating temperature profile and algorithm stopping criteria
    maxIters = 500000; % maximum iterations
    maxTime = 10; % max time in seconds for time based tests
    minBR = 0.075; % min BR target for performance based tests
    T0 = 0.1; % initial max temperature
    Tfinal = 0.001; % final minimum temperature
    Tarray = [T0 zeros(1,maxIters-1)]; % array of temperature valeus
    alpha = nthroot(Tfinal,maxIters); % value of alpha

    % creating temperature values according to specified profile
    for n = 2:maxIters
        Tarray(n) = Tarray(n-1)*alpha;
    end

    tic % start timer

    for n = 1:maxIters

        gauss = Tarray(n)*(randn(1,N) + randn(1,N)*1j); % create a new gaussian random perturbation
        exciteGauss = exciteCurrent + gauss; % addition of random perturbation to currente xcitation current vector

        %----- Current Bounds: Comment Out if not required ---------------%
        exciteGauss((imag(exciteGauss) > 5)) = (real(exciteGauss((imag(exciteGauss) > 5)))+ 5j); % if imag component is greater than +5, set it to + 5
        exciteGauss((real(exciteGauss) > 5)) = (5 + imag(exciteGauss((real(exciteGauss) > 5)))*1j); % if real component is greater than +5, set it to + 5
        exciteGauss((imag(exciteGauss) < -5)) = (real(exciteGauss((imag(exciteGauss) < -5)))+ -5j); % if imag component is less than -5, set it to - 5
        exciteGauss((real(exciteGauss) < -5)) = (-5 + imag(exciteGauss((real(exciteGauss) < -5)))*1j); % if real component is less than -5, set it to - 5
        %----- Current Bounds: Comment Out if not required ---------------%

        %  add random excitation to current vector
        erMetropolis = er*exciteGauss';  % generate beampattern with random excitation

        [ML, maxSL,BR] =  peakFinderULA(erMetropolis); % measure corresponding beampattern characteristics
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
        if ((toc >= maxTime) || (optBR <= minBR)) % if current time is greater than time limit
            iters = n; % store corresponding iteration count
            endTime = toc;  % store corresponding finishing time
            break % bvreak out of loop
        end

    end

    erMetropolis = er*exciteOptimum'; % set optimised excitation current vector

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

    %     % Temperature profile of Metropolis Algorithm
    figure
    plot(Tarray,'k')
    grid on
    xlabel('Iterations (\it{m})','FontSize', 25)
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
    str1 = sprintf('./Metropolis_ULA_matfiles/unboundedPerformance_BR_%d_Metropolis_ULA_iters_%d_time_%d',optBR,iters,endTime);
    str2 = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
    strout = [str1 str2 ('.mat')];
    save(strout,'iters','endTime','seed','optBR','origBR','optMaxSL','origMaxSL','optML','origML','erMetropolis','erInitial','exciteOptimum','thetaArray','Tarray')

    beep % alert to indicate completion
    disp('Simulation ended ...')
    disp(date) % display date
    disp(datestr(now, 'HH:MM:SS')) % display
    % end of simulation

% end
