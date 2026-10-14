% Clear
clear; clc; close all;

%% User input
mode = checkInput("Select mode (1/2/3): ", 0, 0);

% Validate mode selection
while ~ismember(mode, [1 2 3])
    fprintf("Invalid mode. Please enter 1, 2, or 3.\n");
    mode = checkInput("Select mode (1/2/3): ", 0, 0);
end

% Initialise geometry and flow structures
geom = struct();  % Hold geometry variables
flow = struct();  % Hold flow properties

% Geometry inputs
geom.S = checkInput("Enter the area of the wing in m^2 (e.g. 14.86): ", 1, 0);
geom.AR = checkInput("Enter the aspect ratio (e.g. 5.63): ", 1, 0);
geom.lambda = checkInput("Enter the taper ratio (e.g. 0.95): ", 1, 0);
geom.eps_tip = checkInput("Enter the twist angle at the wingtip in degrees(e.g. -3): ", 0, 0);

% Derived geometry parameters
geom.b = sqrt(geom.AR * geom.S);
geom.c_root = 2 * geom.S / (geom.b * (1 + geom.lambda));
geom.c_tip  = geom.lambda * geom.c_root;

% Flow property inputs
flow.Vinf = checkInput("Enter the freestream velocity in m/s (e.g. 45): ", 1, 0);
flow.h    = checkInput("Enter the altitude in m (e.g. 1000): ", 1, 0);

% Airfoil data loading
while true
    filename = input("Enter the name of the airfoil datafile (e.g. N65-2-415.dat): ", "s");

    % Check file exists before loading
    if exist(filename, 'file') ~= 2
        fprintf("File not found. Please enter a valid filename.\n");
        continue;
    end

    airfoil = loadAirfoil(filename);
    break;
end

%% Mode handling

% Switch-case handles different modes selected (and invalid modes which should technically be unreachable)
switch mode
    case 1
        % Gather inputs for both N and AoA
        N = checkInput("Enter N: ", 1, 1);
        alpha_inf = deg2rad(checkInput("Enter AoA (deg): ", 0, 0));  % user types degrees
        flow.alpha_inf = alpha_inf;

        % Discretize wing
        disc = buildDiscretisation(geom, N);

        % Solve lifting line w/ horshoe theory
        liftline = solveLiftingLine(geom, flow, airfoil, disc);

        % Plot lift distribution
        figure(1);
        liftline = plotLiftDistribution(flow, liftline, disc, geom);

        % Save figure
        filename = sprintf("Case1_LiftDistr_AoA=%.1f_N=%d.png", rad2deg(flow.alpha_inf), N);
        saveas(gcf, filename);

        % Integrate for Cl & Cd values
        liftline = simpsonIntegration(disc, liftline, geom, flow);

        % Export results to csv file
        writeCase1Results(flow, N, liftline);

    case 2
        % Gather inputs for both N and AoA
        N = checkInput("Enter N: ", 1, 1);
        alpha_inf = (deg2rad(-20:1:30))';

        % Discretize wing
        disc = buildDiscretisation(geom, N);

        % Initialise Cl and Cd vectors
        Cl_wing = zeros(length(alpha_inf), 1);
        Cd_wing = zeros(length(alpha_inf), 1);

        % Loop through all angles of attack to find Cl and Cd at each
        for i = 1:length(alpha_inf)
            flow.alpha_inf = alpha_inf(i);

            % Solve lifting line problem using current N and integrate for Cl & Cd
            liftline = solveLiftingLine(geom, flow, airfoil, disc);
            liftline = simpsonIntegration(disc, liftline, geom, flow);

            Cl_wing(i) = liftline.Cl_wing;
            Cd_wing(i) = liftline.Cd_wing;
        end

        % Plot lift curve
        figure(1);
        plotCoeff(alpha_inf, Cl_wing, "C_L", 'b', N);
        set(gcf, 'Name', "C_L vs AoA");

        % Save lift curve
        filename = sprintf("Case2_CL_vs_AoA_N=%d.png", N);
        saveas(gcf, filename);
        theme('light');

        % Plot drag polar
        figure(2);
        plotCoeff(alpha_inf, Cd_wing, "C_D", 'r', N);
        set(gcf, 'Name', "C_D vs AoA");

        % Save drag polar
        filename = sprintf("Case2_CD_vs_AoA_N=%d.png", N);
        saveas(gcf, filename);
        theme('light');

        % Export results to csv file
        writeCase2Results(flow, alpha_inf, N, Cl_wing, Cd_wing);

    case 3
        % Gather input for AoA
        alpha_inf = deg2rad(checkInput("Enter AoA (deg): ", 0, 0));    % user types degrees
        flow.alpha_inf = alpha_inf;

        % Define range of discretisation values to compare convergence
        N = [10, 20, 50, 100];

        % Initialise Cl and Cd vectors
        Cl_wing_N = zeros(length(N),1);
        Cd_wing_N = zeros(length(N),1);

        % Create figure for lift distribution comparison
        figure(1); 
        
        % Loop through all N's
        for i = 1:length(N)

            % Discretize wing for current N
            disc = buildDiscretisation(geom, N(i));

            % Solve lifting line problem using current N
            liftline = solveLiftingLine(geom, flow, airfoil, disc);

            % Plot lift distribution in subplot
            subplot(2,2,i);
            liftline = plotLiftDistribution(flow, liftline, disc, geom);

            % Formatting
            title(sprintf('N = %d', N(i)));
            xlabel('Position from center (m)');
            ylabel('Lift per Unit Span (N/m)');
            grid on;

            % Integrate for Cl and Cd
            liftline = simpsonIntegration(disc, liftline, geom, flow);

            Cl_wing_N(i) = liftline.Cl_wing;
            Cd_wing_N(i) = liftline.Cd_wing;
        end

        % Add overall figure title
        sgtitle(sprintf('Lift Distribution Comparison (\\alpha = %.2f deg)', rad2deg(flow.alpha_inf)));

        % Save figure
        filename = sprintf("Case3_LiftDistr_AoA=%.1f.png", rad2deg(flow.alpha_inf));
        saveas(gcf, filename);
        theme('light');

        % Export results to csv file
        writeCase3Results(flow, N, Cl_wing_N, Cd_wing_N);

    otherwise
        % Handle invalid mode
        error("Invalid mode selected"); % Should be unreachable but for good practice
end


%% Verify Input
function val = checkInput(prompt, isPos, isInt)

    while true
        str = input(prompt, "s");   % Read as string
        val = str2double(str);      % Convert to number

        % Check if numeric
        if isempty(val) || ~isnumeric(val) || isnan(val)
            fprintf("Invalid input. Please enter a number.\n");
            continue;
        end

        % Check positivity
        if isPos && val <= 0
            fprintf("Input must be positive.\n");
            continue;
        end

        % Check integer
        if isInt && floor(val) ~= val
            fprintf("Input must be an integer.\n");
            continue;
        end

        % If all checks pass, exit loop
        break;
    end
end

%% Extract airfoil data
function airfoil = loadAirfoil(filename)

    data = readmatrix(filename);

    airfoil.alpha = deg2rad(data(:, 1));
    airfoil.Cl    = data(:,2);
    airfoil.Cd    = data(:,3);
end

%% Discretisation
function disc = buildDiscretisation(geom, N)

    % Wing span and panel width
    b  = geom.b;
    dy = b / N;a

    % Define spanwise edges, panel centres, and edge intervals
    edges = -b/2 : dy : b/2;
    yc    = edges(1:end-1) + dy/2;
    yints = [edges(1:end-1)', edges(2:end)'];

    % Store into disc struct
    disc.N     = N;
    disc.dy    = dy;
    disc.yc    = yc;
    disc.yints = yints;

    % Define chord and twist distributions
    disc.c = geom.c_root - (geom.c_root - geom.c_tip) .* abs(yc)/(b/2);
    disc.eps = deg2rad(geom.eps_tip) .* abs(yc)/(b/2);

end

%% Calculate downwash
function W = calcDownwash(disc, gamma)

    N = disc.N;

    % Initialise downwash vector
    W = zeros(1, N);
    for i = 1:N
        % Estimate downwash components at i-th panel
        w_i = hshoe(disc.yc(i), disc.yints, gamma(:));
        % Total downwash = sum of components
        W(i) = sum(w_i);
    end

end

%% Solve Lifting Line Function
function liftline = solveLiftingLine(geom, flow, airfoil, disc)

    Vinf = flow.Vinf;
    alpha_inf = flow.alpha_inf;

    % Initial elliptical guess for circulation
    gamma_elliptical = 4 * Vinf * geom.S * alpha_inf / geom.b * sqrt(1 - (disc.yc/(geom.b/2)).^2); 
    
    % Set convergence parameters
    maxIter = 1000;
    iter = 0;
    err = Inf;
    tol = 1e-5;
    D = 0.05;
    gamma = gamma_elliptical;
    
    % Loop until convergence or maxIter reached
    while (iter < maxIter) & (err > tol)
    
        iter = iter + 1;
    
        % Compute downwash
        W = calcDownwash(disc, gamma);
    
        % Effective angle of attack
        alpha_e = disc.eps + atan((Vinf * sin(alpha_inf) - W) ...
                                ./(Vinf * cos(alpha_inf)));
    
        % Interpolate Cl and Cd
        Cl = interp1(airfoil.alpha, airfoil.Cl, alpha_e, "linear");
        Cd = interp1(airfoil.alpha, airfoil.Cd, alpha_e, "linear");
    
        % Update circulation (incl. damping factor) 
        gamma_new = 0.5 * flow.Vinf .* disc.c .* Cl;
        gamma_next = gamma + D * (gamma_new - gamma);
    
        % Quantify error
        err = max(abs(gamma_next - gamma));
    
        % Move to next iteration
        gamma = gamma_next;
    end

    % Create liftline structure w/ relevant variables
    liftline.gamma = gamma;
    liftline.Cl    = Cl;
    liftline.Cd    = Cd;
    liftline.W     = W;
end

%% Plot Lift Distribution
function liftline = plotLiftDistribution(flow, liftline, disc, geom)

    % Gather density at current height
    [~,~,~,rho] = atmosisa(flow.h);

    % Calculate lift per unit span
    L = rho * flow.Vinf * liftline.gamma;
    
    % Prepare for plotting by creating wing disp vector
    y_plot = [-geom.b/2, disc.yc, geom.b/2];

    % Enforce endpoint lift = 0 as per theory
    L_plot = [0, L, 0];

    % Plot lift distrubuton
    plot(y_plot, L_plot, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', disc.N));

    % Formatting
    title(sprintf('Wing Lift Distribution (\\alpha = %.1f deg)', rad2deg(flow.alpha_inf)));
    xlabel('Spanwise position y [m]');
    ylabel('Lift per unit span L'' [N/m]');
    theme('light');

    % Update liftline structure
    liftline.L = L;
end

%% Plot Cd or Cl vs AoA
function plotCoeff(alpha_inf, coeff, coeffLabel, colour, N)

    % Plot drag polar or lift curve
    plot(rad2deg(alpha_inf), coeff, [colour '-'], 'LineWidth', 1.5);
    
    % Formatting (customised for Cl or Cd)
    xlabel('\alpha_\infty (deg)');
    ylabel(coeffLabel);
    if coeffLabel == "C_L"
        title(sprintf('Lift Curve (N = %d)', N));
    else
        title(sprintf('Drag Polar (N = %d)', N));
    end

    grid on;
end

%% Simpson Integration
function liftline = simpsonIntegration(disc, liftline, geom, flow)

    % Extract variables from structs for ease of use
    N = disc.N;
    dy = disc.dy;
    Cl = liftline.Cl;
    Cd = liftline.Cd;
    
    % Calculate incidence angle w/ downwash
    alpha_i = atan(liftline.W / flow.Vinf);
    
    % Integrand setup w/ edge case handling
    f = disc.c .* Cl .* cos(alpha_i); % Cl
    f = [0, f, 0];
    g = disc.c .* (Cd + Cl .* sin(alpha_i)); % Cd
    g = [0, g, 0];
    
    % Initialise Cl and Cd
    Cl_wing = 0;
    Cd_wing = 0;

    % Simpson integration over span
    for i = 2:N-1
        Cl_wing = Cl_wing + dy/6 * (f(i+1) + 4*f(i) + f(i-1));
        Cd_wing = Cd_wing + dy/6 * (g(i+1) + 4*g(i) + g(i-1));
    end

    % Normalise by S
    Cl_wing = Cl_wing/geom.S;
    Cd_wing = Cd_wing/geom.S;

    % Update liftline struct
    liftline.Cl_wing = Cl_wing;
    liftline.Cd_wing = Cd_wing;
end

%% Writing results for each case into .csv file
function writeCase1Results(flow, N, liftline)

    % Extract AoA and define output filename
    filename_out = sprintf("Results_aoa=%.1f_N=%d.txt", rad2deg(flow.alpha_inf), N);
    fid = fopen(filename_out, 'w');

    % Write header, case description, and flight condition data
    fprintf(fid, "3D Wing Aerodynamic Results\n");
    fprintf(fid, "Case 1: Single AoA\n");
    fprintf(fid, "Altitude (m): %.2f\n", flow.h);
    fprintf(fid, "Freestream Velocity (m/s): %.2f\n", flow.Vinf);
    fprintf(fid, "Angle of Attack (deg): %.2f\n", rad2deg(flow.alpha_inf));
    fprintf(fid, "Discretisation (N): %d\n\n", N);

    % Write Cl and Cd results
    fprintf(fid, "Lift Coefficient (Cl_wing): %.6f\n", liftline.Cl_wing);
    fprintf(fid, "Drag Coefficient (Cd_wing): %.6f\n", liftline.Cd_wing);

    % Close file
    fclose(fid);
end

function writeCase2Results(flow, alpha_inf, N, Cl_wing, Cd_wing)

    % Define output filename
    filename_out = sprintf("Results_N=%d.txt", N);
    fid = fopen(filename_out, 'w');

    % Write header, case description, and flight condition data
    fprintf(fid, "3D Wing Results\n");
    fprintf(fid, "Case 2\n");
    fprintf(fid, "Altitude (m): %.2f\n", flow.h);
    fprintf(fid, "Freestream Velocity (m/s): %.2f\n", flow.Vinf);
    fprintf(fid, "Discretisation (N): %d\n\n", N);

    % Write Cl and Cd results for each AoA
    fprintf(fid, "%10s %15s %15s\n", "AoA (deg)", "Cl_wing", "Cd_wing");
    for i = 1:length(alpha_inf)
        fprintf(fid, "%10.2f %15.6f %15.6f\n", ...
            rad2deg(alpha_inf(i)), Cl_wing(i), Cd_wing(i));
    end

    % Close file
    fclose(fid);
end

function writeCase3Results(flow, N, Cl_wing, Cd_wing)

    % Define output filename
    filename_out = sprintf("Results_aoa=%.1f.txt", rad2deg(flow.alpha_inf));
    fid = fopen(filename_out, 'w');

    % Write header, case description, and flight condition data
    fprintf(fid, "3D Wing Results\n");
    fprintf(fid, "Case 3\n");
    fprintf(fid, "Altitude (m): %.2f\n", flow.h);
    fprintf(fid, "Freestream Velocity (m/s): %.2f\n", flow.Vinf);
    fprintf(fid, "Angle of Attack (deg): %.2f\n\n", rad2deg(flow.alpha_inf));

    % Write Cl and Cd results for each N
    fprintf(fid, "%10s %15s %15s\n", "N", "Cl_wing", "Cd_wing");
    for i = 1:length(N)
        fprintf(fid, "%10d %15.6f %15.6f\n", ...
            N(i), Cl_wing(i), Cd_wing(i));
    end

    % Close file
    fclose(fid);
end