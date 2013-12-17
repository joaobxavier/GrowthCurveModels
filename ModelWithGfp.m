classdef ModelWithGfp
    % ModelWithGfp - keeps data and model associated with it 
    
    properties (SetAccess = protected)
        p;
%         y0;
        % experimental data (for fitting)
        timehours;
        od;
        gfp;
        % state variables
        x; % Biomass
        c; % Carbon
        n; % Nitrogen
        ie; % Iron
        ni; % Internal nitrogen
        ii; % Internal iron
        gi; % Internal gfp
        % properties for experimental replicates
        hasReplicates = false; % value will be true if there are replicates
        odReplicates;
        gfpReplicates;    
    end
    
    methods
                
        % constructor: creates a model
        % p - model parameters
        function mwg = ModelWithGfp(p)
            mwg.p = p;
        end
        
        % set the array od initial values for simulation
        function mwg = setInitialValues(mwg, y0)
            mwg.x0 = y0(1);
            mwg.C0Value = y0(2);
            mwg.N0Value = y0(3);
            mwg.I0Value = y0(4);
            mwg.ni0 = y0(5);
            mwg.ii0 = y0(6);
            mwg.gi0 = y0(7);
        end
        
        % Add an experimental data time series
        % od and gfp may contain mutliple replicates
        function mwg = setData(mwg, timehours, od, gfp)
            mwg.timehours = timehours;
            if size(od, 2) > 1
                mwg.od = median(od, 2);
                mwg.gfp = median(gfp, 2);
                mwg.hasReplicates = true;
                mwg.odReplicates = od;
                mwg.gfpReplicates = gfp;
            else
                mwg.od = od;
                mwg.gfp = gfp;
            end
        end
        
%         % Change value of a parameter
%         % p - parameter name (field in mwg.p)
%         % v - new value for parameter
%         function mwg = changeParameterValue(mwg, p, v)
%             mwg.p.(p) = v;
%         end
        
        % Change value of many parameters
        % p - cell array with parameter names (field in mwg.p)
        % v - araray with new value for parameters
        function mwg = changeParameterValues(mwg, p, v)
            if ~iscell(p)
                p = {p};
            end
            
            for i = 1:length(p)
                v(i)
                p{i}
                mwg.p.(p{i}) = v(i);
%                 mwg = setfield(mwg, p{i}, v(i));
            end
        end
        
        % Get the value of parameters
        % p - cell array with parameter names (field in mwg.p)
        function v = getParameterValues(mwg, p)
            if ~iscell(p)
                p = {p};
            end
            v = zeros(1, length(p));
            for i = 1:length(p)
                v(i) = mwg.p.(p{i});
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% MATHEMATICAL MODEL %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Implements the biochemical model
        function dy = growthFunctionInternalNRWithGFP(mwg, t, y)
            
            %Growth with Carbon, Nitrogen snd iron. No internal carbon pool
            x  = y(1); % Biomass
            c  = y(2); % Carbon
            n  = y(3); % Nitrogen
            ie = y(4); % Iron
            ni = y(5); % Internal nitrogen
            ii = y(6); % Internal iron
            gi = y(7); % Internal GFP
            
            % if there is carbon available, grow
            % if there is no carbon stop growing
            growth = mwg.growthRate(c, ni, ii);
            
            % biomass
            dx = (growth - mwg.p.kd * (c == 0)) * x;
            
            % carbon
            dc = - 1/mwg.p.Yc * growth * x;
            
            % nitrogen
            dn = - 1/mwg.p.Yn * growth * x * (n > 0) * (ni > 0);
            
            % iron
            di = - 1/mwg.p.Yi * growth * x * (ie > 0);
            
            % internal nitrogen
            dni = - (1/mwg.p.Yn + ni) * growth * (ni > 0) * (n == 0);
            
            % internal iron (subtract 1 because of dilution during growth)
            dii = - (1/mwg.p.Yi + ii) * growth * (ii > 0) * (ie == 0);
            
            % gfp production is implemented as a separate function
            % to facilitate fitting the gfp parameters
            %d = mwg.dgfpi(c, ie, ii, x, gi, n) - dx/x * gi;

            d = mwg.dgfpi(growth, x, gi) - dx/x * gi;
            
            dy = [dx; dc; dn; di; dni; dii; d];
        end
        
        
        % the growth rate function
        function g = growthRate(mwg, c, ni, ii)
            g = mwg.p.qG .* ii .* ni .* (c > 0);            
        end
        
        
        % a piece-wise model for the gfp expression rate
        function d = dgfpi(mwg, growth, x, gi)
            % calculate the fraction of the growth rate
            g = growth;
            f = g ./ (mwg.p.qG .* mwg.p.ii0.* mwg.p.ni0) + eps;
            % density dependent process
            d = mwg.p.qRb .* x; 
            % up regulation due to mild growth limitation
            d = d + mwg.p.qR .* (1./f -1) ;
            % down regulation due to severe nutrient limitation
            d = d .*  1 ./ (1 + (mwg.p.kg./f).^mwg.p.h1) ;
            % decay due to growth arrest
            d = d - mwg.p.kdGFP .* (growth <= 0) .* gi;
        end

        % solve the ordinary differential equations
        function mwg = solveModel(mwg, timehours, y0)
            if nargin < 3
                % if y0 is not provided, use predefined y0
                y0 = [mwg.p.x0, mwg.p.C0Value, mwg.p.N0Value, ...
                    mwg.p.I0Value, mwg.p.ni0, mwg.p.ii0, mwg.p.gi0];
            end
            if nargin < 2
                % if timehours is not provided, use predefined timehours
                timehours = mwg.timehours;
            end
            % set the value of the qG
            mwg.p.qG    = mwg.p.muMax / mwg.p.ni0 / mwg.p.ii0;
%             y0(5) = mwg.p.ni0;
%             y0(6) = mwg.p.ii0;
            %
            f = @(t,y)mwg.growthFunctionInternalNRWithGFP(t, y);
            options = odeset('NonNegative', 1:length(y0));
            tFinal = max(timehours);
            [t, y] = ode45(f, [0 tFinal], y0, options);
            % interpolate so that simulation has same time points as data
            y = interp1(t, y, timehours);
            t = timehours;
            % copy the state variables
            mwg.x  = y(:, 1); % Biomass
            mwg.c  = y(:, 2); % Carbon
            mwg.n  = y(:, 3); % Nitrogen
            mwg.ie = y(:, 4); % Iron
            mwg.ni = y(:, 5); % Internal nitrogen
            mwg.ii = y(:, 6); % Internal iron
            mwg.gi = y(:, 7); % Internal gfp
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% PARAMETER FITTING %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % compares OD from data with OD from simulation. Calculates the
        % sum of squared errors.
        % Ignores OD < 0.01 and NaNs
        function err = calculateErrorOfFit(mwg)
            % find valid times points (where od is at least 0.01)
            v = find(and(mwg.od >= 0.01, ~isnan(mwg.od)));
            % compute error
            err = sum((log(mwg.x(v)./mwg.od(v))).^2);
        end
        
        % compares model with data for GFP. The comparison is done
        % by calculating the rate of Gfp/OD. Present version uses
        % the rate calculated directly from data and compares with
        % a model where expression is only a function of the growth rate.
        % Ignores points with OD < 0.01 and NaNs
        function err = calculateErrorOfGfpRate(mwg)
            % find valid times points (where od is at least 0.01)
            gfpPerCell = mwg.gfp./mwg.od;
            growthData = ModelWithGfp.calculateRate(...
                mwg.timehours, log(mwg.od));
            data = ModelWithGfp.calculateRate(...
                mwg.timehours,gfpPerCell)...
                + gfpPerCell.*growthData;
            v = find(and(mwg.od >= 0.01, ~isnan(data)));
            % THESE TWO LINES CAN BE USED TO COMPARE DIRECTLY WITH
            % MODEL
            %  growthModel = mwg.growthRate(mwg.c, mwg.ni, mwg.ii);
            %  model = mwg.dgfpi(growthModel,...
            %    mwg.x,  mwg.gi);
            model = mwg.dgfpi(growthData,...
                mwg.od,  gfpPerCell);   
            % compute error
            err = sum((data(v) - model(v)).^2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% PLOTTING ROUTINES %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Plot the data and model. 
        % colorLine - RGB color value
        %Accepts following options for
        % variableToPlot: 'od', 'gfp', 'gi', 'dgidt' and 'dlogxdt'
        function [] = plotDataAndModel(mwg, colorLine, variableToPlot)
            switch variableToPlot
                case 'od'
                    dataToPlot = mwg.od;
                    modelToPlot = mwg.x;
                    if mwg.hasReplicates
                        dataToPlotTop = max(mwg.odReplicates, [], 2);
                        dataToPlotBottom = min(mwg.odReplicates, [], 2);
                    end
                case 'gfp'
                    dataToPlot = mwg.gfp;
                    modelToPlot = mwg.gi.*mwg.x;
                    if mwg.hasReplicates
                        dataToPlotTop = max(mwg.gfpReplicates, [], 2);
                        dataToPlotBottom = min(mwg.gfpReplicates, [], 2);
                    end
                case 'gi'
                    dataToPlot = mwg.gfp./mwg.od;
                    modelToPlot = mwg.gi;
                    if mwg.hasReplicates
                        dataToPlot = median(...
                            mwg.gfpReplicates./mwg.odReplicates, 2);
                        dataToPlotTop = max(...
                            mwg.gfpReplicates./mwg.odReplicates, [], 2);
                        dataToPlotBottom = min(...
                            mwg.gfpReplicates./mwg.odReplicates, [], 2);
                    end                    
                case 'dgidt'
                    gfpPerCell = mwg.gfp./mwg.od;
                    growthData = ModelWithGfp.calculateRate(...
                        mwg.timehours, log(mwg.od));
                    dataToPlot = ModelWithGfp.calculateRate(...
                        mwg.timehours,gfpPerCell)...
                        + gfpPerCell.*growthData;
                     growthModel = mwg.growthRate(mwg.c, mwg.ni, mwg.ii);
                     modelToPlot = mwg.dgfpi(growthModel, mwg.x,  mwg.gi);
                     % use this line to plot the rate calculated directly
                     % from data growth rate
                     % modelToPlot = mwg.dgfpi(growthData, mwg.od,  gfpPerCell);
                    if mwg.hasReplicates
                        warning('plotting dgidt does not work with replicates yet')
                        dataToPlotTop = dataToPlot;
                        dataToPlotBottom = dataToPlot;
                    end                    
                    
                case 'dlogxdt'
                    dataToPlot = ModelWithGfp.calculateRate(...
                        mwg.timehours, log(mwg.od));
                    modelToPlot = mwg.growthRate(mwg.c, mwg.ni, mwg.ii);
                    if mwg.hasReplicates
                        warning('plotting dlogxdt does not work with replicates yet')
                        dataToPlotTop = dataToPlot;
                        dataToPlotBottom = dataToPlot;
                    end                    
            end
            % plotting routine
            v = find(and(mwg.od >= 0.01, ~isnan(mwg.od)));
            % plot experimental replicates
            if mwg.hasReplicates
                % use following code to plot replicate range as lines
                % instead of area
                % h = plot(mwg.timehours(v),...
                %  dataToPlotTop(v), 'LineWidth', 1, 'Color', colorLine);
                % hold on;
                % h = plot(mwg.timehours(v),...
                %  dataToPlotBottom(v), 'LineWidth', 1, 'Color', colorLine);
                h = area(mwg.timehours(v), [dataToPlotBottom(v),...
                    dataToPlotTop(v) - dataToPlotBottom(v)], 'EdgeColor', 'None');
                set(h(1), 'FaceColor', 'none');
                set(h(2), 'FaceColor', colorLine);
                hold on;
            end
            % plot median
            h = plot(mwg.timehours(v), dataToPlot(v), 'LineWidth', 2, 'Color', colorLine);
            hold on;
            % plot model
            h = plot(mwg.timehours, modelToPlot, 'LineWidth', 1, 'Color', colorLine);
            hold off;
            xlabel('Time [h]');
            ylabel(variableToPlot);
            set(gca, 'YScale', 'log','YLim', [0.01 1], 'XLim', [0 48]);
            set(gca, 'Color', [0.2 0.2 0.2]);
        end
        
        % plot a time series of model, without the data.
        % accepts any of the state variables as variableToPlot:
        % x; % Biomass
        % c; % Carbon
        % n; % Nitrogen
        % ie; % Iron
        % ni; % Internal nitrogen
        % ii; % Internal iron
        % gi; % Internal gfp
        function [] = plotModel(mwg, colorLine, variableToPlot)
            h = plot(mwg.timehours, mwg.(variableToPlot), 'LineWidth', 2, 'Color', colorLine);
            xlabel('Time [h]');
            ylabel(variableToPlot);
            set(gca, 'XLim', [0 48]);
        end
    end
    
    methods(Static)
        % numeric calculation of time derivative. Uses local averaging to
        % smooth data before calculating the finite difference.
        function r = calculateRate(time, data, windowSize)
            % get the number of replicates
            nReps = size(data, 2);
            
            % filter data
            if nargin == 2
                data = ModelWithGfp.filterData(data);
            else
                data = ModelWithGfp.filterData(data, windowSize);
            end            
            % relicate time array
            timeDifferential = diff(time);
            timeDifferential = repmat(timeDifferential(:), 1, nReps);
            % filter data to reduce noice
            r = data;
            % calculate derivative
            r = diff(r, 1, 1)./timeDifferential;
            % replicate the first data point to make data same size
            % as original data
            r = [r(1, :); r];
        end
        
        % filter data to reduce noise. Optional window size (default: 5)
        function dFiltered = filterData(d, windowSize)
            if nargin == 1
                windowSize = 5;
            end
            dFiltered = filter(ones(1,windowSize)/windowSize,1,d, [], 1);
        end
    end
end
