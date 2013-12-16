classdef ModelsOfManyTimeSeries
    %ModelsOfManyTimeSeries - compiles several models and their data
    % can optimize parameters for many models simulatneously
    % includes plotting routines    
    properties
        % array that will keep all the time series
        % instances of ModelWithGfp
        arrayOfModels = [];
    end
    
    methods
        
        % constructor - build an empty object
        function mwg = ModelsOfManyTimeSeries()
        end
        
        % create a new model and add
        function mwg = addTimeSeries(mwg, t, od, p, gfp)
            model = ModelWithGfp(p);
            model = model.setData(t, od, gfp);
            model = model.solveModel;            
            if isempty(mwg.arrayOfModels)
                mwg.arrayOfModels(1).model = model;
            else
                mwg.arrayOfModels(end+1).model = model;
            end            
        end
        

        function mwg = changeParameterValue(mwg, p, v, listOfModels)
            if nargin == 3
                listOfModels = 1:length(mwg.arrayOfModels);
            end
            for i = listOfModels
                if ~iscell(p)
                    mwg.arrayOfModels(i).model =...
                        mwg.arrayOfModels(i).model.changeParameterValue(p,v);
                else
                    mwg.arrayOfModels(i).model =...
                        mwg.arrayOfModels(i).model.changeParameterValues(p,v);                    
                end
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% COMPUTING ROUTINES %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % call solveModel for all models or the ones listed in the optional
        % variable listOfModels
        function mwg = solveModel(mwg, listOfModels)
            if nargin == 1
                listOfModels = 1:length(mwg.arrayOfModels);
            end
            for i = listOfModels;
                mwg.arrayOfModels(i).model =...
                    mwg.arrayOfModels(i).model.solveModel;
            end
        end
        
        % optimize a single parameter for OD using fminbnd
        % must provide upper and lower bounds
        function mwg = optimizeParameter(mwg, p, vmin, vmax)
            function err = errorAll(x)
                mwg = mwg.changeParameterValue(p, x);
                mwg = mwg.solveModel;
                err = 0;
                for i = 1:length(mwg.arrayOfModels)
                    err = err + mwg.arrayOfModels(i).model.calculateErrorOfFit;
                end
            end 
            v = fminbnd(@errorAll, vmin, vmax, optimset('Display', 'iter'));
            mwg = mwg.changeParameterValue(p, v);
            mwg = mwg.solveModel;
        end
 
        
        % optimize many parameters for OD  simulatenously using fminsearch
        % p is a cell array with all the parameters to optimize
        function mwg = optimizeParameters(mwg, p, listOfModels)
            if nargin == 2
                listOfModels = 1:length(mwg.arrayOfModels);
            end
            % initialize with the values read from the first model
            b0 = mwg.arrayOfModels(1).model.getParameterValues(p);

            function err = errorAll(b)
                mwg = mwg.changeParameterValue(p, b);
                mwg = mwg.solveModel;
                err = 0;
                for i = listOfModels
                    err = err + mwg.arrayOfModels(i).model.calculateErrorOfFit;
                end
            end 
            
            mdl = fminsearch(@errorAll, b0, optimset('Display', 'iter'));
            
            % set the values
            mwg = mwg.changeParameterValue(p, mdl);
            mwg = mwg.solveModel;
        end

        % optimize many parameters for OD simulatenously using fminsearch
        % p is a cell array with all the parameters to optimize
        function mwg = fitParameters(mwg, p, listOfModels, flagForPlot)
            if nargin == 2
                listOfModels = 1:length(mwg.arrayOfModels);
            end
            
            % initialize with the values read from the first model
            b0 = mwg.arrayOfModels(1).model.getParameterValues(p);
            [timehours, ods, ~] = mwg.getAllTimeSeries(listOfModels);
            odModel = zeros(size(ods));

            function modelResult = model(b, t)
                mwg = mwg.changeParameterValue(p, b);
                mwg = mwg.solveModel(listOfModels);
                c = 0;
                for j = listOfModels
                    c = c + 1;
                    odModel(:, c) =...
                        interp1(mwg.arrayOfModels(j).model.timehours,...
                        mwg.arrayOfModels(j).model.x, timehours, 'nearest');
                end
                modelResult = log(odModel(:));
            end 
            
            [mdl,R,J,CovB,MSE] =...
                nlinfit(timehours, log(ods(:)),...
                @model, b0, optimset('Display', 'iter'))
            % set the values
            mwg = mwg.changeParameterValue(p, mdl);
            mwg = mwg.solveModel;
            
            if (nargin > 3) && (flagForPlot)
                figure;
                plot(timehours, ods(1:length(timehours)), 'b-');
                hold on;
                [ypred, delta] = nlpredci(@model,...
                    timehours,mdl,R,'Covar',CovB,...
                    'MSE',MSE,'SimOpt','on');
                ci = nlparci(mdl,R,'covar',CovB)
                lower = ypred - delta;
                upper = ypred + delta;
                plot(timehours,exp(ypred(1:length(timehours))),'k','LineWidth',2);
                plot(timehours,...
                    [exp(lower(1:length(timehours))),...
                    exp(upper(1:length(timehours)))]...
                    ,'r--','LineWidth',1.5);
                hold off;

            end
        end
        
        % compile all data points in a single matrix to use with nlinfit
        function [timehours, ods, gfps] = getAllTimeSeries(mwg, listOfModels)
            if nargin == 1
                listOfModels = 1:length(mwg.arrayOfModels);
            end
            timehours = linspace(0, 40, 240); 
            ods  = zeros(length(timehours), length(listOfModels));
            gfps = zeros(length(timehours), length(listOfModels));
            % compile all times
            c = 0;
            for i = listOfModels
                c = c+1;
                t = mwg.arrayOfModels(i).model.timehours(:);
                o = mwg.arrayOfModels(i).model.od(:);
                g = mwg.arrayOfModels(i).model.gfp(:);
                % interpolate data
                ods(:, c) = interp1(t, o, timehours);
                gfps(:, c) = interp1(t, g, timehours);
            end
            % remove invalid points
            %invalid = or((ods<0.01), isnan(ods));
            invalid = ods<0.01;
            invalid = max(invalid, [], 2);
            timehours = timehours(~invalid);
            ods = ods(~invalid, :);
            gfps = gfps(~invalid, :);
        end
        
        % optimize many parameters for GFP rate  using fminbnd
        function mwg = optimizeParameterGfpRate(mwg, p, vmin, vmax,...
                listOfModelsToFit)
            if nargin == 4
                listOfModelsToFit = 1:length(mwg.arrayOfModels);
            end
            % function that adds up error from all models
            function err = errorAll(x)                
                mwg = mwg.changeParameterValue(p, x);
                err = 0;                
                for i = listOfModelsToFit
                    err = err +...
                        mwg.arrayOfModels(i).model.calculateErrorOfGfpRate;
                end
            end 
            
            v = fminbnd(@errorAll, vmin, vmax, optimset('Display', 'iter'));
            mwg = mwg.changeParameterValue(p, v);
        end

        % optimize many parameters for GFP simulatenously using fminsearch
        % p is a cell array with all the parameters to optimize        
        function mwg = optimizeParametersForGfpRate(mwg, p, listOfModelsToFit)
            if nargin == 2
                listOfModelsToFit = 1:length(mwg.arrayOfModels);
            end

            % initialize with the values read from the first model
            b0 = mwg.arrayOfModels(1).model.getParameterValues(p);

            function err = errorAll(b)
                mwg = mwg.changeParameterValue(p, b);
                err = 0;
                for i = listOfModelsToFit
                    err = err +...
                        mwg.arrayOfModels(i).model.calculateErrorOfGfpRate;
                end
            end 
            
            mdl = fminsearch(@errorAll, b0, optimset('Display', 'iter'));
            
            % set the values
            mwg = mwg.changeParameterValue(p, mdl);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% PLOTTING ROUTINES %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % Plot a model variable by calling plotDataAndModel for all
        % the models or a list of models (listOfModelsToPlot)
        function [] = plotAllModels(mwg, listOfModelsToPlot, variableToPlot)
            if nargin == 1
                listOfModelsToPlot = 1:length(mwg.arrayOfModels);
            end
            if nargin < 3
                variableToPlot = 'od';
            end
            nModles = length(listOfModelsToPlot);
            cmap = jet(nModles);
            for i = 1:nModles
                mwg.arrayOfModels(listOfModelsToPlot(i)).model.plotDataAndModel(...
                    cmap(i, :), variableToPlot);
                hold on;
            end
            hold off;
        end

        % Plot a state variable by calling plotModel for all
        % the models or a list of models (listOfModelsToPlot)
        function [] = plotSimulation(mwg, listOfModelsToPlot, variableToPlot)
            nModles = length(listOfModelsToPlot);
            cmap = jet(nModles);
            for i = 1:nModles
                mwg.arrayOfModels(listOfModelsToPlot(i)).model.plotModel(...
                    cmap(i, :), variableToPlot);
                hold on;
            end
            hold off;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% PRINTING ROUTINES %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % prints a list of all the parameters int he model
        function [] = printParameterSet(mwg)
            p = mwg.arrayOfModels(1).model.p;
            fnames = fieldnames(p);
            for i = 1:length(fnames)
                fprintf('p.%s = %f;\n', fnames{i}, p.(fnames{i}));
            end
                
        end
    end
    
end

