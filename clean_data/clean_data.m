function [xClean, yClean] = clean_data(xRaw, yRaw, varargin)
% CLEAN_DATA    Clean one-dimensional data by removing points
%   [xClean, yClean] = CLEAN_DATA(xRaw, yRaw) opens a GUI for
%   cleaning data
%
%   [xClean, yClean] = CLEAN_DATA(xRaw, yRaw, 'Interactive', false, 'OutputPath', 'mypath.hdf5') cleans
%   data without using GUI and outputs the result to the provided path
%
%   [xClean, yClean] = CLEAN_DATA(xRaw, yRaw, options) uses user-provided
%   options struct as default cleaning parameters. 
% 
%   Available options (defaults highlighted) are:
%   
%   headers:
%       headers for variables (*{'X', 'Y'}*)
%   sort:
%       sort the data (*true*/false)
%   removeOutliersX:
%       remove x outliers (*true*/false)
%   outlierMethodX:
%       which x outlier method to use (*'mean'*/'median')
%   outlierThresholdX:
%       how many stdevs away to remove pts (*2*)
%   removeOutliersY:
%       remove y outliers (*true*/false)
%   outlierMethodY:
%       which y outlier method to use ('movmean'/*'movmedian'*)
%   outlierThresholdY:
%       how many stdevs away to remove pts (*2*)
%   indicesToRemove:
%       indices of points to always remove (*[]*)
%   indicesToKeep:
%       indices of points to always keep (*[]*)
%   cutoffX:
%       matrix containing rows of [xmin, xmax], i.e. bands in which
%       points should be removed (*[]*)
%
%   EXAMPLE: To open an interactive data-cleaning GUI with the 
%   x-axis labelled 'Voltage' and the y-axis labelled 'Current', and
%   no x outlier removal by default, do the following:
%
%   CLEAN_DATA(xRaw, yRaw, struct(...
%       'headers', P{'Voltage, 'Current'}}, ...
%       'removeOutliersX', false
%   ))

    %===================================================================
    %                   Definition of default options
    %===================================================================

    % Default data headers
    defaultOptions.headers = {'X', 'Y'};

    % By default, we sort the data
    defaultOptions.sort = true;

    % By default, we remove NaN and missing values
    defaultOptions.removeMissing = true;
    
    % By default, we do not remove X outliers.
    defaultOptions.removeOutliersX = true;
    defaultOptions.outlierMethodX = "mean";
    defaultOptions.outlierThresholdX = 2;
    defaultOptions.outlierWindowX = 20;

    % By default, we do remove Y outliers, using a moving median detection
    % window with a threshold of 1 and a window size of 5.
    defaultOptions.removeOutliersY = true;
    defaultOptions.outlierMethodY = "movmedian";
    defaultOptions.outlierThresholdY = 2;
    defaultOptions.outlierWindowY = 20;

    % Manual indices to remove and keep. By default, these are empty.
    defaultOptions.indicesToRemove = [];
    defaultOptions.indicesToKeep = [];

    % x bands in which to remove data. By default, this is empty.
    defaultOptions.cutoffX = [];

    %===================================================================
    %                   Argument parsing
    %===================================================================
    
    p = inputParser;
    p.StructExpand = false;

    % User-specified options struct as third positional parameter
    addOptional(p, 'Options', defaultOptions)

    % Variable names
    validString = @(a) ischar(a) || isstring(a);

    % Whether the data cleaning should be interactive or not
    addParameter(p, 'Interactive', true, @islogical)

    % File output path
    addParameter(p, 'OutputPath', 'data.hdf5', validString)

    parse(p, varargin{:})

    interactive = p.Results.Interactive;
    outputPath = p.Results.OutputPath;
    options = defaultOptions;
    userOptions = p.Results.Options;
    optionFields = fields(userOptions);
    for i = 1:length(optionFields)
        field = optionFields{i};
        val = userOptions.(field);
        options.(field) = val;
    end

    %===================================================================
    %       Initial data cleaning and sorting using default options
    %===================================================================

    % Make raw data matrix
    rawData = [xRaw yRaw];
    
    % Perform default data cleaning
    [xClean, yClean, keptInds, removedInds] = clean(xRaw, yRaw, options);
    cleanData = [xClean, yClean];

    % Get sorted data for display
    [xSorted, ySorted] = sortData(xRaw, yRaw);

    % Exit program if not in interactive mode
    if ~interactive
        exportData(struct(), struct());
        return
    end

    %===================================================================
    %               GUI definition (for interactive mode)
    %===================================================================
    
    % Create the GUI
    fig = uifigure;
    fig.Name = "PEPL Data Cleaner v0.1.0";

    % Create the main layout, with a single row and two columns.
    % The first column holds the plots, and the second column holds the
    % options.
    gl = uigridlayout(fig, [1, 2]);
    gl.RowHeight = {'1x'};
    gl.ColumnWidth = {'1x', 'fit'}; 

    % Grid layout for axes. Two rows, one column. First row shows unaltered
    % (raw) data. Second column has cleaned data.
    gl_axes = uigridlayout(gl, [2, 1]);
    gl_axes.Layout.Row = 1;
    gl_axes.Layout.Column = 1;

    ax1 = uiaxes(gl_axes);
    ax1.Layout.Row = 1;
    ax1.Layout.Column = 1;
    
    ax2 = uiaxes(gl_axes);
    ax2.Layout.Row = 2;
    ax2.Layout.Column = 1;

    % Grid layout for options. Two columns.
    % First column holds description, second column holds interactive
    % element.
    numOptions = 8;
    gl_opts = uigridlayout(gl, [numOptions+2, 2]);
    gl_opts.Layout.Row = 1;
    gl_opts.Layout.Column = 2;
    gl_opts.ColumnWidth = {'fit', '1x'};
    rowHeights = cell(numOptions+2, 1);

    defaultHeight = 25.0;
    for i = 1:numOptions+2
        rowHeights{i} = defaultHeight;
    end

    currentRow = 1;
    
    % First row has checkbox for sorting data.
    lbl_sort = uilabel(gl_opts);
    lbl_sort.Layout.Row = currentRow;
    lbl_sort.Layout.Column = 1;
    lbl_sort.Text = "Sort data";
    cbx_sort = uicheckbox(gl_opts, 'Text', '');
    cbx_sort.Layout.Row = currentRow;
    cbx_sort.Layout.Column = 2;
    cbx_sort.Value = true;
    cbx_sort.ValueChangedFcn = @sortDataCallback;
    currentRow = currentRow + 1;

    % Checkbox for removing X outliers
    lbl_outlierX = uilabel(gl_opts, 'Text', "Remove X Outliers");
    lbl_outlierX.Layout.Row = currentRow;
    lbl_outlierX.Layout.Column = 1; 
    cbx_outlierX = uicheckbox(gl_opts, 'Text', '', 'Value', options.removeOutliersX);
    cbx_outlierX.Layout.Row = currentRow;
    cbx_outlierX.Layout.Column = 2;
    cbx_outlierX.ValueChangedFcn = @outlierXCheckboxCallback;
    currentRow = currentRow + 1;

    % Numeric entry for X outlier threshold
    lbl_outlierThresholdX = uilabel(gl_opts, 'Text', "X outlier threshold");
    lbl_outlierThresholdX.Layout.Row = currentRow;
    lbl_outlierThresholdX.Layout.Column = 1; 
    field_outlierThresholdX = uislider(gl_opts, 'Value', options.outlierThresholdX);
    field_outlierThresholdX.Limits = [0, 5];
    field_outlierThresholdX.Layout.Row = currentRow;
    field_outlierThresholdX.Layout.Column = 2;
    field_outlierThresholdX.ValueChangedFcn = @outlierThresholdXCallback;
    rowHeights{currentRow} = 40.0;
    currentRow = currentRow + 1;

    % Checkbox for removing Y outliers
    lbl_outlierY = uilabel(gl_opts, 'Text', "Remove Y Outliers");
    lbl_outlierY.Layout.Row = currentRow;
    lbl_outlierY.Layout.Column = 1; 
    cbx_outlierY = uicheckbox(gl_opts, 'Text', '', 'Value', options.removeOutliersY);
    cbx_outlierY.Layout.Row = currentRow;
    cbx_outlierY.Layout.Column = 2;
    cbx_outlierY.ValueChangedFcn = @outlierYCheckboxCallback;
    currentRow = currentRow + 1;

    % Numeric entry for Y outlier threshold
    lbl_outlierThresholdY = uilabel(gl_opts, 'Text', "Y outlier threshold");
    lbl_outlierThresholdY.Layout.Row = currentRow;
    lbl_outlierThresholdY.Layout.Column = 1; 
    field_outlierThresholdY = uislider(gl_opts, 'Value', options.outlierThresholdY);
    field_outlierThresholdY.Limits = [0, 5];
    field_outlierThresholdY.Layout.Row = currentRow;
    field_outlierThresholdY.Layout.Column = 2;
    field_outlierThresholdY.ValueChangedFcn = @outlierThresholdYCallback;
    rowHeights{currentRow} = 40.0;
    currentRow = currentRow + 1;

    % Numeric entry for Y window size
    lbl_outlierWindowY = uilabel(gl_opts, 'Text', "Y outlier window size");
    lbl_outlierWindowY.Layout.Row = currentRow;
    lbl_outlierWindowY.Layout.Column = 1; 
    field_outlierWindowY = uieditfield(gl_opts, 'numeric', 'Value', options.outlierWindowY);
    field_outlierWindowY.Limits = [1, Inf];
    field_outlierWindowY.RoundFractionalValues = 'on';
    field_outlierWindowY.Layout.Row = currentRow;
    field_outlierWindowY.Layout.Column = 2;
    field_outlierWindowY.ValueChangedFcn = @outlierWindowYCallback;
    currentRow = currentRow + 1;

    % Table for cutoff bands
    gl_cutoff = uigridlayout(gl_opts, [3, 2]);
    gl_cutoff.Layout.Row = currentRow;
    gl_cutoff.Layout.Column = [1, 2];
    gl_cutoff.RowHeight = {defaultHeight, 'fit', defaultHeight};
    rowHeights{currentRow} = 'fit';
    
    % Make label
    lbl_cutoff = uilabel(gl_cutoff, 'Text', 'Remove data between specified x bounds.');
    lbl_cutoff.Layout.Row = 1;
    lbl_cutoff.Layout.Column = [1, 2];
    
    % Make table
    tbl_cutoff = uitable(gl_cutoff, 'Data', options.cutoffX, 'ColumnName', {'Min X', 'Max X'});
    tbl_cutoff.Layout.Column = [1, 2];
    tbl_cutoff.Layout.Row = 2;
    tbl_cutoff.ColumnEditable = [true, true];
    tbl_cutoff.CellEditCallback = @cutoffXTableCallback;
    
    % Make buttons for adding and removing rows
    btn_addRow = uibutton(gl_cutoff, 'Text', 'Add row');
    btn_addRow.Layout.Row = 3;
    btn_addRow.Layout.Column = 1;
    btn_addRow.ButtonPushedFcn = @addRowCallback;
    btn_rmvRow = uibutton(gl_cutoff, 'Text', 'Remove row');
    btn_rmvRow.Layout.Row = 3;
    btn_rmvRow.Layout.Column = 2;
    btn_rmvRow.ButtonPushedFcn = @removeRowCallback;
    rowHeights{currentRow} = 'fit';
    
    currentRow = currentRow + 1;

    % Row for buttons to clear kept and removed points
    gl_clear = uigridlayout(gl_opts, [1, 2]);
    gl_clear.Layout.Row = currentRow;
    gl_clear.Layout.Column = [1, 2];
    gl_clear.ColumnWidth = {'1x', '1x'};
    rowHeights{currentRow} = 'fit';
    btn_clearKeptPts = uibutton(gl_clear, 'Text', 'Clear kept points');
    btn_clearKeptPts.Layout.Row = 1;
    btn_clearKeptPts.Layout.Column = 1;
    btn_clearKeptPts.ButtonPushedFcn = @clearKeptPtsCallback;
    btn_clearRmvdPts = uibutton(gl_clear, 'Text', 'Clear removed points');
    btn_clearRmvdPts.Layout.Row = 1;
    btn_clearRmvdPts.Layout.Column = 2;
    btn_clearRmvdPts.ButtonPushedFcn = @clearRmvdPtsCallback;

    currentRow = currentRow + 1;

    % Penultimate row has button to apply options
    btn_apply = uibutton(gl_opts, "Text", "Apply options");
    btn_apply.Layout.Row = currentRow;
    btn_apply.Layout.Column = [1,2];
    btn_apply.ButtonPushedFcn = @applyOptionsCallback;
    currentRow = currentRow + 1;

    % Final row has button to export data
    btn_export = uibutton(gl_opts, "Text", "Export");
    btn_export.Layout.Row = currentRow;
    btn_export.Layout.Column = [1,2];
    btn_export.ButtonPushedFcn = @exportData;

    % Apply row heihgts
    gl_opts.RowHeight = rowHeights;
    
    % Generate initial plots
    updatePlots(options);

    %===================================================================
    %  End main function body. Next section contains callbacks and
    %  auxilliary functions.
    %===================================================================

    function exportData(~, ~)
    % Export data to one of several data formats
    % Raw data is exported with cleaned data

        % Query the user for an file of one of the listed file types
        [file, path] = uiputfile(...
            {'*.hdf5'; '*.mat'; '*.csv'; '*.dat'; '*.txt'; '*.xlsx'},...
            'Save data', outputPath...
        );

        % If a valid file is selected, then continue. Otherwise, exit.
        if ~isequal(file,0) || ~isequal(path,0)

            % Get the file extension and name
            [~, name, ext] = fileparts(file);
            fullpath = strcat(path, file);
        
            % Decide what to do for each extension
            switch ext
                case '.hdf5'
                    % For hdf5, we write the raw data and the cleaned data
                    % (both with metadata) to the file path
                    if exist(fullpath, 'file') == 2
                      delete(fullpath);
                    end
                    % Write raw data to file with metadata
                    rawDataSet = "/RawData";
                    h5create(fullpath, rawDataSet, [size(rawData, 2), size(rawData, 1)]);
                    h5write(fullpath, rawDataSet, rawData');
                    h5writeatt(fullpath, rawDataSet, 'headers', options.headers);
                    
                    % Write clean data to file
                    cleanDataSet = "/CleanedData";
                    h5create(fullpath, cleanDataSet, [size(cleanData, 2) size(cleanData, 1)]);
                    h5write(fullpath, cleanDataSet, cleanData');
                    h5writeatt(fullpath, cleanDataSet, 'headers', options.headers);

                    % Write metadata, including options the data cleaning
                    % was performed with, to aid reproducability
                    h5writeatt(fullpath, cleanDataSet, 'sort', tochar(options.sort));
                    h5writeatt(fullpath, cleanDataSet, 'removeOutliersX', tochar(options.removeOutliersX));
                    h5writeatt(fullpath, cleanDataSet, 'outlierMethodX', options.outlierMethodX);
                    h5writeatt(fullpath, cleanDataSet, 'outlierThresholdX', options.outlierThresholdX);
                    h5writeatt(fullpath, cleanDataSet, 'outlierWindowX', int64(options.outlierWindowX));
                    h5writeatt(fullpath, cleanDataSet, 'removeOutliersY', tochar(options.removeOutliersY));
                    h5writeatt(fullpath, cleanDataSet, 'outlierMethodY', options.outlierMethodY);
                    h5writeatt(fullpath, cleanDataSet, 'outlierThresholdY', options.outlierThresholdY);
                    h5writeatt(fullpath, cleanDataSet, 'outlierWindowY', int64(options.outlierWindowY));
                    h5writeatt(fullpath, cleanDataSet, 'indicesToKeep', int64(options.indicesToKeep));
                    h5writeatt(fullpath, cleanDataSet, 'indicesToRemove', int64(options.indicesToRemove));
                    h5writeatt(fullpath, cleanDataSet, 'cutoffX', options.cutoffX);
                    h5writeatt(fullpath, cleanDataSet, 'keptIndices', int64(keptInds));
                case '.mat'
                    % For .mat files, we used MATLAB's buit-in data saving
                    % Save the clean data, the raw data, and the options
                    % structure.
                    save(fullpath, "cleanData", "rawData", "options")
                otherwise
                    % If we're just writing ASCII data, we just save the
                    % cleaned data.
                    writematrix(cleanData, fullpath);
            end

            % Write cleaning options to file for reproducability if not
            % .hdf5 or .mat
            if ~isequal(ext, ".hdf5") && ~isequal(ext, ".mat")
                metadata_file = strcat(path, name, ext, ".json");
                options.keptIndices = keptInds;
                json_data = jsonencode(options, "PrettyPrint", true);
                fid = fopen(metadata_file, 'w');
                fprintf(fid, json_data);
                fclose(fid);
            end
        end
    end

    function [x2, y2, keptInds, rmvdInds] = clean(x, y, options)
        % Main cleaning function
        % Initialize results
        keptInds = (1:length(x))';
        x2 = x;
        y2 = y;
    
        % Find and remove NaNs/missing values
        [~, I] = rmmissing([x2 y2]);
        x2 = x2(~I);
        y2 = y2(~I);
        keptInds = keptInds(~I);
       
        % Remove X outliers.    
        if(options.removeOutliersX)
            dx = diff(x);
    
            switch options.outlierMethodX
                case {"mean", 'mean', "median", 'median'} 
                    [~, I] = rmoutliers(dx, options.outlierMethodX, ...
                    "ThresholdFactor", options.outlierThresholdX);
                otherwise
                    [~, I] = rmoutliers(dx, options.outlierMethodX, ...
                    options.outlierWindowX, ...
                    "ThresholdFactor", options.outlierThresholdX);
            end
            
            I = [true; ~I];
            x2 = x(I);
            y2 = y(I);
            keptInds = keptInds(I);
        end
    
        % Next, we sort the remaining data
        if (options.sort)
            [x2, I] = sort(x2);
            y2 = y2(I);
            keptInds = keptInds(I);
        end
    
        % Remove all data inside select x bands
        nbands = size(options.cutoffX, 1);
        for ii = 1:nbands
            x_low = options.cutoffX(ii, 1);
            x_high = options.cutoffX(ii, 2);
            mask = or(x2 < x_low, x2 > x_high);
            x2 = x2(mask);
            y2 = y2(mask);
            keptInds = keptInds(mask);
        end
    
        % Next, we remove Y outliers
        if(options.removeOutliersY)
            switch options.outlierMethodY
                case {"mean", 'mean', "median", 'median'} 
                    [~, I] = rmoutliers(y2, options.outlierMethodY, ...
                    "ThresholdFactor", options.outlierThresholdY);
                    keptInds = keptInds(~I);
                otherwise
                    [~, I] = rmoutliers(y2, options.outlierMethodY, ...
                    options.outlierWindowY, ...
                    "ThresholdFactor", options.outlierThresholdY);
                    keptInds = keptInds(~I);
            end
        end

        % Apply indices to keep and indices to remove
        keptInds = setdiff(unique([keptInds; options.indicesToKeep']), options.indicesToRemove);
        x2 = xRaw(keptInds);
        y2 = yRaw(keptInds);

        % Sort again if needed.
        if options.sort
            [x2, I] = sort(x2);
            y2 = y2(I);
            keptInds = keptInds(I);
        end

        % Find indices of removed points
        rmvdInds = setdiff(1:length(x), keptInds)';
    end

    function rawDataClickCallback(src, event)
        % Callback for when the plots are clicked
        xInt = event.IntersectionPoint(1);
        yInt = event.IntersectionPoint(2);
         
        [k, ~] = dsearchn(rawData, [xInt yInt]);
         
        % If k is already in the list of points to keep, remove it. 
        % Otherwise, add it. 
        ind_keep = 0;
        for ii = 1:length(options.indicesToKeep)
            if options.indicesToKeep(ii) == k
                ind_keep = ii;
                break
            end
        end

        if ind_keep == 0
            options.indicesToKeep(end+1) = k;
        else
            options.indicesToKeep = [...
                options.indicesToKeep(1:ind_keep-1) ... 
                options.indicesToKeep(ind_keep+1:end) ...
            ];
        end

        % If k is already in the list of points to remove, remove it. 
        % Otherwise, add it. 
        ind_remove = 0;
        for ii = 1:length(options.indicesToRemove)
            if options.indicesToRemove(ii) == k
                ind_remove = ii;
                break
            end
        end

        if ind_remove == 0 && ind_keep ~= 0
            options.indicesToRemove(end+1) = k;
        else
            options.indicesToRemove = [...
                options.indicesToRemove(1:ind_remove-1) ... 
                options.indicesToRemove(ind_remove+1:end) ...
            ];
        end

        % Make sure indices to keep and remove are unique and sorted
        options.indicesToKeep = unique(options.indicesToKeep);
        options.indicessToRemove = unique(options.indicesToRemove);

        % Apply the new options and refresh plots
        applyOptionsCallback(src, event);
    end

    function [x_sorted, y_sorted] = sortData(x, y)
        % Sorting helper function
        [x_sorted, I] = sort(x);
        y_sorted = y(I);
    end

    function sortDataCallback(~, event)
        % Callback to turn sorting on and off
        options.sort = event.Value;
    end

    function applyOptionsCallback(~, ~)
        % Callback for applying options. Also updates cleanData and
        % refreshes the plots.
        [xClean, yClean, keptInds, removedInds] = clean(xRaw, yRaw, options);
        cleanData = [xClean yClean];
        updatePlots(options)
    end

    function updatePlots(options)
        % Callback for updating plots. 
        if ~interactive
            return
        end
        xs = xRaw;
        ys = yRaw;
        if options.sort
            xs = xSorted;
            ys = ySorted;
        end
        cla(ax1)

        % Color excluded regions
        min_y = min(ys);
        max_y = max(ys);
        min_x = min(xs);
        max_x = max(xs);
        legend_objects = [];

        x_cutoff = xs;
        y_cutoff = ys;
        numBands = size(options.cutoffX, 1);
        for ii = 1:numBands
            x0 = options.cutoffX(ii, 1);
            x1 = options.cutoffX(ii, 2);
            x_fill = [x0, x1, x1, x0];

            ymax = 2 * max(abs(min_y), max_y);
            ymin = -2 * max(abs(min_y), max_y);
            y_fill = [ymin, ymin, ymax, ymax];
            color = [1, 0, 0];
            a = fill(ax1, x_fill, y_fill, color, 'FaceAlpha', 0.25);
            if ii == 1
                legend_objects(end+1) = a;
                hold(ax1, 'on');
            end
            % Remove points in cutoff region
            mask = or(x_cutoff < x0, x_cutoff> x1);
            x_cutoff = x_cutoff(mask);
            y_cutoff = y_cutoff(mask);
        end
        
        line = plot(ax1, xs, ys, 'k', 'LineWidth', 2);
        if numBands == 0
            hold(ax1, 'on');
        end
        line.ButtonDownFcn = @rawDataClickCallback;
        sc_excluded = scatter(ax1, xRaw(removedInds), yRaw(removedInds), 60, 'r', 'x', 'LineWidth', 2);
        sc_excluded.ButtonDownFcn = @rawDataClickCallback;

        sc_kept = scatter(ax1, xRaw(options.indicesToKeep), yRaw(options.indicesToKeep), 60, 'g', 'o', 'LineWidth', 2);
        sc_kept.ButtonDownFcn = @rawDataClickCallback;

        legend_objects(end+1) = sc_excluded;
        legend_objects(end+1) = sc_kept;

        pad = 0.1 * (max_y - min_y);
        xlim(ax1, [min_x, max_x])
        ylim(ax1, [min_y-pad, max_y+pad])
        xlabel(ax1, options.headers{1})
        ylabel(ax1, options.headers{2})
        title(ax1, "Raw data")
        box(ax1, "on")
        hold(ax1, 'off')
        
        if numBands > 0
            legend(ax1, legend_objects, {'Cutoff bands', 'Removed points', 'Kept points'}, 'Location', 'northoutside');
        else
            legend(ax1, legend_objects, {'Removed points', 'Kept points'}, 'Location', 'northoutside');
        end

        cla(ax2)
        sc_raw = scatter(ax2, x_cutoff, y_cutoff, 30, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
        sc_raw.ButtonDownFcn = @rawDataClickCallback;
        hold(ax2, 'on');
        line_raw = plot(ax2, x_cutoff, y_cutoff, 'r--', 'Color', [1.0, 0.0, 0.0, 0.5]);
        line_raw.ButtonDownFcn = @rawDataClickCallback;
        line2 = plot(ax2, xClean, yClean, 'k', 'LineWidth', 2);
        line2.ButtonDownFcn = @rawDataClickCallback;
        xlim(ax2, [min_x, max_x])
        min_yc = min(yClean);
        max_yc = max(yClean);
        pad = 0.1 * (max_yc - min_yc);
        ylim(ax2, [min_yc - pad, max_yc + pad])
        title(ax2, "Cleaned data")
        xlabel(ax2, options.headers{1})
        ylabel(ax2, options.headers{2})
        box(ax2, "on")
        hold(ax2, 'off')
    end

    function outlierXCheckboxCallback(~, event)
        % Callback for enabling/disabling x outlier removal
        options.removeOutliersX = event.Value;
    end

    function outlierThresholdXCallback(~, event)
        % Callback for setting the x outlier threshold
        options.outlierThresholdX = max(eps, event.Value);
    end

    function outlierYCheckboxCallback(~, event)
        % Callback for enabling/disabling y outlier removal
        options.removeOutliersY = event.Value;
    end

    function outlierThresholdYCallback(~, event)
        % Callback for updating the Y outlier threshold
        options.outlierThresholdY = max(eps, event.Value);
    end

    function outlierWindowYCallback(~, event)
        % Callback for updating the Y outlier window size
        options.outlierWindowY = event.Value;
    end

    function cutoffXTableCallback(~, ~)
        % Callback for refreshing cutoff table.
        options.cutoffX = tbl_cutoff.Data;
    end

    function addRowCallback(~, ~)
        % Callback for adding a row to cutoff table.
        tbl_cutoff.Data(end+1, :) = [-Inf, -Inf];
    end

    function removeRowCallback(src, event)
        % Callback for removing row from cutoff table.
        sz = size(tbl_cutoff.Data);
        if sz(1) > 0
            tbl_cutoff.Data = tbl_cutoff.Data(1:end-1, :);
        else
            tbl_cutoff.Data = [];
        end
        cutoffXTableCallback();
        applyOptionsCallback(src, event);
    end

    function clearKeptPtsCallback(src, event)
        % Callback for clearing selected points to keep.
        options.indicesToKeep = [];
        applyOptionsCallback(src, event);
    end

    function clearRmvdPtsCallback(src, event)
        % Callback for clearing selected points to remove.
        options.indicesToRemove = [];
        applyOptionsCallback(src, event);
    end

    function str = tochar(bool)
        % Convert a logical to a string 'true' or 'false'
        % Used in HDF5 output since HDF5 doesn't have logicals.
        if bool
            str = 'true';
        else
            str = 'false';
        end
    end
end
