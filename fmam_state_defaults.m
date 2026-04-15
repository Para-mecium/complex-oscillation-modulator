classdef fmam_state_defaults
    properties (Constant)
        assemblySampleCount = 500
        timeResampleCount = 20000
        phaseSampleCount = 500
    end

    methods (Static)
        function discretization = defaultDiscretization()
            discretization = struct( ...
                'assemblySampleCount', fmam_state_defaults.assemblySampleCount, ...
                'reconstruction', struct( ...
                    'timeResampleCount', fmam_state_defaults.timeResampleCount, ...
                    'phaseSampleCount', fmam_state_defaults.phaseSampleCount));
        end

        function discretization = normalizeDiscretization(discretization)
            defaults = fmam_state_defaults.defaultDiscretization();
            if nargin < 1 || isempty(discretization)
                discretization = defaults;
                return
            end
            if ~isstruct(discretization)
                error('fmam_state_defaults:InvalidDiscretization', ...
                    'discretization must be a struct.')
            end

            allowedFields = {'assemblySampleCount', 'reconstruction'};
            fieldNames = fieldnames(discretization);
            for i = 1:numel(fieldNames)
                if ~ismember(fieldNames{i}, allowedFields)
                    error('fmam_state_defaults:InvalidDiscretization', ...
                        'Unsupported discretization field ''%s''.', fieldNames{i})
                end
            end

            merged = defaults;
            for i = 1:numel(fieldNames)
                fieldName = fieldNames{i};
                if strcmp(fieldName, 'reconstruction')
                    continue
                end
                merged.(fieldName) = discretization.(fieldName);
            end
            if isfield(discretization, 'reconstruction')
                merged.reconstruction = fmam_state_defaults.normalizeReconstruction( ...
                    discretization.reconstruction, defaults.reconstruction);
            end

            validateattributes(merged.assemblySampleCount, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'finite'}, ...
                'fmam_state_defaults', 'discretization.assemblySampleCount');
            merged.assemblySampleCount = double(merged.assemblySampleCount);
            discretization = merged;
        end

        function reconstruction = normalizeReconstruction(reconstruction, defaults)
            if nargin < 2 || isempty(defaults)
                defaults = fmam_state_defaults.defaultDiscretization().reconstruction;
            end
            if nargin < 1 || isempty(reconstruction)
                reconstruction = defaults;
                return
            end
            if ~isstruct(reconstruction)
                error('fmam_state_defaults:InvalidReconstruction', ...
                    'reconstruction must be a struct.')
            end

            allowedFields = {'timeResampleCount', 'phaseSampleCount'};
            fieldNames = fieldnames(reconstruction);
            for i = 1:numel(fieldNames)
                if ~ismember(fieldNames{i}, allowedFields)
                    error('fmam_state_defaults:InvalidReconstruction', ...
                        'Unsupported reconstruction field ''%s''.', fieldNames{i})
                end
            end

            merged = defaults;
            for i = 1:numel(fieldNames)
                merged.(fieldNames{i}) = reconstruction.(fieldNames{i});
            end

            validateattributes(merged.timeResampleCount, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'finite'}, ...
                'fmam_state_defaults', 'reconstruction.timeResampleCount');
            validateattributes(merged.phaseSampleCount, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'finite'}, ...
                'fmam_state_defaults', 'reconstruction.phaseSampleCount');

            merged.timeResampleCount = double(merged.timeResampleCount);
            merged.phaseSampleCount = double(merged.phaseSampleCount);
            reconstruction = merged;
        end
    end
end
