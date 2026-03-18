function extract_fit

defaultDirectory = cd;

selPath = uigetdir(defaultDirectory, 'Select folder');
if selPath == 0
    warning("No folder was selected. Analysis is cancelled");
    return
end
cd(selPath);

filesMAT = dir('*SpotSize*.mat');
fileNameMAT = {filesMAT(:).name};

if isempty(fileNameMAT)
    warning("No file was found. Analysis is cancelled");
    cd(defaultDirectory);
    return
end

table_xls = cell(1,1);
L = 1;

for i = 1:numel(fileNameMAT)

    close all

    currentfile = fileNameMAT{i};

    % Load into struct (robust to missing variables)
    S = load(currentfile);

    % If file doesn't contain expected content, skip gracefully
    if ~isfield(S, 'spot_record')
        warning("Skipping %s (missing 'spot_record')", currentfile);
        continue
    end

    % ---- Resolve / default all variables needed for the "main" table ----
    cum_arrival    = getfield_default(S, 'cum_arrival', []);
    cum_duration   = getfield_default(S, 'cum_duration', []);
    fit_association    = getfield_default(S, 'fit_association', []);
    fit_dissociation   = getfield_default(S, 'fit_dissociation', []);
    fit_dissociation_2exp = getfield_default(S, 'fit_dissociation_2exp', []);
    gof_analysis   = getfield_default(S, 'gof_analysis', struct());
    gof_2exp       = getfield_default(S, 'gof_2exp', struct());
    surf_coverage  = getfield_default(S, 'surf_coverage', struct());
    surf_cov_slope = getfield_default(S, 'surf_cov_slope', []);

    % Determine which fits exist
    if isempty(cum_arrival)
        fit_ass = [];
    else
        fit_ass = fit_association;
    end

    if isempty(cum_duration)
        fit_dis = [];
    else
        fit_dis = fit_dissociation;
    end

    if ~isempty(fit_dissociation_2exp)
        fit_2exp = fit_dissociation_2exp;
        gof_2exp_local = gof_2exp;
    else
        fit_2exp = [];
        gof_2exp_local = struct();
    end

    % Surface coverage slope (only if provided)
    if ~isempty(surf_cov_slope)
        if ~isstruct(surf_coverage) || isempty(surf_coverage)
            surf_coverage = struct();
        end
        surf_coverage.slope = surf_cov_slope;
    end

    [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp, ...
        gof_analysis, gof_2exp_local, surf_coverage, currentfile, L);
    L = L + 2;

    % ---- Fractions block (robust) ----
    if isfield(S, 'fraction') && ~isempty(S.fraction)
        fraction = S.fraction;

        intensity_fraction = getfield_default(S, 'intensity_fraction', []);
        if ~isempty(intensity_fraction)
            boundary = ["min", string(intensity_fraction), "max"];
        else
            boundary = ["min", "max"]; % fallback
        end

        if ~isempty(fraction) && isempty(fraction{1})
            fraction(1) = []; % compatibility
        end

        for j = 1:numel(fraction)

            Fj = fraction{j};
            if isempty(Fj), continue; end

            if ~isfield(Fj,'cum_arrival') || isempty(Fj.cum_arrival)
                fit_ass = [];
            else
                fit_ass = getfield_default(Fj,'fit_association',[]);
            end

            if ~isfield(Fj,'cum_duration') || isempty(Fj.cum_duration)
                fit_dis = [];
            else
                fit_dis = getfield_default(Fj,'fit_dissociation',[]);
            end

            if isfield(Fj,'fit_dissociation_2exp') && ~isempty(Fj.fit_dissociation_2exp)
                fit_2exp = Fj.fit_dissociation_2exp;
                gof_2exp_local = getfield_default(Fj,'gof_2exp',struct());
            else
                fit_2exp = [];
                gof_2exp_local = struct();
            end

            if numel(boundary) >= j+1
                title_sec = strcat("Intensity fraction between ", boundary(j), " and ", boundary(j+1));
            else
                title_sec = strcat("Intensity fraction #", string(j));
            end

            % surf coverage inside fraction
            surf_cov = getfield_default(Fj,'surf_coverage',struct());
            if isfield(Fj,'surf_cov_slope') && ~isempty(Fj.surf_cov_slope)
                if ~isstruct(surf_cov) || isempty(surf_cov)
                    surf_cov = struct();
                end
                surf_cov.slope = Fj.surf_cov_slope;
            end

            [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp, ...
                gof_analysis, gof_2exp_local, surf_cov, title_sec, L);
            L = L + 2;
        end
    end

    L = L + 2;
end

writecell(table_xls, fullfile(selPath, 'analysis cumulative.xlsx'));
cd(defaultDirectory);

end


function [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp, ...
    gof_analysis, gof_2exp, surf_coverage, currentfile, L)

% Safe getters for nested fields (return [] if missing)
getA = @(path) get_nested(gof_analysis, path);
getB = @(path) get_nested(gof_2exp,    path);

table_xls{L,1} = currentfile;
L = L + 1;

% ---- k_on ----
table_xls{L,1} = 'k_on';
if isempty(fit_ass)
    table_xls{L,2} = [];
else
    table_xls{L,2} = getfield_default(fit_ass,'p1',[]);
    table_xls{L,3} = getA({'association','arr_rate_95conf',1});
    table_xls{L,4} = getA({'association','arr_rate_95conf',2});
    table_xls{L,5} = 'R_sq';
    table_xls{L,6} = getA({'association','rsquare'});
end
L = L + 1;

% ---- k_off ----
table_xls{L,1} = 'k_off';
if isempty(fit_dis)
    table_xls{L,2} = [];
    L = L + 1;
    table_xls{L,1} = 'Irreversible fraction';
    table_xls{L,2} = [];
else
    table_xls{L,2} = getfield_default(fit_dis,'b',[]);
    table_xls{L,3} = getA({'dissociation','koff_95conf',1});
    table_xls{L,4} = getA({'dissociation','koff_95conf',2});
    table_xls{L,5} = 'R_sq';
    table_xls{L,6} = getA({'dissociation','rsquare'});

    L = L + 1;
    table_xls{L,1} = 'Irreversible fraction';

    a = getfield_default(fit_dis,'a',[]);
    c = getfield_default(fit_dis,'c',[]);
    if ~isempty(a) && ~isempty(c) && (a+c) ~= 0
        irr_fraction = c/(a+c)*100;
    else
        irr_fraction = [];
    end
    table_xls{L,2} = irr_fraction;
    table_xls{L,3} = getA({'dissociation','irr_frac_95conf',1});
    table_xls{L,4} = getA({'dissociation','irr_frac_95conf',2});

    % ---- Optional 2-exp block ----
    if ~isempty(fit_2exp)
        L = L + 1;
        table_xls{L,1} = 'k_off_1';
        b = getfield_default(fit_2exp,'b',[]);
        table_xls{L,2} = b;
        table_xls{L,3} = ternary(~isempty(b), b - getB({'dissociation','koff_95conf',1}), []);
        table_xls{L,4} = ternary(~isempty(b), b + getB({'dissociation','koff_95conf',1}), []);

        a2 = getfield_default(fit_2exp,'a',[]);
        c2 = getfield_default(fit_2exp,'c',[]);
        e2 = getfield_default(fit_2exp,'e',[]);
        rev_value = a2 + c2 + e2;

        if ~isempty(a2) && ~isempty(rev_value) && rev_value ~= 0
            fraction_1 = a2/rev_value*100;
        else
            fraction_1 = [];
        end
        table_xls{L,5} = 'Frac_1';
        table_xls{L,6} = fraction_1;

        frac1_conf = getB({'dissociation','frac1_95conf',1});
        if ~isempty(frac1_conf) && ~isempty(fraction_1) && ~isempty(a2) && a2 ~= 0
            perc_error = frac1_conf * fraction_1 / a2;
            table_xls{L,7} = fraction_1 - perc_error;
            table_xls{L,8} = fraction_1 + perc_error;
        else
            table_xls{L,7} = [];
            table_xls{L,8} = [];
        end
        table_xls{L,9}  = 'R_sq';
        table_xls{L,10} = getB({'dissociation','rsquare'});

        L = L + 1;
        table_xls{L,1} = 'k_off_2';
        d = getfield_default(fit_2exp,'d',[]);
        table_xls{L,2} = d;
        table_xls{L,3} = ternary(~isempty(d), d - getB({'dissociation','koff_95conf',2}), []);
        table_xls{L,4} = ternary(~isempty(d), d + getB({'dissociation','koff_95conf',2}), []);

        if ~isempty(c2) && ~isempty(rev_value) && rev_value ~= 0
            fraction_2 = c2/rev_value*100;
        else
            fraction_2 = [];
        end
        table_xls{L,5} = 'Frac_2';
        table_xls{L,6} = fraction_2;

        frac1_conf2 = getB({'dissociation','frac1_95conf',2});
        if ~isempty(frac1_conf2) && ~isempty(fraction_2) && ~isempty(c2) && c2 ~= 0
            perc_error = frac1_conf2 * fraction_2 / c2;
            table_xls{L,7} = fraction_2 - perc_error;
            table_xls{L,8} = fraction_2 + perc_error;
        else
            table_xls{L,7} = [];
            table_xls{L,8} = [];
        end
        table_xls{L,9}  = 'R_sq';
        table_xls{L,10} = getB({'dissociation','rsquare'});

        L = L + 1;
        table_xls{L,1} = 'Irreversible fraction';
        if ~isempty(e2) && ~isempty(rev_value) && rev_value ~= 0
            irr_fraction_2exp = e2/rev_value*100;
        else
            irr_fraction_2exp = [];
        end
        table_xls{L,2} = irr_fraction_2exp;
        table_xls{L,3} = getB({'dissociation','irr_frac_95conf',1});
        table_xls{L,4} = getB({'dissociation','irr_frac_95conf',2});
    end
end

L = L + 1;
table_xls{L,1} = 'Surface coverage';
if isstruct(surf_coverage) && ~isempty(surf_coverage)
    table_xls{L,2} = getfield_default(surf_coverage,'mean',[]);
    table_xls{L,3} = getfield_default(surf_coverage,'slope',[]);
else
    table_xls{L,2} = [];
    table_xls{L,3} = [];
end

end


% ------------ Helpers (robust field handling) ------------

function v = getfield_default(S, fname, defaultVal)
if isstruct(S) && isfield(S, fname)
    v = S.(fname);
else
    v = defaultVal;
end
end

function v = get_nested(S, path)
% path is a cell array like {'association','arr_rate_95conf',1}
v = [];
try
    cur = S;
    for k = 1:numel(path)
        key = path{k};
        if isnumeric(key)
            if isempty(cur) || numel(cur) < key
                v = [];
                return
            end
            cur = cur(key);
        else
            if ~isstruct(cur) || ~isfield(cur, key)
                v = [];
                return
            end
            cur = cur.(key);
        end
    end
    v = cur;
catch
    v = [];
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
