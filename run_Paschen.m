function [out] = run_Paschen(in)
    arguments
        in.T = 300;
        in.gamma = 0.01;
        in.fileList = {};
        in.fit_plot = false;
        in.pd_mark_val = [];
        in.pd_mark_label = {};
    end

    lnwdth = 2;

    % Optional colors (auto-repeat if fewer than files)
    colors = lines(numel(in.fileList));
    
    data = struct();
    paschen_curve = struct();
    
    %%  Load and process all files 
    for i = 1:numel(in.fileList)
        data(i).name = in.fileList{i};
        data(i).shortName = regexprep(in.fileList{i}, '\.dat$', ''); % for legend
        data(i).out = read_Bolsig(in.fileList{i});
        data(i).alpha = data(i).out.Transp.A18_Townsend_ioniz_coef_alphaN_m2;
        
        % Some files may not have attachment coefficients
        if isfield(data(i).out.Transp, 'A19_Townsend_attach_coef_etaN_m2')
            data(i).eta = data(i).out.Transp.A19_Townsend_attach_coef_etaN_m2;
        else
            data(i).eta = zeros(size(data(i).alpha));
        end
        
        % Effective Townsend coefficient (α - η)
        data(i).eff = data(i).alpha - data(i).eta;
        data(i).EN = data(i).out.Transp.E_N_Td;
        data(i).meanE = data(i).out.Transp.A1_Mean_energy_eV;
    end
    
    %%  Find breakdown point (E/N where α = η) 
    breakdown = struct();
    
    for i = 1:numel(data)
        % Compute difference between alpha and eta
        diffAE = data(i).alpha - data(i).eta;
    
        % Find sign changes (where α - η crosses zero)
        signChangeIdx = find(diffAE(1:end-1).*diffAE(2:end) < 0, 1, 'first');
        
        if ~isempty(signChangeIdx)
            % Linear interpolation between the two points around the crossing
            x1 = data(i).EN(signChangeIdx);
            x2 = data(i).EN(signChangeIdx+1);
            y1 = diffAE(signChangeIdx);
            y2 = diffAE(signChangeIdx+1);
            EN_eq = x1 - y1*(x2-x1)/(y2-y1); % interpolated E/N where α = η
    
            % Store results
            breakdown(i).EN_eq = EN_eq;
            breakdown(i).alpha_eq = interp1(data(i).EN, data(i).alpha, EN_eq);
            breakdown(i).eta_eq = breakdown(i).alpha_eq;
            breakdown(i).name = abbreviateName(data(i).shortName);
        else
            breakdown(i).EN_eq = NaN;
            breakdown(i).alpha_eq = NaN;
            breakdown(i).eta_eq = NaN;
            breakdown(i).name = abbreviateName(data(i).shortName);
        end
    end
    
    %%  Print summary table 
    fprintf('\n%-25s | %-10s | %-10s\n', 'Mixture', 'E/N (Td)', 'α=η (m²)');
    fprintf('%s\n', repmat('-',1,50));
    for i = 1:numel(breakdown)
        fprintf('%-25s | %10.3f | %10.3e\n', breakdown(i).name, ...
            breakdown(i).EN_eq, breakdown(i).alpha_eq);
    end
    
    %%  Plot α and η 
    figure; grid on
    for i = 1:numel(data)
        semilogy(data(i).EN, data(i).alpha, '--', 'Color', colors(i,:), 'LineWidth', lnwdth);
        if i==1, hold on; end
        semilogy(data(i).EN, data(i).eta, ':', 'Color', colors(i,:), 'LineWidth', lnwdth);
    end
    % Highlight α=η breakdown points
    for i = 1:numel(data)
        if ~isnan(breakdown(i).EN_eq)
            scatter(breakdown(i).EN_eq, breakdown(i).alpha_eq,"filled",'MarkerFaceColor',colors(i,:));
        end
    end
    legend('Location', 'best');
    xlabel('E/N (Td)');
    ylabel('Townsend coefficients (m^2)');
    hold off
    % xlim([0 35]); ylim([1E-30 inf]);
    
    % Build legend automatically (force char type for MATLAB legend)
    leg = cell(1, 2*numel(data));
    for i = 1:numel(data)
        short = abbreviateName(data(i).shortName);
        short = char(short); % ensure char type
        leg{2*i-1} = ['\alpha ' short];
        leg{2*i}   = ['\eta ' short];
    end
    legend(leg, 'Location', 'best');
    title('\alpha and \eta');
    
    %%
    %  Compute Paschen curves (Bolsig+ and classical fit) 
    k = 1.380649E-23; % Boltzmann constant [J/K]
    
    figure; grid on
    for i = 1:numel(data)
        %  Extract Bolsig+ data 
        EoverN = data(i).EN * 1E-21;   % Td → V·m^2 (E/N)
        alphaN = data(i).alpha - data(i).eta;          % α_eff / N
    
        %  Compute Paschen curve from effective alpha 
        pd = in.T*k*log(1+1/in.gamma)./alphaN;     % Eq.(5)
        Vb = EoverN./alphaN*log(1+1/in.gamma);     % Eq.(4)

        % trim negative values
        pd = pd(pd>0);
        Vb = Vb(Vb>0);

        loglog(pd, Vb, '-', 'Color', colors(i,:), 'LineWidth', lnwdth);
        if i==1, hold on; end
    
        paschen_curve(i).pd = pd;
        paschen_curve(i).Vb = Vb;
    
        %  Fit alpha/p = A exp(-B p/E) ==
    
        % α/p = (α/N) / (p/N) = αN / (kT)
        alpha_over_p = alphaN / (k*in.T);
    
        % (p/E) = (p/N) / (E/N) = (kT) / (E/N)
        p_over_E = (k*in.T) ./ EoverN;
    
        % Keep physically meaningful values
        valid = alpha_over_p > 0 & p_over_E > 0 & isfinite(alpha_over_p) & isfinite(p_over_E);
        x = p_over_E(valid);
        y = log(alpha_over_p(valid));
    
        if numel(x) > 5
            coeffs = polyfit(x, y, 1);
            lnA = coeffs(2);
            B   = -coeffs(1);
            A   = exp(lnA);
        else
            A = NaN; B = NaN;
            warning("Not enough valid points for fit in dataset %d", i);
        end
    
        fitted_params(i).A = A;
        fitted_params(i).B = B;
    
        %  Classical Paschen curve from fitted A,B =
        if ~isnan(A) && ~isnan(B)
            pd_fit  = logspace(-1, 3.5, 500);
            Vb_fit  = (B .* pd_fit) ./ ( log(A .* pd_fit) - log(log(1+1/in.gamma)) );

            i_neg = Vb_fit<0;

            pd_fit(i_neg) = [];
            Vb_fit(i_neg) = [];

            if in.fit_plot
                loglog(pd_fit, Vb_fit, '--', 'Color', colors(i,:), ...
                       'LineWidth', lnwdth);
            end
        end
    
    end
    xlabel('pd (Pa·m)');
    ylabel('Breakdown Voltage (V)');
    xlim([min(pd_fit)*0.8 max(pd_fit)*1.2]);
    ylim([min(Vb_fit)*0.8 max(Vb_fit)*1.2]);
    legend(cellfun(@abbreviateName, {data.shortName}, 'UniformOutput', false), 'Location', 'southeast');
    % title('Paschen Curves: Bolsig+ vs Exponential Fit');
    grid on
    if not(isempty(in.pd_mark_val))
        h = xline(in.pd_mark_val, '--r',in.pd_mark_label, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % ← key line
    end
        fontsize(16,"points")

    %%  Plot mean energy 
    figure; hold on; grid on
    for i = 1:numel(data)
        plot(data(i).EN, data(i).meanE, '-', 'Color', colors(i,:), 'LineWidth', lnwdth);
    end
    xlabel('E/N (Td)');
    ylabel('Mean electron energy (eV)');
    legend(cellfun(@abbreviateName, {data.shortName}, 'UniformOutput', false), 'Location', 'best');
    title('Mean Electron Energy');

    % output section
    out.pd = pd;
    out.Vb = Vb;
    out.A = A;
    out.B = B;
    out.pd_fit = pd_fit;
    out.Vb_fit = Vb_fit;

% Helper function 
function s = abbreviateName(fname)
    % Convert filename like "02_out_He_Air_01perc.dat" → "He-Air 0.1%"
    % and "05_out_He_Air_CO2_05perc.dat" -> "He-Air CO₂ 0.5%"
    s = regexprep(fname, '^(\d+_out_)', '');     % remove leading index
    s = regexprep(s, '\.dat$', '');              % remove extension
    s = strrep(s, '_', ' ');                     % underscores → spaces
    s = strtrim(s);

    % Detect the numeric value before "perc"
    tok = regexp(s, '(\d+)\s*perc', 'tokens', 'once');
    if ~isempty(tok)
        numStr = tok{1}; % e.g. '00', '01', '05', '1', '100'
        if startsWith(numStr, '0') && numel(numStr) > 1
            % leading zero: interpret as fractional percent
            denom = 10^(numel(numStr)-1);
            val = str2double(numStr) / denom;
        else
            val = str2double(numStr);
        end
        pctStr = sprintf('%g%%', val);
        s = regexprep(s, '\d+\s*perc', pctStr);
    end

    s = strtrim(s);
end

end