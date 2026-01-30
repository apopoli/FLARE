function Bolsig = read_Bolsig(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    Bolsig = struct();

    Bolsig.Rate_Coeffs = read_rate_block(fid, 'Rate coefficients (m3/s)', true);

    fseek(fid, 0, 'bof');  % rewind file
    try
        Bolsig.Inv_Rate_Coeffs = read_rate_block(fid, 'Inverse rate coefficients (m3/s)', true);
    catch
        warning('Inverse rate coefficients block not found.');
    end

    fseek(fid, 0, 'bof');  % rewind again
    try
        Bolsig.Transp = read_rate_block(fid, 'Transport coefficients', false);
    catch
        warning('Transport coefficients block not found.');
    end

    fclose(fid);
end

function block_data = read_rate_block(fid, block_title, has_energy_column)
    block_data = struct();
    found_block = false;
    header_lines = {};

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if contains(line, block_title)
            found_block = true;

            while true
                line = strtrim(fgetl(fid));
                if startsWith(line, 'R#') || contains(line, 'E/N (Td)')
                    column_headers_line = line;
                    break;
                end
                header_lines{end+1} = line;
            end
            break;
        end
    end

    if ~found_block
        error('Block "%s" not found.', block_title);
    end

    % Parse column headers
    raw_tokens = strsplit(strtrim(column_headers_line));
    col_names = {};
    i = 1;
    while i <= numel(raw_tokens)
        token = raw_tokens{i};
        if i < numel(raw_tokens) && strcmp(token, 'E/N') && strcmp(raw_tokens{i+1}, '(Td)')
            col_names{end+1} = 'E/N (Td)';
            i = i + 2;
        elseif has_energy_column && i < numel(raw_tokens) && strcmp(token, 'Energy') && strcmp(raw_tokens{i+1}, '(eV)')
            col_names{end+1} = 'Energy (eV)';
            i = i + 2;
        else
            col_names{end+1} = token;
            i = i + 1;
        end
    end
    n_cols = numel(col_names);

    % Sanitize headers
    field_names = {'E_N_Td'};
    for i = 1:numel(header_lines)
        sanitized = regexprep(strtrim(header_lines{i}), '\s+', '_');
        sanitized = regexprep(sanitized, '[^a-zA-Z0-9_]', '');
        field_names{end+1} = sanitized;
    end

    % Read data
    data = textscan(fid, repmat('%f', 1, n_cols));

    % Assign values
    block_data.E_N_Td = data{2};  % E/N (Td) is 2nd column
    if strcmp(block_title,'Transport coefficients')
        k = 3;
    else
        k = 4;
    end
    for i = k:n_cols  % Skip R#, E/N (Td), and optionally Energy (eV)
        if ~all(isnan(data{i}))
            if k==3
                field = field_names{i - 1};  % shift by 1: we already assigned E/N
            else
                field = field_names{i - 2};  % shift by 2: we already assigned E/N, Energy if present
            end
            block_data.(field) = data{i};
        end
    end
end
