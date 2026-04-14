function propname = fmam_target_canonicalize(propname)
%FMAM_TARGET_CANONICALIZE Normalize supported FMAM modulation target names.

    if isstring(propname)
        propname = char(propname);
    end
    if ~ischar(propname)
        error('Target property names must be character vectors or strings.')
    end

    switch lower(propname)
        case 'params'
            propname = 'params';
        case 'p_psi'
            propname = 'p_Psi';
        case 'q_psi'
            propname = 'q_Psi';
        case 'p_var'
            propname = 'p_var';
        case 'q_var'
            propname = 'q_var';
        case 'varphimax'
            propname = 'varPhiMax';
        case 'varphimin'
            propname = 'varPhiMin';
        case 'obsphimax'
            propname = 'obsPhiMax';
        case 'obsphimin'
            propname = 'obsPhiMin';
        case 'varamp'
            propname = 'varAmp';
        case 'varmax'
            propname = 'varMax';
        case 'varmin'
            propname = 'varMin';
        case 'obsamp'
            propname = 'obsAmp';
        case 'obsmax'
            propname = 'obsMax';
        case 'obsmin'
            propname = 'obsMin';
        case 'varphase'
            propname = 'varPhase';
        case 'obsphase'
            error('obsPhase targets are not supported by FMAM_ODE.')
        otherwise
            error('Unsupported modulation target property ''%s''.',propname)
    end
end
