classdef cassegrain_spectrum
% Top-hat wavelength bands shared by the still and video simulators.
methods (Static)
    function short = short_rel(b)
        if ~(isscalar(b) && isreal(b) && b >= 0 && b < 2)
            error('cassegrain:spectrum', 'Bandwidth must be in [0, 2) for a positive lower wavelength');
        end
        short = 1 - b/2;
    end

    function [lambda_nm, b] = from_um(lower, upper)
        if ~(isscalar(lower) && isscalar(upper) && isreal(lower) && isreal(upper) && ...
             isfinite(lower) && isfinite(upper) && 0 < lower && lower <= upper)
            error('cassegrain:spectrum', 'Wavelengths must be finite, with 0 < lower <= upper (um)');
        end
        centre = lower + (upper - lower)/2;
        lambda_nm = centre * 1000;
        b = (upper - lower)/centre;
        if ~isfinite(lambda_nm) || b >= 2
            error('cassegrain:spectrum', 'Wavelength range is too large to represent');
        end
    end

    function [lower, upper] = edges_um(lambda_nm, b)
        centre = lambda_nm * 1e-3;
        lower = centre * cassegrain_spectrum.short_rel(b);
        upper = centre * (1 + b/2);
    end

    function nl = n_lambda(w, b)
        short = cassegrain_spectrum.short_rel(b);
        if b == 0, nl = 1; return; end
        % Limit relative bin width and blue-edge rim phase change per bin.
        bins = max(10*b/short, 2*abs(w)*b/short^2);
        if ~isfinite(bins) || bins > 128
            error('cassegrain:sampling', ...
                  'Band needs more than 128 wavelength samples; narrow the band or reduce defocus');
        end
        nl = max(floor(3 + 20*b + 0.7*b*abs(w)), ceil(bins));
    end
end
end
