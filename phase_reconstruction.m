function phase_recon = phase_reconstruction(phase_matrix, kc, technique)
    phase_low = phase_matrix(1:kc,:);
    phase_high = phase_matrix(kc+1:end,:);
    phase_recon = phase_low;
    %Phase reconstruction technique lifting the narrowband into higherband.
    if technique == 1; %strcmp(technique,'Lift')
        while size(phase_matrix,1) - size(phase_recon,1) > 0
            delta = size(phase_matrix,1) - size(phase_recon,1)
            if delta >= kc
                phase_recon = [phase_recon; phase_low(1:end,:)];
            else 
                phase_recon = [phase_recon; phase_low(1:delta,:)];
            end
        end
    elseif technique == 2; %strcmp(technique,'Mirror')
        %Phase reconstruction technique mirroring the narrowband into highband around cutoff frequency.
        while size(phase_matrix,1) - size(phase_recon,1) > 0
            delta = size(phase_matrix,1) - size(phase_recon,1);
            if delta >= kc
                phase_recon = [phase_recon; flip(phase_low,1)];
            else 
                phase_recon = [phase_recon; flip(phase_low(1:delta,:),1)];
            end
        end
    elseif technique == 3; %strcmp(technique, 'Xcorr')
        %Phase reconstruction technique lifting the narrowband into
        %highband with appropriate amount of lag using cross-correlation.
        phi_array = [];
        for i = 1:size(phase_recon,2)
            [Rmm , lags] = xcorr(phase_low(:,i), phase_high(:,i));
            [peak, index] = max(Rmm);
            delta = lags(index) + 1;
                if delta > 0
                    phi_recon = [phase_low(:,i); phase_low(delta:end,i)];
                else
                    phi_recon = [phase_low(:,i); phase_low(1:end+delta,i)];
                end
            phi_recon = [phi_recon; flip(phi_recon,1)];
            phi_recon = phi_recon(1:size(phase_matrix,1));
            phi_recon = [phi_recon; zeros(size(phase_matrix,1)-size(phi_recon,1),1)];
            phi_array = [phi_array, phi_recon];
        end
        phase_recon = phi_array;
    else
        disp('Please provide correct phase reconstruction technique');
    end
end