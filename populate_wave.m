function [data] = populate_wave(wave_data, x_grid, y_grid, times)
    data = zeros(size(x_grid,1), size(y_grid,2), length(times));
    
    X = x_grid;
    Y = y_grid;
    
    timesteps = wave_data.timesteps;
    
    switch wave_data.type
        case 'plane'
            for i = 1:length(timesteps)
                %vars
                theta = pi/4;
                spatial_freq = wave_data.spatial_freq(i);
                A = wave_data.amplitude(i);
                freq = wave_data.temp_freq(i);

                %phase distribution
                phase = spatial_freq*(-cos(theta)*X + sin(theta)*Y);

                %wave distribution
                data(:,:,timesteps(i)) = A*cos(freq*times(i) + phase);
            end
            
        case 'rotational'
            for i = 1:length(timesteps)
                %vars
                x_center = wave_data.x_center(i);
                y_center = wave_data.y_center(i);
                spatial_freq = wave_data.spatial_freq(i);
                A = wave_data.amplitude(i);
                freq = wave_data.temp_freq(i);

                %phase distribution
                phase = spatial_freq*atan2(x_grid-y_center, y_grid-x_center);

                %wave distribution
                data(:,:,timesteps(i)) = A*cos(freq*times(i) + phase);
            end

        case 'target'
            for i = 1:length(timesteps)
                %vars
                x_center = wave_data.x_center(i);
                y_center = wave_data.y_center(i);
                spatial_freq = wave_data.spatial_freq(i);
                A = wave_data.amplitude(i);
                freq = wave_data.temp_freq(i);

                %phase distribution
                phase = spatial_freq*atan2(Y-y_center, X-x_center);

                %wave distribution
                data(:,:,timesteps(i)) = A*cos(freq*times(i) + phase);
            end

        end
end