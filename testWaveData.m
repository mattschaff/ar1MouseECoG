    %% Simulation Data: simple planer model (one wave with noise)
    wave_array = struct();
    for i=1:10
        wave_array(i).type = 'plane'; 
        wave_array(i).y_center = ones(1,5000);
        wave_array(i).x_center = ones(1,5000);
        wave_array(i).theta = ones(1,5000);
        wave_array(i).temp_freq = ones(1,5000);
        wave_array(i).spatial_freq = ones(1,5000);
        wave_array(i).amplitude = ones(1,5000);
        wave_array(i).timesteps = [1:5000];
    end
    disp(wave_array);
    
    