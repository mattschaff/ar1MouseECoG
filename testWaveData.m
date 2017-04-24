    %% Simulation Data: simple planer model (one wave with noise)
    wave_array = struct();
    for i=1:10
        wave_array(i).type = 'target'; 
        wave_array(i).y_center = ones(1,5000);
        wave_array(i).x_center = ones(1,5000);
        wave_array(i).theta = ones(1,5000);
        wave_array(i).temp_freq = ones(1,5000)*0.1;
        wave_array(i).spatial_freq = ones(1,5000);
        wave_array(i).amplitude = ones(1,5000);
        wave_array(i).timesteps = [1:5000];
    end
    disp(wave_array);
    
    x = -1:0.1:1;
    [X, Y] = meshgrid(x, x);
    times = 1:500;
    
    data = populate_wave(wave_array(1), X, Y, times);
    
    figure(1);
    clf;
    
    for i = 1:size(data,3)
        imagesc(data(:,:,i));
        title(num2str(i));
        
        pause(0.01);
    end
    
    %% Simulated Planar Waves
    wave_array = struct();
    for i=1:2
        wave_array(i).type = 'plane'; 
        wave_array(i).y_center = ones(1,5000);
        wave_array(i).x_center = ones(1,5000);
        wave_array(i).theta = ones(1,5000);
        wave_array(i).temp_freq = ones(1,5000);
        wave_array(i).spatial_freq = ones(1,5000)*5;
        wave_array(i).amplitude = ones(1,5000);
        wave_array(i).timesteps = [1:5000];
    end
    wave_array(1).theta = ones(1,5000).*0;
    wave_array(2).theta = ones(1,5000).*pi;
    wave_array(2).temp_freq = ones(1,5000).*1.1;
    %initialize grid
    x = -1:0.1:1;
    y = 0;
    [X, Y] = meshgrid(x, y);
    times = (1:500)*.001;
    
    data1 = populate_wave(wave_array(1), X, Y, times);
    data2 = populate_wave(wave_array(2), X, Y, times);
    combined_data = data1 + data2;
    combined_data_rand = combined_data + normrnd(0,.3,size(data1));
%     figure;
%     subplot(1,3,1);
%     imagesc(squeeze(data1)');
%     subplot(1,3,2);
%     imagesc(squeeze(data2)');
%     subplot(1,3,3);
%     imagesc(squeeze(combined_data)');
    figure(1);
    clf;
    
    for i = 1:size(combined_data,3)
        imagesc(combined_data(:,:,i));
        title(num2str(i));
        caxis([-1,1]);
        pause(0.1);        
    end
    %% asdf
    figure(2);
    clf;
    hold on;
    %plot(squeeze(data2(1,1,:)));
    %plot(squeeze(data1(1,1,:)));
    plot(squeeze(combined_data(1,1,:)));
    plot(squeeze(combined_data_rand(1,1,:)));
    hold off;