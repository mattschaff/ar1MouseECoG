%% Set X,Y Arrays
    
    x = [-100:0.1:100];
    [X,Y] = meshgrid(x,x);

%% Planar Wave Phase Distribution

    %vars
    theta = pi/4;
    spatial_freq = 0.1;

    %phase distribution
    phase_planar = spatial_freq*(-cos(theta)*X + sin(theta)*Y);
    figure(1);
    imagesc(phase_planar);
    colorbar;
    title('Planar Wave Phase Distribution');

    %movie
    movie_length = 100; % units = frames
    freq = 2*(2*pi/movie_length);
    figure(2);
    for t=1:movie_length
        imagesc(cos(t*freq + phase_planar));
        title(['Time: ' num2str(t) '; Movie Length: ' num2str(movie_length)]);
        pause(0.001);
    end

%% Target Wave Phase Distribution

    %vars
    x_center = 0;
    y_center = 0;
    spatial_freq = 0.1;

    %phase distribution
    phase_target = -spatial_freq*sqrt((X-x_center).^2 + (Y-y_center).^2);
    figure(3);
    imagesc(phase_target);
    colorbar;
    title('Target Wave Phase Distribution');

    %movie
    movie_length = 100; % units = frames
    freq = 2*(2*pi/movie_length);
    figure(4);
    for t=1:movie_length
        imagesc(cos(t*freq + phase_target));
        title(['Time: ' num2str(t) '; Movie Length: ' num2str(movie_length)]);
        colorbar;
        pause(0.001);
    end

%% Rotational Wave Phase Distribution

    %vars
    x_center = 0;
    y_center = 0;
    spatial_freq = 10;
    
    %phase distribution
    phase_rotational = spatial_freq*atan2(Y-y_center, X-x_center);
    figure(5);
    imagesc(phase_rotational);
    colorbar;
    title('Rotational Wave Phase Distribution');
    
    %movie
    movie_length = 100; % units = frames
    freq = 2*(2*pi/movie_length);
    figure(6);
    for t=1:movie_length
        imagesc(cos(t*freq + phase_rotational));
        title(['Time: ' num2str(t) '; Movie Length: ' num2str(movie_length)]);
        colorbar;
        pause(0.001);
    end
