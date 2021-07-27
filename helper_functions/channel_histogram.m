function channel_histogram(data)
    
    x_dim = size(data,1);
    
    for i = 1:x_dim
        title(sprintf('Channel %s',i));
        plt = histogram(data(i,:));
        disp(i);
        path = sprintf('C:/Users/Jordan/Documents/cls_thesis/neuro_thesis/gc_hierarchies/jordan_results/data_exploration/pele_mean_channel_check/channel%d.png', i);
        saveas(plt, path);
    end
end