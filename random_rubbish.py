
# an attempt at sigma-cliping functions
def calcualte_noise(mz_vals, int_vals, mz_window_halfsize, 
                    theta_threshold =0.001, alpha = 1 ):
    
    medians = np.zeros(len(mz_vals))
    stds = np.zeros(len(mz_vals))
    # get window
    for i, mass in enumerate(mz_vals):
        lower_lim = mass - mz_window_halfsize
        upper_lim = mass + mz_window_halfsize
        
        # bitwise addition for final mask
        mask = np.bitwise_and(mz_vals>lower_lim, mz_vals<upper_lim)
        
        # cut the windows
        mz_window = mz_vals[mask]
        int_window = mz_vals[mask]
        
        # test here if mz_window actually contains n>1 elements, else everything is 0
        if len(mz_window) < 2:
            medians[i] = 0
            stds[i] = 0
            continue # continue with next iteration
        
        # get median and std.dev
        median_old = stat.median(int_window)
        std_dev_old = stat.stdev(int_window)
        
        # instance sigma_threshold
        theta = 1
        
        # sigma-clip in action
        while theta > theta_threshold:
            
            # get the interval 
            l_lim = median_old - alpha*std_dev_old
            u_lim = median_old + alpha*std_dev_old
        
            # bitwise addition for final mask
            mask = np.bitwise_and(int_window>l_lim, int_window<u_lim)
            
            # keep only values in interval
            mz_window = mz_window[mask]
            int_window = int_window[mask]
            
            median_new = stat.median(int_window)
            std_dev_new = stat.stdev(int_window)
            
            theta = (std_dev_old - std_dev_new)/std_dev_new
            median_old = median_new
            std_dev_old = std_dev_new
            
        # get the old_values
        medians[i] = median_old
        stds[i] = std_dev_old
        
    return medians, stds
