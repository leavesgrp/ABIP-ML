classdef (Abstract) model_interface
    methods (Abstract)
        A_times(obj, rhs)
        AT_times(obj, rhs)
        Q_times(obj, rhs)
        solve_lin_sys(obj, w, rhs, setting, init, warm_start, k)
        normalize_data(obj, w, setting)
        no_normalize_settings(obj, work, setting)
        individual_sparsity_ratio(obj)
        individual_data_work_settings(obj,work,setting,part)
 
    end
end