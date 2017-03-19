function var = return_variance(experiment,model,p)
     var = sum((experiment-model).^2)/(length(experiment)-p);
end