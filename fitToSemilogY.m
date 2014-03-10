function y = fitToSemilogY(x_sample, y_sample, input)
    fit = polyfit(x_sample, log(y_sample), 1);
    y = exp(fit(2)).*exp(input.*fit(1));
end
