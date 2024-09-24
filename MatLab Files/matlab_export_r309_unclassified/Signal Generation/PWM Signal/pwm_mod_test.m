function outScale = pwm_mod_test(inFrac)
    % test function to scale from 25% to 75% of the input
    x = inFrac;
    x1 = 0;
    y1 = 0.10;
    x2 = 1;
    y2 = 0.90;
    b = 0.10;
    
    outScale = ((y2-y1)/(x2-x1))*(x-x1) + b;
    
end