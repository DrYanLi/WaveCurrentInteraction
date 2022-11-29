function [px,py] = shear_data(sgn)
switch sgn
    case 1  % The one in the DispR manuscript, the shear profile 1-1
        p1 = 0;
        p2 = 0;
        p3 = 2.684;
        p4 = 8.784;
        p5 = 10.48;
        p6 = 5.367;
        p7 = 0.9884;
        px = [p1 p2 p3 p4 p5 p6 p7];
        py = [0 0 0 0 0 0 0];
    case 2  % The one in the DispR manuscript, the shear profile 1-2
        p1 = 0;
        p2 = 0;
        p3 = 0.1212;
        p4 = -0.0086;
        p5 = 3.041;
        p6 = 4.275;
        p7 = 1.098;
        px = [p1 p2 p3 p4 p5 p6 p7];
        py = [0 0 0 0 0 0 0];
    case 3 % The one in the DispR manuscript, the shear profile 1-3 
        p1 = 0;
        p2 = 0;
        p3 = 0.4921;
        p4 = 2.172;
        p5 = 3.811;
        p6 = 2.999;
        p7 = 1.509;
        px = [p1 p2 p3 p4 p5 p6 p7];
        py = [0 0 0 0 0 0 0];
    case 5 % 
        p1 = 0;
        p2 = 0;
        p3 = 0;
        p4 = 0;
        p5 = 0;
        p6 = 1;
        p7 = 0;
        px = [p1 p2 p3 p4 p5 p6 p7];
        py = [0 0 0 0 0 0 0];
end

end

