function dxds = shape(s3, x, C2)

dxds = [ C2(s3);
        cos(x(1));
        sin(x(1))];
end

