obs = cell(1,1);
obs{1} = str2func('f1');

function output = f1(y)
    output = y(:,1) + y(:,2);
end