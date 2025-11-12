function result = ternary(condition, trueFunc, falseFunc)
    if condition
        result = trueFunc();
    else
        result = falseFunc();
    end
end